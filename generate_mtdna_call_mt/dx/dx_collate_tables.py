# Pipeline for the merging of mtSwirl pipeline outputs on DNANexus.
# This pipeline also merges QC files from DNANexus. Note that these files were not
# produced for all samples, though it seems like the set of samples for which QC files
# are available corresponds to the set of files for which there was sufficient
# freemix contamination to warrant a second run of contamination estimation.

import hail as hl
import pandas as pd
import numpy as np
import argparse
import dxpy
import pyspark
import psutil, multiprocessing
import json
import sys, os, glob, re
import collections

from datetime import datetime
from dxpy.utils import file_load_utils
from dxpy.bindings.download_all_inputs import _parallel_file_download, _get_num_parallel_threads, _create_dirs, \
                                              _sequential_file_download
from functools import reduce, partial

FUSE_PREFIX = '/mnt/project/'
AUTOSOMES = ['chr' + str(x) for x in range(1, 23)]


def make_input_json(d, fl_out="sample.json"):
    f = open(fl_out, "w")
    json.dump(d, f)
    f.close()
    return fl_out


def make_keyed_df(key, lst, key_to_suffix):
    df = pd.DataFrame({key: lst})
    df['s'] = df[key].apply(os.path.basename)
    df['s'] = df['s'].str.replace(key_to_suffix[key]+'$',"")
    return df


def make_keyed_df_multi(key, lst, key_to_suffix):
    df = pd.DataFrame({key: lst})
    df['batch'] = df[key].apply(lambda x: os.path.basename(os.path.dirname(x)))
    df['batch'] = df['batch'].str.replace(key_to_suffix[key]+'$',"")
    return df


def reader1(args):
    idx, f = args
    df = pd.read_csv(f, index_col=0, header=None, sep='\t').transpose().assign(newcolthis=idx)
    return df


def reader2(args, filter_to=None):
    idx, f = args
    tab = pd.read_csv(f, index_col=None, header=0, sep='\t').assign(newcolthis=idx)
    if filter_to is not None:
        tab = tab[tab['chr'].isin(filter_to)]
    return tab


def custom_get_job_input_filenames(input_dict):
    """Extract list of files, returns a set of directories to create, and
    a set of files, with sources and destinations. The paths created are
    relative to the input directory.

    Note: we go through file names inside arrays, and create a
    separate subdirectory for each. This avoids clobbering files when
    duplicate filenames appear in an array.
    """

    files = collections.defaultdict(list)  # dictionary, with empty lists as default elements
    dirs = []  # directories to create under <idir>

    # Local function for adding a file to the list of files to be created
    # for example:
    #    iname == "seq1"
    #    subdir == "015"
    #    value == { "$dnanexus_link": {
    #       "project": "project-BKJfY1j0b06Z4y8PX8bQ094f",
    #       "id": "file-BKQGkgQ0b06xG5560GGQ001B"
    #    }
    # will create a record describing that the file should
    # be downloaded into seq1/015/<filename>
    def add_file(iname, subdir, value, descr):
        if not dxpy.is_dxlink(value):
            return
        if descr is None:
            handler = dxpy.get_handler(value)
            if not isinstance(handler, dxpy.DXFile):
                return
            thisnm = handler.name
            thisid = handler.id
        else:
            handler = descr
            thisnm = handler['name']
            thisid = handler['id']
        filename = file_load_utils.make_unix_filename(thisnm)
        print('ding2')
        trg_dir = iname
        if subdir is not None:
            trg_dir = os.path.join(trg_dir, subdir)
        files[iname].append({'trg_fname': os.path.join(trg_dir, filename),
                             'handler': handler,
                             'src_file_id': thisid})
        dirs.append(trg_dir)

    # An array of inputs, for a single key. A directory
    # will be created per array entry. For example, if the input key is
    # FOO, and the inputs are {A, B, C}.vcf then, the directory structure
    # will be:
    #   <idir>/FOO/00/A.vcf
    #   <idir>/FOO/01/B.vcf
    #   <idir>/FOO/02/C.vcf
    def add_file_array(input_name, link_tuples):
        link_tuples = list(link_tuples)
        num_files = len(link_tuples)
        if num_files == 0:
            return
        num_digits = len(str(num_files - 1))
        dirs.append(input_name)
        for i, (link, descr) in enumerate(link_tuples):
            subdir = str(i).zfill(num_digits)
            add_file(input_name, subdir, link, descr)

    for input_name, value in list(input_dict.items()):
        if isinstance(value, zip):
            # This is a file array
            add_file_array(input_name, value)
        else:
            add_file(input_name, None, value, None)

    ## create a dictionary of the all non-file elements
    rest_hash = {key: val for key, val in list(input_dict.items()) if key not in files}
    return dirs, files, rest_hash


def custom_download_all_inputs(input_links, suffix_key, exclude=None, parallel=False, max_threads=8, enforce_nonmissing=True):
    '''
    :param exclude: List of input variables that should not be downloaded.
    :type exclude: Array of strings
    :param parallel: Should we download multiple files in parallel? (default: False)
    :type filename: boolean
    :param max_threads: If parallel is True, how many threads should be used
        to download files? (default: 8)
    :type append: int
    :returns: dict of lists of strings where each key is the input variable
                and each list element is the full path to the file that has
                been downloaded.

    A custom version of the download_all_inputs function that does not assume an input JSON file.
    By convention, if an input parameter "FOO" has value

        {"$dnanexus_link": "file-xxxx"}

    and filename INPUT.TXT, then the linked file will be downloaded into the
    path:

        $HOME/in/FOO/INPUT.TXT

    If an input is an array of files, then all files will be placed into
    numbered subdirectories under a parent directory named for the
    input. For example, if the input key is FOO, and the inputs are {A, B,
    C}.vcf then, the directory structure will be:

        $HOME/in/FOO/0/A.vcf
                     1/B.vcf
                     2/C.vcf

    Zero padding is used to ensure argument order. For example, if there are
    12 input files {A, B, C, D, E, F, G, H, I, J, K, L}.txt, the directory
    structure will be:

        $HOME/in/FOO/00/A.vcf
                     ...
                     11/L.vcf
    '''

    # Input directory, where all inputs are downloaded
    idir = file_load_utils.get_input_dir()
    try:
        dirs, inputs, rest = custom_get_job_input_filenames(input_links)
    except IOError:
        msg = 'Error: Could not find the input json file: {0}.\n'.format('<<DICT INPUT>>')
        msg += '       This function should only be called from within a running job.'
        print(msg)
        raise

    print('Filenames obtained.')
    # Exclude directories
    # dirs contain all folders (e.g. $HOME/in/FOO) and their sub folders (e.g. $HOME/in/FOO/1, $HOME/in/FOO/2, etc.)
    # If the main folder is excluded, its sub-folder would also be excluded from dirs_to_create
    dirs_to_create = []
    for d in dirs:
        keep = True
        if (exclude is not None) and (d is not None):
            if (d.split('/')[0] in exclude):
                keep = False
        if keep:
            dirs_to_create.append(d)

    # Create the directory structure, in preparation for download.
    # Allows performing the download in parallel.
    _create_dirs(idir, dirs_to_create)
    print('Directories created.')

    # Remove excluded inputs
    if exclude:
        inputs = file_load_utils.filter_dict(inputs, exclude)

    # Convert to a flat list of elements to download
    to_download = []
    for ival_list in inputs.values():
        to_download.extend(ival_list)

    # Download the files
    if parallel:
        total_mem = psutil.virtual_memory().total >> 20  # Total RAM in MB
        num_cores = multiprocessing.cpu_count()
        max_num_parallel_downloads = _get_num_parallel_threads(max_threads, num_cores, total_mem)
        sys.stderr.write("Downloading files using {} threads".format(max_num_parallel_downloads))
        _parallel_file_download(to_download, idir, max_num_parallel_downloads)
    else:
        _sequential_file_download(to_download, idir)

    # output a pandas table with columns for all downloaded samples with sample IDs
    trimmed_inputs = [make_keyed_df(k, [idir + '/' + subdict['trg_fname'] for subdict in v], suffix_key) for k, v in inputs.items()]
    joint_table = reduce(lambda x, y: pd.merge(x, y, on = 's', how='outer'), trimmed_inputs)
    if enforce_nonmissing:
        # verify that all tables are non-missing
        tf_null = joint_table.isnull().any().any()
        if tf_null:
            raise ValueError('ERROR: all df values should be non-missing.')
    
    return joint_table


def import_and_cat_tables(directory_df, id_col, path_col, new_id_col, append_ids_and_t=False, max_threads=8, filter_by=None, enforce_nonmiss=False, filter_chr=None):
    ids = [x for _, x in directory_df[id_col].iteritems()]
    directories = [x for _, x in directory_df[path_col].iteritems()]
    p = multiprocessing.Pool(processes=max_threads)

    to_subset_to = filter_by if filter_by is not None else ids
    empty_df = pd.DataFrame({new_id_col:[]})

    if append_ids_and_t:
        df_from_each_file = p.map(reader1, [(idx, f) for idx, f in zip(ids, directories) if idx in to_subset_to])
    else:        
        df_from_each_file = p.map(partial(reader2, filter_to=filter_chr), [(idx, f) for idx, f in zip(ids, directories) if idx in to_subset_to])
    
    if len(df_from_each_file) == 0:
        concatenated_df = empty_df
    else:
        concatenated_df = pd.concat(df_from_each_file, ignore_index=True, axis=0)
        concatenated_df = concatenated_df.rename({'newcolthis': new_id_col}, axis=1)

    # ensure that all ids are found in the new df
    if enforce_nonmiss:
        new_ids = [x for _, x in concatenated_df[new_id_col].iteritems()]
        tf_found_new = all([x_old in new_ids for x_old in ids])
        tf_found_old = all([x_new in ids for x_new in new_ids])
        if not (tf_found_new and tf_found_old):
            raise ValueError('ERROR: the list of individuals in directory_df must be the same as that obtained from the read tables.')
    
    return concatenated_df.reset_index(drop=True)


def fuse_find_data_objects(folder, suffix, recursive=True):
    folders = folder.split(',')
    identified_objects = []
    for this_folder in folders:
        path_search = FUSE_PREFIX + this_folder + '**/*' + suffix
        identified_objects+=glob.glob(path_search, recursive=recursive)
    return identified_objects


def dx_find_data_objects(folder, suffix, recursive=True):
    folders = folder.split(',')
    identified_objects = []
    for this_folder in folders:
        gen_item = dxpy.find_data_objects(classname='file', name=suffix, name_mode='glob', describe=True, folder='/'+this_folder, recurse=recursive)
        identified_objects+=[FUSE_PREFIX + re.sub('^/', '', x['describe']['folder']) + '/' + x['describe']['name'] for x in gen_item]
    return identified_objects


def subset_to_data_objects(lst, suffix):
    return [x for x in lst if re.search(suffix, x)]


def produce_fuse_file_table(input_links, suffix_key, enforce_nonmissing=True, single_sample=False):
    # output a pandas table with columns for all downloaded samples with sample IDs
    fun_use = make_keyed_df if single_sample else make_keyed_df_multi
    trimmed_inputs = [fun_use(k, v, suffix_key) for k, v in input_links.items()]
    joint_table = reduce(lambda x, y: pd.merge(x, y, on = 'batch', how='outer'), trimmed_inputs)
    if enforce_nonmissing:
        # verify that all tables are non-missing
        tf_null = joint_table.isnull().any().any()
        if tf_null:
            raise ValueError('ERROR: all df values should be non-missing.')
    
    return joint_table


def run_describe(lst):
    """
    dxpy.describe() only allows up to 1000 elements; here we loop across the size of list
    """
    n = len(lst)
    n_split = 999
    split_lists = [lst[i * n_split:(i + 1) * n_split] for i in range((n + n_split - 1) // n_split )] 
    return [y for x in split_lists for y in dxpy.describe(x)]


def compute_nuc_coverage(df, column_name):
    df[column_name] = ((df.total_mapped_reads - df.singletons - df.mate_diff_chr - df.duplicates) * df.READ_LENGTH ) / df.genome_length
    return df


def main(pipeline_output_folder, vcf_suffix, coverage_suffix, mtstats_suffix, yield_suffix, idxstats_suffix, qc_stats_folder, qc_suffix,
         file_paths_table_output, per_sample_stats_output, dx_init, avoid_filtering_idxstats_chr, unified_prefix):

    # start SQL session
    my_database = dxpy.find_one_data_object(name=dx_init.lower())["id"]
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    hl.init(sc=sc, tmp_dir=f'dnax://{my_database}/tmp2/')
    hl._set_flags(no_whole_stage_codegen='1')

    # download mito pipeline data
    print(f'{datetime.now().strftime("%H:%M:%S")}: Finding all relevant data objects...')
    all_batch_files = dx_find_data_objects(pipeline_output_folder, unified_prefix + '*', True)
    data_dict = {'vcf': subset_to_data_objects(all_batch_files, vcf_suffix), 
                 'coverage': subset_to_data_objects(all_batch_files, coverage_suffix),
                 'stats': subset_to_data_objects(all_batch_files, mtstats_suffix),
                 'yield': subset_to_data_objects(all_batch_files, yield_suffix),
                 'idxstats': subset_to_data_objects(all_batch_files, idxstats_suffix)}
    
    # checks on dict
    print(f'{datetime.now().strftime("%H:%M:%S")}: Running checks...')
    batches = {k: [os.path.basename(os.path.dirname(x)) for x in v] for k,v in data_dict.items()}
    if not all([len(set(x)) == len(x) for _, x in batches.items()]):
        for k,v in batches.items():
            print(f'Duplicate batches in {k}:')
            multiple_found = {x: v.count(x) for x in set(v) if v.count(x) > 1}
            for item, idx in multiple_found.items():
                print(item + ' has ' + str(idx) + ' items.')
        raise ValueError('ERROR: there are duplicate batches (or multiple files per batch).')
    
    # obtain paths
    print(f'{datetime.now().strftime("%H:%M:%S")}: Obtaining paths...')
    downloaded_files = produce_fuse_file_table(data_dict, {'vcf':vcf_suffix, 'coverage':coverage_suffix, 'stats':mtstats_suffix, 'yield':yield_suffix, 'idxstats':idxstats_suffix})

    # checks on downloaded data
    print(f'{datetime.now().strftime("%H:%M:%S")}: Checking file paths...')
    per_row_un = downloaded_files.apply(lambda x: x.apply(lambda y: os.path.basename(os.path.dirname(y))).unique(), 1)
    per_row_un = per_row_un.apply(lambda x: [item for item in x if len(item) != 0])
    if not all(per_row_un.apply(len) == 1):
        raise ValueError('ERROR: all imported files must be ordered in the same way.')
    if not all(downloaded_files['batch'] == per_row_un.map(lambda x: x[0])):
        raise ValueError('ERROR: all imported files must have the same batch order as the specified batch name.')

    # read qc metrics
    print(f'{datetime.now().strftime("%H:%M:%S")}: Downloading QC metrics...')
    data_dict_qc = {'qc': dx_find_data_objects(qc_stats_folder, '*'+qc_suffix, True)}
    downloaded_qc_files = produce_fuse_file_table(data_dict_qc, {'qc':qc_suffix}, single_sample=True)

    # import stats and qc and merge all into table
    print(f'{datetime.now().strftime("%H:%M:%S")}: Importing all statistics...')
    print(f'{datetime.now().strftime("%H:%M:%S")}: Importing run statistics...')
    stats_table = import_and_cat_tables(downloaded_files, 'batch', 'stats', 'batch', enforce_nonmiss=False)
    print(f'{datetime.now().strftime("%H:%M:%S")}: Importing yield statistics...')
    yield_table = import_and_cat_tables(downloaded_files, 'batch', 'yield', 'batch', enforce_nonmiss=False)
    print(f'{datetime.now().strftime("%H:%M:%S")}: Importing idxstats statistics...')
    if avoid_filtering_idxstats_chr:
        idxstats_table = import_and_cat_tables(downloaded_files, 'batch', 'idxstats', 'batch', enforce_nonmiss=False)
    else:
        idxstats_table = import_and_cat_tables(downloaded_files, 'batch', 'idxstats', 'batch', enforce_nonmiss=False, filter_chr=AUTOSOMES)

    # produce munged idxstats table
    idxstats_summary = idxstats_table.groupby(idxstats_table['s']).agg(
        total_mapped_reads=pd.NamedAgg(column='mapped_reads', aggfunc='sum'),
        total_unmapped_reads=pd.NamedAgg(column='unmapped_reads', aggfunc='sum'),
        genome_length=pd.NamedAgg(column='len', aggfunc='sum')
    )

    print(f'{datetime.now().strftime("%H:%M:%S")}: Importing all QC...')
    middle_id_item = set(stats_table.s.str.split('_').map(lambda x: x[1]))
    if len(middle_id_item) != 1:
        raise ValueError('ERROR: only a single ID second term is supported.')
    else:
        middle_id_item = list(middle_id_item)[0]
    downloaded_qc_files['s_mod'] = downloaded_qc_files.s.str.split('_').map(lambda x: x[0] + '_' + middle_id_item + '_' + x[2] + '_' + x[3])
    qc_table = import_and_cat_tables(downloaded_qc_files, 's_mod', 'qc', 's', append_ids_and_t=True, filter_by=list(stats_table['s']))
    final_stats_table = stats_table.merge(qc_table, how='outer', on='s'
                                  ).merge(yield_table, how='inner', on=['s','batch']
                                  ).merge(idxstats_summary, how='inner', on='s')
    final_stats_table = compute_nuc_coverage(final_stats_table, 'nuc_mean_coverage')

    # output
    print(f'{datetime.now().strftime("%H:%M:%S")}: Outputting flat files...')
    downloaded_files.to_csv(re.sub('ht$', 'tsv', file_paths_table_output), sep='\t', index=False)
    final_stats_table.to_csv(re.sub('ht$', 'tsv', per_sample_stats_output), sep='\t', index=False)

    # output to sql
    print(f'{datetime.now().strftime("%H:%M:%S")}: Generating Hail tables and outputting to SQL database...')
    ht_files = hl.Table.from_pandas(downloaded_files).repartition(5).key_by('batch')
    #ht_files = hl.import_table('file://' + os.getcwd() + '/' + re.sub('ht$', 'tsv', file_paths_table_output), impute=True, min_partitions=5)
    ht_files.write(f'dnax://{my_database}/{file_paths_table_output}', overwrite=True)
    
    #final_stats_table = pd.read_csv(os.getcwd() + '/' + re.sub('ht$', 'tsv', per_sample_stats_output), sep='\t')
    final_stats_table['genetic_sex'] = final_stats_table.genetic_sex.astype(str)
    for field in ['discordance_prc', 'read_haps_error_percentage', 'freemix_percentage', 'yield', 'prc_proper_pairs', 'prc_auto_ge_15x']:
        if field in final_stats_table.columns:
            final_stats_table[field] = final_stats_table[field].map(lambda x: np.nan if x == 'None' else x).astype(float)
    for field in ['alias']:
        if field in final_stats_table.columns:
            final_stats_table[field] = final_stats_table[field].map(lambda x: np.nan if x == 'None' else x).astype(str)
    ht_stats = hl.Table.from_pandas(final_stats_table).checkpoint(f'dnax://{my_database}/tmp2/stats_output_temp.ht')
    ht_stats = ht_stats.repartition(50).key_by('s')
    #ht_stats = hl.import_table('file:///' + os.getcwd() + '/' + re.sub('ht$', 'tsv', per_sample_stats_output), impute=True, min_partitions=50)
    ht_stats.write(f'dnax://{my_database}/{per_sample_stats_output}', overwrite=True)


parser = argparse.ArgumentParser()
parser.add_argument('--pipeline-output-folder', type=str, required=True, 
                    help="Folder containing folders, each of which should contain pipeline outputs. Do not include project name. Should contain VCF, coverage, and diagnostic files. Can be a comma-delimited list.")
parser.add_argument('--file-paths-table-output', type=str, required=True, 
                    help="Local path to ht to output which will contain all paths for coverages and vcfs.")
parser.add_argument('--per-sample-stats-output', type=str, required=True, 
                    help="Local path to ht to output which contains all per-sample data.")
parser.add_argument('--dx-init', type=str, required=True,
                    help='SQL database path for use in DNAnexus.')

parser.add_argument('--unified-prefix', type=str, default='batch_', help='Prefix for all other files.')
parser.add_argument('--vcf-suffix', type=str, default='batch_merged_mt_calls.vcf.bgz',
                    help="Suffix of each final VCF to import. Expects this to be multi-sample.")
parser.add_argument('--coverage-suffix', type=str, default='batch_merged_mt_coverage.tsv.bgz',
                    help="Suffix of each coverage tsv to import. Expects multi-sample tables.")
parser.add_argument('--mtstats-suffix', type=str, default='batch_analysis_statistics.tsv',
                    help="Suffix of each mtPipeline statistics file to import. Expects multi-sample analysis.")
parser.add_argument('--yield-suffix', type=str, default='batch_yield_metrics.tsv.gz',
                    help="Suffix of each coverage tsv to import. Expects multi-sample analysis.")
parser.add_argument('--idxstats-suffix', type=str, default='batch_idxstats_metrics.tsv.gz',
                    help="Suffix of each mtPipeline statistics file to import. Expects multi-sample analysis.")
parser.add_argument('--qc-stats-folder', type=str, default='Bulk/Whole genome sequences/Concatenated QC Metrics/',
                    help="Folder containing folders, each of which should contain QC data from WGS. Also supports a single folder with relevant files in it. Do not include project name.")
parser.add_argument('--qc-suffix', type=str, default='.qaqc_metrics',
                    help="Suffix of each WGS QC file to import. Assumes that this file contains multiple rows for a single sample.")

parser.add_argument('--avoid-filtering-idxstats-chr', action='store_true',
                    help='If enabled, this flag will prevent filtering to autosomes when producing nucDNA coverage estimates.')


# defaults for debugging
pipeline_output_folder = '220618_MitochondriaPipelineSwirl/v2.5_Multi_first50/,220618_MitochondriaPipelineSwirl/20k/'
pipeline_output_folder = '220618_MitochondriaPipelineSwirl/20k/,220618_MitochondriaPipelineSwirl/next24500/,220618_MitochondriaPipelineSwirl/next10k_upto55k/,220618_MitochondriaPipelineSwirl/next45k_upto100k/,220618_MitochondriaPipelineSwirl/next8500_upto110k/,220618_MitochondriaPipelineSwirl/next8500_upto110k_p2/,220618_MitochondriaPipelineSwirl/next8500_upto110k_p3/,220618_MitochondriaPipelineSwirl/next28k_upto137k/,220618_MitochondriaPipelineSwirl/next9k_upto146k/,220618_MitochondriaPipelineSwirl/next27k_upto173k/,220618_MitochondriaPipelineSwirl/final27k_upto200k/'
vcf_suffix = 'batch_merged_mt_calls.vcf.bgz'
coverage_suffix = 'batch_merged_mt_coverage.tsv.bgz'
mtstats_suffix = 'batch_analysis_statistics.tsv'
yield_suffix = 'batch_yield_metrics.tsv.gz'
idxstats_suffix = 'batch_idxstats_metrics.tsv.gz'
qc_stats_folder = 'Bulk/Whole genome sequences/Concatenated QC Metrics/'
qc_suffix = '.qaqc_metrics'
file_paths_table_output = 'tab_batch_file_paths.ht'
per_sample_stats_output = 'tab_per_sample_stats.ht'
dx_init = '220619_MitochondriaPipelineSwirl_v2_5_Multi_20k'
dx_init = '220722_MitochondriaPipelineSwirl_v2_5_Multi_200k'
avoid_filtering_idxstats_chr = False
unified_prefix = 'batch_'


if __name__ == '__main__':
    args = parser.parse_args()
    main(**vars(args))

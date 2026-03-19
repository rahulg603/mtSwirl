# Pipeline for the merging of mtSwirl pipeline outputs in AllofUs.
# Imports and appends QC file information as well.
# see QC field explanations here:
# https://aousupporthelp.zendesk.com/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized

import hail as hl
import pandas as pd
import argparse
import multiprocessing
import re, os


def reader(args):
    idx, f = args
    tab = pd.read_csv(f, index_col=None, header=0, sep='\t').assign(newcolthis=idx)
    return tab


def import_and_cat_tables(directory_df, id_col, path_col, new_id_col, max_threads=2, filter_by=None, enforce_nonmiss=False):
    ids = [x for _, x in directory_df[id_col].iteritems()]
    directories = [x for _, x in directory_df[path_col].iteritems()]
    p = multiprocessing.Pool(processes=max_threads)

    to_subset_to = filter_by if filter_by is not None else ids
    empty_df = pd.DataFrame({new_id_col:[]})
       
    #df_from_each_file = p.map(reader, [(idx, f) for idx, f in zip(ids, directories) if idx in to_subset_to])
    df_from_each_file = [reader((idx, f)) for idx, f in zip(ids, directories) if idx in to_subset_to]
    
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


def main(pipeline_output_path, qc_stats_path, file_paths_table_output, 
         per_sample_stats_output, file_paths_table_flat_output, per_sample_stats_flat_output):

    # import pipeline output file
    pipeline_output_file = pd.read_csv(pipeline_output_path, sep='\t')
    pipeline_output_file = pipeline_output_file.rename({'merged_calls': 'vcf', 
                                                        'merged_coverage':'coverage', 
                                                        'merged_statistics': 'stats', 
                                                        'subpath_id':'batch'}, axis=1)

    # download mito pipeline data
    print('Obtaining QC stats...')
    qc_stats = hl.import_table(qc_stats_path, impute=True).to_pandas() 
    qc_stats = qc_stats.rename({'research_id': 's',
                                'verify_bam_id2_contamination': 'freemix_percentage',
                                'mean_coverage': 'nuc_mean_coverage'}, axis=1)
    
    # checks on imported pipeline results
    print('Running checks...')
    batches = list(pipeline_output_file.batch)
    if len(set(batches)) != len(batches):
        raise ValueError('ERROR: there are duplicate batches.')
    basename_log = [os.path.dirname(x) for x in pipeline_output_file.merging_log]
    basename_vcf = [os.path.dirname(x) for x in pipeline_output_file.vcf]
    basename_cov = [os.path.dirname(x) for x in pipeline_output_file.coverage]
    basename_stats = [os.path.dirname(x) for x in pipeline_output_file.stats]
    if not all([(a == b) and (a == c) and (a == d) for a,b,c,d in zip(basename_log, basename_vcf, basename_cov, basename_stats)]):
        raise ValueError('ERROR: all files should have the same path within a batch.')
    if not all([re.search(search, x) for x, search in zip(basename_vcf, pipeline_output_file.batch)]):
        raise ValueError('ERROR: all files should have paths matching the batch.')

    # import stats and qc and merge all into table
    print('Importing run statistics...')
    stats_table = import_and_cat_tables(pipeline_output_file, 'batch', 'stats', 'batch', enforce_nonmiss=False)

    print('Running checks on run statistics...')
    samples = list(stats_table.s)
    if len(set(samples)) != len(samples):
        raise ValueError('ERROR: there are duplicate samples.')

    # join with qc stats
    final_stats_table = stats_table.merge(qc_stats, how='left', on='s')
    samples2 = list(final_stats_table.s)
    if len(set(samples2)) != len(samples2):
        raise ValueError('ERROR: there are duplicate samples in the merged stats / QC table.')

    # output
    print('Outputting flat files...')
    pipeline_output_file.to_csv(file_paths_table_flat_output, sep='\t', index=False)
    final_stats_table.to_csv(per_sample_stats_flat_output, sep='\t', index=False)

    # output to hail tables
    print('Generating Hail tables and outputting to Google Cloud...')
    ht_files = hl.Table.from_pandas(pipeline_output_file).repartition(5).key_by('batch')
    ht_files.write(file_paths_table_output, overwrite=True)
    
    ht_stats = hl.Table.from_pandas(final_stats_table).repartition(50).key_by('s')
    ht_stats.write(per_sample_stats_output, overwrite=True)


parser = argparse.ArgumentParser()
parser.add_argument('--pipeline-output-path', type=str, required=True, 
                    help="This should come from cromwell_run_monitor.py, which outputs a list of successfully analyzed batches.")
parser.add_argument('--file-paths-table-output', type=str, required=True, 
                    help="gs:// path to ht to output which will contain all paths for coverages and vcfs.")
parser.add_argument('--per-sample-stats-output', type=str, required=True, 
                    help="gs:// path to ht to output which contains all per-sample data.")
parser.add_argument('--file-paths-table-flat-output', type=str, required=True, 
                    help="gs:// path to tsv to output which will contain all paths for coverages and vcfs.")
parser.add_argument('--per-sample-stats-flat-output', type=str, required=True, 
                    help="gs:// path to tsv to output which contains all per-sample data.")

parser.add_argument('--qc-stats-path', type=str, default='gs://fc-aou-datasets-controlled/v6/wgs/vcf/aux/qc/genomic_metrics.tsv',
                    help="Folder containing WGS QC statistics. Contains all relevant data.")


# defaults for debugging
pipeline_output_path = 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/220717_test_runs_120_take2/success.tsv'
qc_stats_path = 'gs://fc-aou-datasets-controlled/v6/wgs/vcf/aux/qc/genomic_metrics.tsv'
file_paths_table_output = 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/ht/test/tab_batch_file_paths.ht'
per_sample_stats_output = 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/ht/test/tab_per_sample_stats.ht'
file_paths_table_flat_output = 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/tsv/test/tab_batch_file_paths.tsv'
per_sample_stats_flat_output = 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/tsv/test/tab_per_sample_stats.tsv'


if __name__ == '__main__':
    args = parser.parse_args()
    main(**vars(args))

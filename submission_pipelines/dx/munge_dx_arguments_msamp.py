from copy import deepcopy
import json
import pandas as pd
import argparse
import os, re
from math import ceil
from utils.util_funcs import read_file_at_path

def split_dataframe_by_position(df, n_each, 
                                m2_ignore_path='./utils/mito_02_completed/all_finished_sample_names.txt'):
    """
    Takes a dataframe and an integer of the number of splits to create.
    Returns a list of dataframes.
    """
    m2_ignore_path = "./samples_2320_comparison.txt"
    m2_ignore = set([f"{s}_23372_0_0" for s in read_file_at_path(m2_ignore_path)])
    # m3_ignore = set([f"{s}_23372_0_0" for s in read_file_at_path(m3_ignore_path)])
    ignore_list = m2_ignore
    print(ignore_list)
    filtered_df = df[df['stage-common.sample_name'].isin(ignore_list)].reset_index(drop=True)
    
    # Calculate number of splits needed
    n_splits = ceil(len(filtered_df) / n_each)
    
    # List to store smaller DataFrames
    dataframes = []
    
    # Split the DataFrame into smaller DataFrames
    for i in range(n_splits):
        start = i * n_each
        end = min((i + 1) * n_each, len(filtered_df))  # Ensure end doesn't exceed the length of the DataFrame
        temp_df = filtered_df.iloc[start:end, :]  # Slice the DataFrame
        dataframes.append(temp_df)
    return dataframes

def link_to_dct(vec):
    matches = [re.search('(project-.+):(file-.+)', str(x)) for x in vec]
    if not all(matches):
        raise ValueError("ERROR: all DNAnexus links should have project-...:file- format.")
    return [{'$dnanexus_link': {'project': x[1], 'id': x[2]}} for x in matches]


def dct_to_link(vec):
    return [x['$dnanexus_link']['project'] + ':' + x['$dnanexus_link']['id'] for x in vec]


def format_jsontable_for_processing(table):
    for col in table.columns:
        this_col = list(table[col])
        is_dict = [type(x) is dict for x in this_col]
        if all(is_dict):
            tf_dx = all(['$dnanexus_link' in y for y in this_col])
            if(tf_dx):
                table[col] = dct_to_link(this_col)
    return(table)


def table_to_dx_json(table):
    tf_cols_dx_link = table.apply(lambda x: x.str.contains('project-.+:file-.+$'), 1).apply(lambda x: all(x), 0)
    cols_dx_link = list(tf_cols_dx_link[tf_cols_dx_link].index)
    cols_not_dx = list(tf_cols_dx_link[~tf_cols_dx_link].index)
    dct_out = {x: list(table[x]) for x in cols_not_dx}
    print(cols_dx_link)
    dct_dx = {x: link_to_dct(list(table[x])) for x in cols_dx_link}
    dct_out.update(dct_dx)
    return(dct_out)


parser = argparse.ArgumentParser()
parser.add_argument('--files', type=str, help='Path to folder with files. Will use columns with stage_common.+ID.', default=None)
parser.add_argument('--specific-file', type=str, default=None, help='If rerunning a particular file, this will further split that file according to specified parameters.')
parser.add_argument('--reshape-folder', type=str, help='This will allow input of a folder of .json files reshaping to a run size of interest.')
parser.add_argument('--out', type=str, required=True, help='Path for outputting batch files.')
parser.add_argument('--n-each', type=int, default=10, help='Number of samples per job.')
parser.add_argument('--n-max', type=int, default=20, help='Number of samples to analyze total.')
parser.add_argument('--base-json', type=str, required=True, help='Path to json file to append to.')
parser.add_argument('--make-name-col', type=str, required=True, help='Name of field to use for sample name.')
parser.add_argument('--target-name-col', type=str, default='sample_name', help='Name of field to make for sample name.')
parser.add_argument('--job-base-name', type=str, default='MitochondriaPipelineMSamp_', help='Base name for each job.')

if __name__ == '__main__':
    args = parser.parse_args()
    
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    with open(args.base_json, 'r') as json_file:    
        base_json = json.load(json_file)

    if args.files is not None:
        if not os.path.exists(args.files):
            raise argparse.ArgumentError('ERROR: --files must point to an existing directory.')

        files_to_read = [x for x in os.listdir(args.files) if re.search('.+[.]tsv$', x)]
        tables_loaded = [pd.read_csv(args.files + '/' + x, sep='\t') for x in files_to_read]
        tables_cat = pd.concat(tables_loaded, axis=0)
        to_use = {x: re.sub(' ID$', '', x) for x in tables_cat.columns if re.search('stage-common.+ ID$', x)}
        to_use.update({args.make_name_col: f'stage-common.{args.target_name_col}'})
        col_filtered_table = tables_cat[to_use.keys()].rename(to_use, axis=1)
        col_subset_table = col_filtered_table.head(min(args.n_max, col_filtered_table.shape[0]))

    elif args.specific_file is not None:
        # Here we read in a JSON and then convert to "col_subset_table" which is a data frame
        with open(args.specific_file, 'r') as import_json:
            imported_json = json.load(import_json)
        col_subset_table = pd.DataFrame({k:v for k, v in imported_json.items() if k not in base_json.keys()})
        col_subset_table = format_jsontable_for_processing(col_subset_table)
        print(col_subset_table)

    else:
        # Now we read in a folder of JSONs
        if not os.path.exists(args.reshape_folder):
            raise argparse.ArgumentError('ERROR: --reshape-folder must point to an existing directory.')
        table_list = []
        for j_fl in os.listdir(args.reshape_folder):
            if re.search('.json$', j_fl):
                with open(args.reshape_folder + j_fl, 'r') as import_json:
                    imported_json = json.load(import_json)
                this_col_subset_table = pd.DataFrame({k:v for k, v in imported_json.items() if k not in base_json.keys()})
                this_col_subset_table = format_jsontable_for_processing(this_col_subset_table)
                table_list.append(this_col_subset_table)
        col_subset_table = pd.concat(table_list)
        print(col_subset_table.count())

    if len(set(col_subset_table[f'stage-common.{args.target_name_col}'])) != col_subset_table.shape[0]:
        raise ValueError('ERROR: duplicate sample names detected.')

    per_job_tables = split_dataframe_by_position(col_subset_table, args.n_each)
    
    k_in_file = list(base_json.keys())
    if any([k in k_in_file for k in col_subset_table.columns]):
        raise argparse.ArgumentError('ERROR: keys to add to JSON are already in the JSON.')

    for idx, table in enumerate(per_job_tables):
        this_json = deepcopy(base_json)
        this_json.update(table_to_dx_json(table))
        this_path = args.out + '/' + args.job_base_name + str(idx) + '.json'
        with open(this_path, 'w') as new_file:
            json.dump(this_json, new_file)


    # Evaluating performance
    # Note that --reshape-folder does not preserve order, but it does hit all of the samples of interest
    # sample_nm = []
    # cram = []
    # crai = []
    # for fl in os.listdir('/Users/rahulgupta/Desktop/UKB_DX_Runners/dx_jobs_100_220618_updatedjson20_next10k_reshape90/'):
    #     with open('/Users/rahulgupta/Desktop/UKB_DX_Runners/dx_jobs_100_220618_updatedjson20_next10k_reshape90/' + fl) as this_j:
    #         jfile = json.load(this_j)
    #         sample_nm.extend(jfile['stage-common.sample_name'])
    #         cram.extend(jfile['stage-common.wgs_aligned_input_bam_or_cram'])
    #         crai.extend(jfile['stage-common.wgs_aligned_input_bam_or_cram_index'])

    # sample_nm_old = []
    # cram_old = []
    # crai_old = []
    # for fl in os.listdir('/Users/rahulgupta/Desktop/UKB_DX_Runners/dx_jobs_100_220618_updatedjson20_next10k/'):
    #     with open('/Users/rahulgupta/Desktop/UKB_DX_Runners/dx_jobs_100_220618_updatedjson20_next10k/' + fl) as this_j:
    #         if re.search('.json$', fl):
    #             jfile = json.load(this_j)
    #             sample_nm_old.extend(jfile['stage-common.sample_name'])
    #             cram_old.extend(jfile['stage-common.wgs_aligned_input_bam_or_cram'])
    #             crai_old.extend(jfile['stage-common.wgs_aligned_input_bam_or_cram_index'])

    # new = pd.DataFrame({'samp': sample_nm, 'cram': cram, 'crai': crai})
    # news = new.sort_values(by = 'samp', axis=0).reset_index(drop=True)
    # old = pd.DataFrame({'samp': sample_nm_old, 'cram': cram_old, 'crai': crai_old})
    # olds = old.sort_values(by = 'samp', axis=0).reset_index(drop=True)
    # all(news==olds)

    # fn_sort = lambda x: x['$dnanexus_link']['id']
    # sample_nm.sort(); sample_nm_old.sort(); cram.sort(key=fn_sort); cram_old.sort(key=fn_sort); crai.sort(key=fn_sort); crai_old.sort(key=fn_sort)
    # sample_nm == sample_nm_old
    # cram == cram_old
    # crai == crai_old

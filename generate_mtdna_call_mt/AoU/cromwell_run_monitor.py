from google.cloud import storage
from datetime import datetime
import pandas as pd
import os, re
import argparse


SHARD_SEARCH = 'shard-[0-9]{1,}'
NMERGE = 3

# get main order
known_order = ['call-SubsetBamToChrMAndRevert', 'call-AlignAndCallR1', 'call-ProduceSelfRefFiles', 
                'call-AlignAndCallR2', 'call-LiftOverAfterSelf', 'call-LiftoverSelfCoverage']
has_subtasks = dict(zip(known_order, [False, True, True, True, False, False]))
map_task_to_order = {x: idx for idx, x in enumerate(known_order)}
subtask_folders = {'call-ProduceSelfRefFiles':'ProduceSelfReferenceFiles', 'call-AlignAndCallR1':'AlignAndCallR1', 'call-AlignAndCallR2':'AlignAndCallR2'}

# subtasks order
subtasks_order = {'call-AlignAndCallR1':['call-CallNucHCIntegrated', 'call-CallNucM2Integrated', 'call-CallMt', 'call-GetContamination', 'call-FilterContamination'], 
                  'call-ProduceSelfRefFiles':['call-ProduceSelfReference', 'call-ChainSwapLiftoverBed'], 
                  'call-AlignAndCallR2':['call-AlignToMtRegShiftedAndMetrics', 'call-CallMtAndShifted', 'call-LiftoverCombineMergeFilterContamSplit']}
subtasks_map_to_order = {k: {x: idx for idx, x in enumerate(lst)} for k, lst in subtasks_order.items()}


def list_gcs_directories(client, bucket, prefix):
    # from https://github.com/GoogleCloudPlatform/google-cloud-python/issues/920
    iterator = client.list_blobs(bucket, prefix=prefix, delimiter='/')
    prefixes = set()
    for page in iterator.pages:
        prefixes.update(page.prefixes)
    return prefixes


def update_path_for_attempts(client, bucket, path, blobs):
    # check for any folders labeled attempt
    dirs = list_gcs_directories(client, bucket, path)
    attempts = [os.path.basename(os.path.dirname(x)) for x in dirs if re.search('attempt-[0-9]{1,}$', os.path.basename(os.path.dirname(x)))]
    
    if len(attempts) > 0:
        # if found, find the folder with the latest attempt
        latest_attempt = max([int(re.search('attempt-([0-9]{1,})$',x).group(1)) for x in attempts])

        # obtain a new path containing the latest attempt
        path = path + 'attempt-' + str(latest_attempt) + '/'

        # check for new blobs
        blobs = list(bucket.list_blobs(prefix=path, delimiter='/'))
    
    return blobs, path, len(attempts) > 0


def test_success(client, bucket, path):
    """
    A completed job has an RC file with code 0 and a finished log file.
    If RC is present and the log file has failed that will explain why.
    If RC is missing, then the job likely is in progress and we will allow multiple log blobs. We will still check for failure.
    Added routine to check for preempted runs.
    """
    blobs = list(bucket.list_blobs(prefix=path, delimiter='/'))
    rc_found = 'rc' in [os.path.basename(x.name) for x in blobs]
    if rc_found:
        rc_blob = [x for x in blobs if os.path.basename(x.name) == 'rc'][0]
        rc_out = (rc_blob.download_as_string().decode('utf8').split('\n'))[0]
    else:
        blobs, path, tf_update = update_path_for_attempts(client, bucket, path, blobs)
        if tf_update and ('rc' in [os.path.basename(x.name) for x in blobs]):
            rc_found = True
            rc_blob = [x for x in blobs if os.path.basename(x.name) == 'rc'][0]
            rc_out = (rc_blob.download_as_string().decode('utf8').split('\n'))[0]
        else:
            rc_out = 'NA'

    log_blobs = [(x, x.name) for x in bucket.list_blobs(prefix=path, delimiter='/') if re.search('.log$', os.path.basename(x.name))]
    if rc_found and (len(log_blobs) != 1):
        raise ValueError('ERROR: there should be exactly 1 log file in each folder for testing when rc was returned.')
    elif len(log_blobs) > 1:
        raise ValueError('ERROR: there should not be >1 log file in the search folder.')
    
    if (len(log_blobs) == 0):
        log_path = 'NA'
        pass_log = False
        fail_log = False
    else:
        this_log = log_blobs[0][0].download_as_string().decode('utf8').split('\n')
        log_path = 'gs://' + bucket.name + '/' + log_blobs[0][0].name
        item_to_check = this_log[len(this_log)-2]
        pass_log = re.search('Done delocalization.$', item_to_check)
        fail_log = re.search('^Required file output.+does not exist.$', item_to_check)
    pass_rc = rc_found and (rc_out == '0')
    if pass_log and pass_rc and not fail_log:
        return 'PASS', len(list(bucket.list_blobs(prefix=path+'out/', delimiter='/'))), log_path, rc_out
    elif fail_log and (not pass_log or not pass_rc):
        return 'FAIL', 0, log_path, rc_out
    elif not fail_log and not pass_log:
        # fine for the pass_rc to be True here, the log is what matters
        return 'IN PROGRESS', 0, log_path, rc_out
    else:
        raise ValueError('Does not make sense for the localization file to have both PASS and FAIL flags.')


def process_single_run(storage_client, bucket, this_folder, run_folder, success_only):
    run_prefix = f'{run_folder}/MitochondriaPipelineWrapper/{this_folder}/'

    # get all subfolders
    subfolders = list_gcs_directories(storage_client, bucket, run_prefix)
    dirnames = [os.path.basename(os.path.dirname(x)) for x in subfolders]

    # check for merged files
    nonshard_folders = [path for path, main_dir in zip(subfolders, dirnames) if not re.search('^call-MitochondriaPipeline_v2_5$', main_dir)]

    if (len(nonshard_folders) == 0) | (not any([re.search('call-MergeMitoMultiSampleOutputsInternal.$', x) for x in nonshard_folders])):
        merge_not_run = True
        merge_log = 'NA'
        merge_res = 'NA'
        merge_items = 0
    else:
        merge_not_run = False
        merge_res, merge_items, merge_log, _ = test_success(storage_client, bucket, run_prefix + 'call-MergeMitoMultiSampleOutputsInternal/')

    if success_only:
        shards = [None]
    else:
        # get shard subfolders
        shard_prefix = f'{run_prefix}call-MitochondriaPipeline_v2_5/'
        shard_subfolders = list_gcs_directories(storage_client, bucket, shard_prefix)
        shard_dirnames = [os.path.basename(os.path.dirname(x)) for x in shard_subfolders]
        shards = [path for path, main_dir in zip(shard_subfolders, shard_dirnames) if re.search(SHARD_SEARCH, main_dir)]

    return shards, merge_not_run, merge_items, merge_res, merge_log


def process_merging_run(bucket, this_folder, merge_items, merge_log):
    if merge_items != 0:
        raise ValueError(f'ERROR: there should be {str(0)} items in the output folder from merging.')
    merged_prefix = '/'.join(os.path.dirname(merge_log).split('/')[3:])
    files_output_merge = [os.path.basename(x.name) for x in bucket.list_blobs(prefix=merged_prefix+'/', delimiter='/')]
    gs_prefix = f'gs://{bucket.name}'
    
    if len([x for x in files_output_merge if re.search('^batch_', x)]) != NMERGE:
        raise ValueError(f'ERROR: there should be {str(NMERGE)} items with the batch_ prefix as products from merging.')
    if 'batch_merged_mt_calls.vcf.bgz' not in files_output_merge:
        raise ValueError(f'ERROR: variant calls not found in batch {this_folder}')
    if 'batch_merged_mt_coverage.tsv.bgz' not in files_output_merge:
        raise ValueError(f'ERROR: coverage file not found in batch {this_folder}')
    if 'batch_analysis_statistics.tsv' not in files_output_merge:
        raise ValueError(f'ERROR: run statistics file not found in batch {this_folder}')
    
    df = pd.DataFrame({'merging_log': [merge_log], 'merged_calls': f'{gs_prefix}/{merged_prefix}/batch_merged_mt_calls.vcf.bgz', 
                        'merged_coverage': f'{gs_prefix}/{merged_prefix}/batch_merged_mt_coverage.tsv.bgz', 
                        'merged_statistics': f'{gs_prefix}/{merged_prefix}/batch_analysis_statistics.tsv'})

    return df


def obtain_latest_shard_run(storage_client, bucket, shards):

    # round 1: make sure there are the same number of elements per shard
    prefixes1 = {os.path.basename(os.path.dirname(x)): x+'MitochondriaPipeline/' for x in list(shards)}
    subpaths1 = {k: list(list_gcs_directories(storage_client, bucket, v)) for k,v in prefixes1.items()}
    if not all([len(v) == 1 for _, v in subpaths1.items()]):
        raise ValueError('ERROR: all sub-shards should have only 1 path.')

    # round 2: check tasks completed per shard
    prefixes2 = {k: v[0] for k,v in subpaths1.items()}
    subpaths2 = {k: list(list_gcs_directories(storage_client, bucket, v)) for k,v in prefixes2.items()}
    tasks_run = {k: [os.path.basename(os.path.dirname(x)) for x in v] for k,v in subpaths2.items()}
    latest_run = {x: known_order[max([map_task_to_order[spec_task] for spec_task in tasks])] for x, tasks in tasks_run.items()}

    # round 3: within each latest run task for each shard, check which was successful
    dct_success_test = {}
    for shard_id, latest_task in latest_run.items():
        if has_subtasks[latest_task]:
            this_path = prefixes2[shard_id] + latest_task + '/' + subtask_folders[latest_task] + '/'
            this_order_mapping = subtasks_map_to_order[latest_task]
            this_ordering = subtasks_order[latest_task]
            subpaths_subtask = list(list_gcs_directories(storage_client, bucket, this_path))
            if len(subpaths_subtask) != 1:
                raise ValueError('ERROR: all sub-shards should have only 1 path.')
    
            subtasks_paths = list(list_gcs_directories(storage_client, bucket, subpaths_subtask[0]))
            subtasks = {os.path.basename(os.path.dirname(x)):x for x in subtasks_paths}
            subtask_latest_run = this_ordering[max([this_order_mapping[spec_task] for spec_task, _ in subtasks.items()])]
            dct_success_test.update({shard_id: {subtask_latest_run: test_success(storage_client, bucket, subpaths_subtask[0] + subtask_latest_run + '/')}})

        else:
            dct_success_test.update({shard_id: {latest_task: test_success(storage_client, bucket, prefixes2[shard_id] + latest_task + '/')}})
    
    return dct_success_test


def check_success_single(storage_client, bucket, this_folder, run_folder, success_only):
    shards, merge_not_run, merge_items, merge_res, merge_log = process_single_run(storage_client, bucket, this_folder, run_folder, success_only)

    if (not merge_not_run) and (merge_res == 'PASS'):
        failed = False       
        df = process_merging_run(bucket, this_folder, merge_items, merge_log)
    elif (not merge_not_run) and (merge_res == 'FAIL'):
        failed = True
        df = pd.DataFrame({'shard': ['NA'], 'latest_task': ['call-MergeMitoMultiSampleOutputsInternal'], 'status': ['FAIL'], 'log_file': [merge_log]})
    elif (not merge_not_run) and (merge_res == 'IN PROGRESS'):
        failed = True
        df = pd.DataFrame({'shard': ['NA'], 'latest_task': ['call-MergeMitoMultiSampleOutputsInternal'], 'status': ['IN PROGRESS'], 'log_file': [merge_log]})
    elif merge_not_run:
        if success_only:
            failed = True
            df = pd.DataFrame({'shard': ['NA'], 'status': ['NOT YET MERGING'], 'log_file': ['NA']})
        else:
            failed = True
            dct_success_test = obtain_latest_shard_run(storage_client, bucket, shards)

            # filter to in progress
            in_progress = [(k, call, item[2]) for k, v in dct_success_test.items() for call, item in v.items() if item[0] == 'IN PROGRESS']

            # filter to failures
            failures = [(k, call, item[2], item[3]) for k, v in dct_success_test.items() for call, item in v.items() if item[0] == 'FAIL']

            # filter to successes
            passes = [(k, call, item[2], item[3]) for k, v in dct_success_test.items() for call, item in v.items() if item[0] == 'PASS']

            # generate table
            df_prog = pd.DataFrame({'shard': [x for x, _, _ in in_progress], 'latest_task': [x for _, x, _ in in_progress], 
                                    'log_file': [x for _, _, x in in_progress], 'status': ['IN PROGRESS' for _ in in_progress], 
                                    'return_code': ['NA' for _ in in_progress]})
            df_fail = pd.DataFrame({'shard': [x for x, _, _, _ in failures], 'latest_task': [x for _, x, _, _ in failures], 
                                    'log_file': [x for _, _, x, _ in failures], 'status': ['FAIL' for _ in failures], 
                                    'return_code': [x for _, _, _, x in failures]})
            df_pass = pd.DataFrame({'shard': [x for x, _, _, _ in passes], 'latest_task': [x for _, x, _, _ in passes], 
                                    'log_file': [x for _, _, x, _ in passes], 'status': ['PASS' for _ in passes], 
                                    'return_code': [x for _, _, _, x in passes]})
            df = pd.concat([df_prog, df_fail, df_pass], axis=0).reset_index(drop=True)
    
    df['subpath_id'] = this_folder
    return df, not failed


def count_shards_single(storage_client, bucket, this_folder, run_folder):
    shards, _, _, _, _ = process_single_run(storage_client, bucket, this_folder, run_folder)
    return len(shards)


def produce_sample_lists(sub_id_file, sample_list_placeholder):
    with open(sub_id_file, 'r') as sub_ids:
        ids = sub_ids.read().splitlines()
    
    df_list = []
    batch_run_number = 0
    for id in ids:
        this_sample_loc = sample_list_placeholder.format(str(batch_run_number))
        with open(this_sample_loc, 'r') as sample_id_list:
            sample_ids = sample_id_list.read().splitlines()
        df = pd.DataFrame({'batch': id, 'batch_run_number': batch_run_number, 's': sample_ids})
        df['shard'] = [f'shard-{str(x)}' for x in range(0, len(sample_ids))]
        df_list.append(df)
        batch_run_number+=1
    
    df_final = pd.concat(df_list, axis=0).reset_index(drop=True)
    return df_final


def print_log(success_only):
    print('')
    print('HACKY CROMWELL VISUALIZER')
    print('July 6th 2022, Rahul Gupta')
    print('Updated August 18th 2022')
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f'Start time: {current_time}')
    if success_only:
        print(f'NOTE: running in a shortcut mode, outputting only successfully completed batches.')
    print('---------------------------------')


parser = argparse.ArgumentParser()
parser.add_argument('--run-folder', type=str, help='Path to Cromwell run folder.')
parser.add_argument('--check-success', action='store_true', help='If enabled, will return information on run completion. Files are saved to disk.')
parser.add_argument('--success-only', action='store_true', help='If enabled will only find folders that were successfully run. Speeds up the run monitor.')
parser.add_argument('--get-shard-count', action='store_true', help='Returns the number of shards generated.')
parser.add_argument('--hide-pass', action='store_true', help='Hides pass (but stalled) shards from output.')
parser.add_argument('--subfolder', type=str, default=None, help='Run path. If provided, will only analyze a single folder. If not enabled, will iterate through all folders in the --run-folder.')
parser.add_argument('--sub-ids', type=str, help='Path to workflow IDs. Assumes this has the correct count and is ordered by the workflow number.')
parser.add_argument('--sample-lists', type=str, help='Path to sample lists for analysis. Place {} where the workflow ID is.')
parser.add_argument('--output', type=str, help='Prefix for output files. Only used if --check-success is enabled.')


if __name__ == '__main__':
    args = parser.parse_args()
    print_log(args.success_only)
    if args.success_only and args.get_shard_count:
        raise argparse.ArgumentError('ERROR: cannot enable both --success-only and --get-shard-count.')
    if args.success_only and args.hide_pass:
        raise argparse.ArgumentError('ERROR: cannot enable both --success-only and --hide-pass.')
    
    bucket_id = os.path.basename(os.getenv("WORKSPACE_BUCKET"))
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_id)

    if args.subfolder is None:
        subfolders_to_test = [os.path.basename(os.path.dirname(x)) for x in list_gcs_directories(storage_client, bucket, f'{args.run_folder}/MitochondriaPipelineWrapper/')]
    else:
        subfolders_to_test = [args.subfolder]

    df_sample = produce_sample_lists(args.sub_ids, args.sample_lists)
    df_sample.to_csv(f'{args.output}.sample_mapping.tsv', sep='\t', index=False)
    
    if args.check_success:
        for suffix in ['success', 'failure', 'running', 'stalled']:
            if os.path.exists(f'{args.output}.{suffix}.tsv'):
                os.remove(f'{args.output}.{suffix}.tsv')

        res = [check_success_single(storage_client, bucket, x, args.run_folder, args.success_only) for x in subfolders_to_test]
        print(f'Of {str(len(res))} runs, {str(sum([success for _, success in res]))} are finished.')
        if sum([success for _, success in res]) > 0:
            for batch, this_res in zip(subfolders_to_test, res):
                if this_res[1]:
                    print(batch)
            df_success = pd.concat([tb for tb, success in res if success], axis=0)
            df_success.to_csv(f'{args.output}.success.tsv', sep='\t', index=False)

        if sum([not success for _, success in res]) > 0:
            print(f'{str(sum([not success for _, success in res]))} are incomplete.')
            if not args.success_only:
                df_incomplete = pd.concat([tb for tb, success in res if not success], axis=0)
                if sum([not success for _, success in res]) > 0:
                    print('The current status of ongoing batches is as follows:')
                    running_batches = [(batch, this_res) for batch, (this_res, success) in zip(subfolders_to_test, res) if not success]
                    print(df_incomplete.groupby(['subpath_id','status'])['shard'].count().to_string())

                    df_fail = df_incomplete[df_incomplete.status == 'FAIL']
                    df_fail = df_fail.merge(df_sample, left_on=['subpath_id', 'shard'], right_on=['batch', 'shard'], how='left')
                    df_fail.to_csv(f'{args.output}.failure.tsv', sep='\t', index=False)
                    if df_fail.shape[0] > 0:
                        print('')
                        print(f'{str(df_fail.shape[0])} have failed across the following stages:')
                        print(df_fail[['subpath_id','shard', 'latest_task']].to_string())
                    
                    df_running = df_incomplete[df_incomplete.status == 'IN PROGRESS']
                    df_running.to_csv(f'{args.output}.running.tsv', sep='\t', index=False)
                    if df_running.shape[0] > 0:
                        print('')
                        print(f'{str(df_running.shape[0])} are running:')
                        print(df_running[['subpath_id', 'shard', 'latest_task']].to_string())

                    if not args.hide_pass:
                        df_stalled = df_incomplete[df_incomplete.status == 'PASS']
                        df_stalled.to_csv(f'{args.output}.stalled.tsv', sep='\t', index=False)
                        if df_stalled.shape[0] > 0:
                            print('')
                            print(f'{str(df_stalled.shape[0])} are listed as PASS (but are in incomplete batches):')
                            print(df_stalled[['subpath_id', 'shard', 'latest_task']].to_string())

    if args.get_shard_count:
        ct = {x: count_shards_single(storage_client, bucket, x, args.run_folder) for x in subfolders_to_test}
        for fold, shardcount in ct.items():
            print(f'For batch {fold}:')
            print(f'Shard count: {str(shardcount)}')

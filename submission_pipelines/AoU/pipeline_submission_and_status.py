import copy
from datetime import datetime, timedelta
import dateutil
import gcsfs
from google.cloud import storage
import json
import numpy
import os
import pandas
from pathlib import Path
import requests
import subprocess
import sys
import time
import traceback
from zoneinfo import ZoneInfo

#local module import
import cromwell_run_monitor as crm

def submit_cromwell_workflows(samples_df, run_name='test_pipeline2', batch_size=10, 
                              mtSwirl_root='/home/jupyter/mtSwirl_fork/mtSwirl/', 
                              num_to_submit=None, num_concurrent_crams=None,
                              addl_sub_interval_sec=0, cromwell_timeout=60, 
                              submission_retries=2):
    '''Takes a dataframe based on the srWGS manifest.csv file, and submits the CRAMs in that file for
    processing by the mtSwirl pipeline on Cromwell. Creates the inputs json file for each batch of 
    batch_size, saves the inputs to the workspace bucket in a directory called run_name,
    submits the batches to cromwell, and saves an expanded version of samples_df to the same directory
    that records the cromwell workflow IDs, status of each run, and the directory in the cloud bucket
    where to find the inputs for the job processing each sample. Returns the gsutil URI to that 
    modified samples_df, which in turn can be used as input to update_cromwell_status() to both update
    the run_status attributes of all submitted samples and return the updated samples_df.
    
    Required column names in samples_df:
    person_id: stores the id number that uniquely identifies each sample in the All of Us project (MUST BE STR DTYPE)
    cram_uri: stores the URI of the srWGS CRAM file corresponding to each person_id
    cram_index_uri: stores the URI of the index of the srWGS CRAM file corrsponding to each person_id
    
    Other parameters:
    run_name: A name for this set of pipeline runs. This will be the root directory for the pipeline inputs
              and results in the workspace bucket.
    batch_size: How many CRAMs should be grouped together into each Cromwell workflow. The mtSwirl pipeline
                has a final step that aggregates the stats and vcf URIs for the samples in that batch, so 
                this step controls how much work this intermediate merging step does. An important consideration
                is that if a batch workflow fails, the pipeline must be re-run for all samples in that batch.
                So, this shouldn't be too large.
    mtSwirl_root: The location on the local file system where the mtSwirl GitHub has been cloned.
    num_to_submit: An integer value giving the maximum number of CRAMs to submit from samples_df. 
                   The minimum value between this parameter and num_to_submit will control how many 
                   workflows are submitted to Cromwell during this function call.
    num_concurrent_crams: An integer value giving the maximum number of CRAMs that should be processed in 
                          parallel. The minimum value between this parameter and num_to_submit will control
                          how many workflows are submitted to Cromwell during this function call.
    addl_sub_interval_sec: Number of seconds to wait in between submission of workflows, passed to time.sleep(). 
                           This can be helpful to avoid synchrony in the running jobs so that you don't have
                           many jobs reaching the same point in a pipeline all at once and potentially
                           creating a bottleneck with GCP resource quotas.
    '''
    #set up run locations in the cloud
    bucket = os.getenv("WORKSPACE_BUCKET")
    project = os.getenv("GOOGLE_PROJECT")
    cloud_fs = gcsfs.GCSFileSystem(project=project, requester_pays=True)

    # set up directory to keep track of the inputs generated during this submission
    output_bucket = os.path.join(bucket, run_name)
    samples_df_uri = os.path.join(output_bucket, 'samples_df.csv')

    #get template of pipeline JSON input
    pipeline_wdl = os.path.join(mtSwirl_root, 'WDL/v2.5_MongoSwirl_Single/scatterWrapper_MitoPipeline_v2_5.wdl')
    json_template_path = os.path.join(mtSwirl_root, 'WDL/files/prepopulated_inputs.json')
    with open(json_template_path) as json_in:
        json_template = json.load(json_in)

    json_template['MitochondriaPipelineWrapper.requester_pays_project'] = project

    #upload to GCS and add the auxiliary workflow script locations to the json_template
    aux_script_params = {
        "MitochondriaPipelineWrapper.nuc_interval_list": os.path.join(mtSwirl_root, 'WDL/files/NUMTv3_all385.hg38.interval_list'),
        "MitochondriaPipelineWrapper.CheckVariantBoundsScript": os.path.join(mtSwirl_root, 'WDL/scripts/check_variant_bounds.R'),
        "MitochondriaPipelineWrapper.FaRenamingScript": os.path.join(mtSwirl_root, 'WDL/scripts/compatibilify_fa_intervals_consensus.R'),
        "MitochondriaPipelineWrapper.CheckHomOverlapScript": os.path.join(mtSwirl_root, 'WDL/scripts/check_overlapping_homoplasmies.R'),
        "MitochondriaPipelineWrapper.HailLiftover": os.path.join(mtSwirl_root, 'WDL/scripts/fix_liftover.py'),
        "MitochondriaPipelineWrapper.JsonTools": os.path.join(mtSwirl_root, 'WDL/scripts/jsontools.py'),
        "MitochondriaPipelineWrapper.MergePerBatch": os.path.join(mtSwirl_root, 'WDL/scripts/merge_per_batch.py')}
    for k, aux_path in aux_script_params.items():
        aux_uri = os.path.join(output_bucket, '/'.join(aux_path.split('/')[-2:]))
        subprocess.run(['gsutil', '-q', 'cp', aux_path, aux_uri], check=True)
        aux_script_params[k] = aux_uri
    json_template.update(aux_script_params)

    #initialize the workflow tracking columns if they do not yet exist in the submitted samples_df
    if numpy.sum(samples_df.columns == 'cromwell_id') == 0:
        samples_df['cromwell_id'] = 'Submission pending'
    if numpy.sum(samples_df.columns == 'inputs_uri') == 0:
        samples_df['inputs_uri'] = 'Submission pending'
    if numpy.sum(samples_df.columns == 'run_status') == 0:
        samples_df['run_status'] = 'Submission pending'
    
    #subset to just the samples with a status of "Submission pending" for batching and submission
    to_submit_df = samples_df.loc[samples_df['run_status'] == 'Submission pending'].copy()
    if to_submit_df.shape[0] == 0:
        print('No samples with a a status of "Submission pending". Nothing new to submit.')
        return samples_df_uri
    
    #check the number of running samples
    num_running = numpy.sum(samples_df['run_status'] == 'Running')

    #group the samples ready for submission into batches of size batch_size, record the inputs, and do the submission
    submission_records_root = os.path.join(output_bucket, 'cromwell_submissions')

    #compute offset for the batch ID (so that inputs for previous batches are not overwritten)
    submitted_batch_list = samples_df.loc[samples_df['run_status'] != 'Submission pending', 'inputs_uri'].unique()
    if len(submitted_batch_list) > 0:
        batch_idx_offset = max([int(os.path.basename(elt).split('_')[-1]) for elt in submitted_batch_list]) + 1
    else:
        batch_idx_offset = 0

    #compute the ending value of batch idx range
    if num_to_submit is None and num_concurrent_crams is None:
        range_end = to_submit_df.shape[0]
    elif num_to_submit is None:
        range_end = min(to_submit_df.shape[0], num_concurrent_crams - num_running)
    elif num_concurrent_crams is None:
        range_end = min(to_submit_df.shape[0], num_to_submit)
    else:
        range_end = min(to_submit_df.shape[0], min(num_to_submit, num_concurrent_crams - num_running))
    range_end = max(range_end, 0)
    print(f'Found {to_submit_df.shape[0]} samples awaiting submission, and {num_running} '
          'currently running samples. Given the provided submission and concurrency limits, '
          f'will submit {range_end} of them in batches of {batch_size}, starting from batch {batch_idx_offset}.')

    for idx in range(0, range_end, batch_size):
        #extract batch samples
        batch_df = to_submit_df.iloc[idx:idx+batch_size]

        #create the lists of samples to use as input to the mtSwirl job
        batch_root = os.path.join(submission_records_root, f'batch_{(idx//batch_size) + batch_idx_offset}')
        print(batch_root)

        sample_list_uri = os.path.join(batch_root, 'sample_name_list.txt')
        with cloud_fs.open(sample_list_uri, 'w') as out:
            out.write('\n'.join(batch_df['person_id'].to_list()) + '\n')

        cram_file_list_uri = os.path.join(batch_root, 'cram_uri_list.txt')
        with cloud_fs.open(cram_file_list_uri, 'w') as out:
            out.write('\n'.join(batch_df['cram_uri'].to_list()) + '\n')

        crai_file_list_uri = os.path.join(batch_root, 'crai_uri_list.txt')
        with cloud_fs.open(crai_file_list_uri, 'w') as out:
            out.write('\n'.join(batch_df['cram_index_uri'].to_list()) + '\n')

        #update the input parameters json with the locations of the batch sample lists
        batch_json_update = {'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_list': cram_file_list_uri,
                             'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_index_list': crai_file_list_uri,
                             'MitochondriaPipelineWrapper.sample_name_list': sample_list_uri,
                             'MitochondriaPipelineWrapper.force_manual_download': 'false'}
        batch_json = copy.deepcopy(json_template)
        batch_json.update(batch_json_update)

        #save the input parameters to both the cloud and the local directory (cromshell submit needs a local json file)
        batch_input_json_uri = os.path.join(batch_root, 'input_parameters.json')
        with cloud_fs.open(batch_input_json_uri, 'w') as out:
            json.dump(batch_json, out)

        batch_input_json = './input_parameters.json'
        with open(batch_input_json, 'w') as out:
            json.dump(batch_json, out)

        #submit the batch as a cromwell workflow
        cmd_arg_list = ['cromshell', '-t', str(cromwell_timeout), '--no_turtle', '--machine_processable', 'submit', pipeline_wdl, batch_input_json]
        batch_submission_cmd_uri = os.path.join(batch_root, 'batch_submission_cmd.txt')
        with cloud_fs.open(batch_submission_cmd_uri, 'w') as out:
            out.write(' '.join(cmd_arg_list) + '\n')
        while True:
            try:
                sub_resp = subprocess.run(cmd_arg_list, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                if submission_retries == 0:
                    raise
                submission_retries -= 1
            else:
                break
        sub_info = json.loads(sub_resp.stdout.decode())
        sub_id_txt_uri = os.path.join(batch_root, 'sub_id.txt')
        with cloud_fs.open(sub_id_txt_uri, 'w') as out:
            out.write(sub_info['id'] + '\n')

        #write the info about the submitted jobs back to the full samples_df dataframe
        samples_df.loc[batch_df.index, 'cromwell_id'] = sub_info['id']
        samples_df.loc[batch_df.index, 'run_status'] = sub_info['status']
        samples_df.loc[batch_df.index, 'inputs_uri'] = batch_root
        samples_df.to_csv(samples_df_uri, storage_options={'project':project, 'requester_pays':True})
        print(f'Batch {idx//batch_size + batch_idx_offset} submitted.')
        
        #wait between batch submissions, if requested
        time.sleep(addl_sub_interval_sec)
    print('Submission complete.')
    return samples_df_uri

def get_app_details(env, app_name):
    '''From the AoU Researcher Workbench Cromwell initiation code snippet.
    Collects metadata about a running Cromwell instance on Researcher Workbench.
    '''
    get_app_url = f'{env["leonardo_url"]}/api/google/v1/apps/{env["google_project"]}/{app_name}'
#    print('start')
    r = requests.get(
        get_app_url,
        params={
            'includeDeleted': 'true',
            'role': 'creator'
        },
        headers={
            'Authorization': f'Bearer {env["token"]}'
        }
    )
    if r.status_code == 404:
        return 'DELETED', None, None, None
    else:
        r.raise_for_status()
    result_json = r.json()
    custom_environment_variables = result_json['customEnvironmentVariables']
    return result_json['status'], custom_environment_variables['WORKSPACE_NAMESPACE'], result_json.get('proxyUrls')

def check_for_app(env):
    '''From the AoU Researcher Workbench Cromwell initiation code snippet.
    Checks for the running Cromwell instance on Researcher Workbench.
    '''
    list_apps_url = f'{env["leonardo_url"]}/api/google/v1/apps/{env["google_project"]}'
    r = requests.get(
        list_apps_url,
        params={
          'includeDeleted': 'false'
        },
        headers = {
            'Authorization': f'Bearer {env["token"]}'
        }
    )
    r.raise_for_status()

    for potential_app in r.json():
        if potential_app['appType'] == 'CROMWELL' and (
                str(potential_app['auditInfo']['creator']) == env['owner_email']
                or str(potential_app['auditInfo']['creator']) == env['user_email']
        ) :
            potential_app_name = potential_app['appName']
            potential_app_status = potential_app['status']

            # We found a CROMWELL app in the correct google project and owned by the user. Now just check the workspace:
            _, workspace_namespace,  proxy_url = get_app_details(env, potential_app_name)
            if workspace_namespace == env['workspace_namespace']:
                return potential_app_name, potential_app_status, proxy_url['cromwell-service']

    return None, None, None

def get_cromwell_url():
    env = {
        'workspace_namespace': os.environ['WORKSPACE_NAMESPACE'],
        'workspace_bucket': os.environ['WORKSPACE_BUCKET'],
        'user_email': os.environ.get('PET_SA_EMAIL', default = os.environ['OWNER_EMAIL']),
        'owner_email': os.environ['OWNER_EMAIL'],
        'google_project': os.environ['GOOGLE_PROJECT'],
        'leonardo_url': os.environ['LEONARDO_BASE_URL']
    }

    # Fetch the token:
    token_fetch_command = subprocess.run(['gcloud', 'auth', 'print-access-token', env['user_email']], capture_output=True, check=True, encoding='utf-8')
    env['token'] = str.strip(token_fetch_command.stdout)
    
    #get cromwell app info
    app_name, app_status, proxy_url = check_for_app(env)
    if app_status != 'RUNNING':
        raise Exception('Error: Cromwell doesn\'t appear to be running.')
    return app_name, app_status, proxy_url, env

def assemble_status_urls(workflow_ids):
    app_name, app_status, app_url, env = get_cromwell_url()
    workflow_status_urls = [app_url + f'/api/workflows/v1/{elt}/status' for elt in workflow_ids]
    return workflow_status_urls, env['token']

def check_status(workflow_status_url, token):
    r = requests.get(
        workflow_status_url,
        timeout=5,
        verify=True,
        headers={
            'Referer': 'https://notebooks.firecloud.org',
            'Authorization': f'Bearer {token}'
        }
    )
    r.raise_for_status()
    return r.json()

def update_cromwell_status(sample_df_uri, verbose=True, quiet=False, cromwell_id_col='cromwell_id',
                           status_col='run_status', no_write=False):
    if quiet is True:
        verbose = False
    sample_df = pandas.read_csv(sample_df_uri, index_col=0,
                                storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    to_check = sample_df.loc[~sample_df[status_col].isin(['Submission pending', 'Succeeded', 'Failed', 'Manual review']), cromwell_id_col].unique()
    if not quiet:
        print(f'Updating status for {to_check.shape[0]} workflow ids.')
    workflow_status_urls, token = assemble_status_urls(to_check)
    for idx, (crom_id, crom_url) in enumerate(zip(to_check, workflow_status_urls)):
        if (verbose is True) and idx and not idx%20:
            print(f'{idx} workflows checked.')
        if verbose is True:
            print(crom_id)
        status_json = check_status(crom_url, token)
        sample_df.loc[sample_df[cromwell_id_col] == crom_id, status_col] = status_json['status']
    if not no_write:
        sample_df.to_csv(sample_df_uri,
                         storage_options={'project':os.getenv('GOOGLE_PROJECT'),
                                          'requester_pays':True})
    if verbose is True:
        print('\n')
    if not quiet:
        print(sample_df[status_col].value_counts().to_string())
    return sample_df

def get_detailed_workflow_status(workflow_id, cromwell_run_prefix='cromwell-execution'):
    bucket_id = os.path.basename(os.getenv("WORKSPACE_BUCKET"))
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_id)
#    cromwell_run_prefix = 'cromwell-execution'
#    workflow_subfolder = 'da52343a-976f-4bf5-a60d-e20683a40581'

    return crm.check_success_single(storage_client, bucket, workflow_id, cromwell_run_prefix, False)

def resubmit_failed_workflow(samples_df_uri, cromwell_id_to_resubmit):
#    cromwell_id_to_resubmit = 'da52343a-976f-4bf5-a60d-e20683a40581'

    #get the samples_df to update
#    samples_df_uri = 'gs://fc-secure-34d99fb0-3748-4367-8197-b01069a4a7f9/medium_test_pipeline/samples_df.csv'
    samples_df = pandas.read_csv(samples_df_uri, index_col=0, 
                                 storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})

    #get the batch corresponding to this cromwell ID
    batch_uri = samples_df.loc[samples_df['cromwell_id'] == cromwell_id_to_resubmit, 'inputs_uri'].values[0]
    print(f'Resubmitting job for {os.path.basename(batch_uri)} with cromwell ID {cromwell_id_to_resubmit}.')

    #set up the cloud_fs for writing input files
    bucket = os.getenv("WORKSPACE_BUCKET")
    project = os.getenv("GOOGLE_PROJECT")
    cloud_fs = gcsfs.GCSFileSystem(project=project, requester_pays=True)

    #download the input_parameters.json for the batch we want to resubmit
    batch_inputs_uri = os.path.join(batch_uri, 'input_parameters.json')
    subprocess.run(['gsutil', '-q', '-u', project, 'cp', batch_inputs_uri, './input_parameters.json'], check=True)

    #preserve the old sub_id.txt in preparation for resubmission and making a new one
    prev_sub_id_uri_glob = os.path.join(batch_uri, 'prev_sub_id*')
    try:
        res = subprocess.run(['gsutil', '-u', project, 'ls', prev_sub_id_uri_glob], check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        if err.stderr.decode().strip() == 'CommandException: One or more URLs matched no objects.':
            prev_sub_idx = 0
        else:
            raise err
    else:
        prev_sub_id_list = res.stdout.decode().strip().split('\n')
        prev_sub_idx = len(prev_sub_id_list)
    prev_sub_id_uri = os.path.join(batch_uri, f'prev_sub_id{prev_sub_idx}.txt')
    sub_id_uri = os.path.join(batch_uri, 'sub_id.txt')
    subprocess.run(['gsutil', '-q', '-u', project, 'mv', sub_id_uri, prev_sub_id_uri], check=True)

    #get the submission command to use
    batch_cmd_uri = os.path.join(batch_uri, 'batch_submission_cmd.txt')
    with cloud_fs.open(batch_cmd_uri) as cmd_in:
        cmd_txt = cmd_in.read().decode().strip()

    #submit the batch as a cromwell workflow
    cmd_arg_list = cmd_txt.split()
    sub_resp = subprocess.run(cmd_arg_list, check=True, capture_output=True)
    sub_info = json.loads(sub_resp.stdout.decode())
    sub_id_txt_uri = os.path.join(batch_uri, 'sub_id.txt')
    with cloud_fs.open(sub_id_txt_uri, 'w') as out:
        out.write(sub_info['id'] + '\n')

    #write the info about the submitted jobs back to the full samples_df dataframe
    samples_df.loc[samples_df['cromwell_id'] == cromwell_id_to_resubmit, 'run_status'] = sub_info['status']
    samples_df.loc[samples_df['cromwell_id'] == cromwell_id_to_resubmit, 'cromwell_id'] = sub_info['id']
    samples_df.to_csv(samples_df_uri, storage_options={'project':project, 'requester_pays':True})
    print(f'Batch {os.path.basename(batch_uri)} resubmitted with workflow ID {sub_info["id"]}.')
    return samples_df

def isolate_failed_shards(samples_df_uri, failed_workflows_list, cromwell_run_prefix='cromwell-execution'):
    #get samples_df
    samples_df = pandas.read_csv(samples_df_uri, index_col=0,
                                 storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    
    # iterate over workflow IDs and change any samples without a FAIL status in the Cromwell run
    # to "Submission pending" status in the samples_df so that they can be packaged into a new
    # batch and re-submitted
    for crom_id in failed_workflows_list:
        #check that this is a valid workflow id
        if crom_id not in samples_df['cromwell_id'].values:
            print(f'Workflow ID {crom_id} is not found. Skipping.')
            continue

        #get relevant subset of the samples_df
        workflow_samples = samples_df.loc[samples_df['cromwell_id'] == crom_id].copy()

        #check to be sure these are failed samples
        workflow_status = workflow_samples['run_status'].unique()
        if len(workflow_status) != 1 or workflow_status[0] != 'Failed':
            print(f'Workflow {crom_id} has a status other than "Failed" (status: {", ".join(workflow_status)}). Skipping.')
            continue

        #get sample submission order for the requested workflow
        samples_order_uri = os.path.join(workflow_samples['inputs_uri'].unique()[0], 'sample_name_list.txt')
        samples_order = pandas.read_csv(samples_order_uri, header=None, 
                                        storage_options={'project':os.getenv('GOOGLE_PROJECT'),
                                                         'requester_pays':True})
        
        #get the table of metadata about the shards of this workflow
        run_info = get_detailed_workflow_status(crom_id, cromwell_run_prefix=cromwell_run_prefix)[0]
        
        #identify the failed shards and update the status of the non-failed ones so they will 
        # be resubmitted in a new batch
        failed_idx = [int(elt.split('-')[-1]) for elt in run_info.loc[run_info['status'] == 'FAIL', 'shard'].values]
        to_resub_idx = workflow_samples.loc[~workflow_samples['person_id'].isin(samples_order.iloc[failed_idx, 0].to_list())].index
        samples_df.loc[to_resub_idx, 'run_status'] = 'Submission pending'
        to_rev_idx = workflow_samples.loc[workflow_samples['person_id'].isin(samples_order.iloc[failed_idx, 0].to_list())].index
        samples_df.loc[to_rev_idx, 'run_status'] = 'Manual review'

    #save the results back to the cloud
    samples_df.to_csv(samples_df_uri, storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    return samples_df

def get_succeeded_job_metrics(samples_df_uri, run_name, force_reload=False):
    run_metrics_uri = os.path.join(os.getenv("WORKSPACE_BUCKET"), f'{run_name}/run_metrics.csv')
    run_metrics = {'cromwell_id':[],
                   'status':[],
                   'start_time':[],
                   'end_time':[],
                   'runtime':[],
                   'merging_log':[],
                   'merged_calls':[],
                   'merged_coverage':[],
                   'merged_statistics':[]}
    #get samples_df
    if isinstance(samples_df_uri, str):
        samples_df = pandas.read_csv(samples_df_uri, index_col=0,
                                     storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
        workflow_ids = samples_df.loc[samples_df['run_status'] == 'Succeeded', 'cromwell_id'].unique()
    else:
        workflow_ids = samples_df_uri.loc[samples_df_uri['run_status'] == 'Succeeded', 'cromwell_id'].unique()
    #get run_metrics (if it already exists), otherwise just make it an empty version of the one we are building
    if not force_reload:
        try:
            existing_run_metrics = pandas.read_csv(run_metrics_uri, index_col=0, sep='\t',
                                                   storage_options={'project':os.getenv('GOOGLE_PROJECT'), 
                                                                    'requester_pays':True})
        except:
            existing_run_metrics = pandas.DataFrame(run_metrics)
        else:
            #if there is an existing run_metrics file, filter for just the workflow_ids that haven't been 
            # recorded there yet.
            workflow_ids = [elt for elt in workflow_ids if elt not in existing_run_metrics['cromwell_id'].values]
    else:
        existing_run_metrics = pandas.DataFrame(run_metrics)

    # now get the results for all of the newly-succeeded workflows
    for idx, workflow_id in enumerate(workflow_ids):
        if idx and not idx%5:
            print(f'Processed {idx} workflows')
        attempts = 4
        to_raise = None
        while attempts >= 0:
            try:
                run_meta_resp = subprocess.run(['cromshell', '-t', '60', '--no_turtle', '--machine_processable', 'metadata', 
                                                '--dont-expand-subworkflows', workflow_id], 
                                                check=True, capture_output=True)
            except Exception as err:
                print(f'Error retrieving info about workflow {workflow_id}. Retrying {attempts} more times.')
                to_raise = err
                attempts -= 1
                time.sleep(15)
            else:
                break
        else:
            raise to_raise
        run_meta = json.loads(run_meta_resp.stdout.decode())
        run_metrics['cromwell_id'].append(workflow_id)
        run_metrics['status'].append(run_meta['status'])
        run_metrics['start_time'].append(run_meta['start'])
        run_metrics['end_time'].append(run_meta['end'])
        run_metrics['runtime'].append(str(dateutil.parser.isoparse(run_meta['end'])
                                          - dateutil.parser.isoparse(run_meta['start'])))
        try:
            run_metrics['merging_log'].append(run_meta['calls']['MitochondriaPipelineWrapper.MergeMitoMultiSampleOutputsInternal'][-1]['backendLogs']['log'])
        except KeyError:
            run_metrics['merging_log'].append('Not found')
        run_metrics['merged_calls'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_calls', 'Not found'))
        run_metrics['merged_coverage'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_coverage', 'Not found'))
        run_metrics['merged_statistics'].append(run_meta['outputs'].get('MitochondriaPipelineWrapper.merged_statistics', 'Not found'))

    run_metrics = pandas.concat([existing_run_metrics, pandas.DataFrame(run_metrics)]).reset_index(drop=True)
    run_metrics.to_csv('./run_metrics.csv', sep='\t')
    run_metrics.to_csv(run_metrics_uri, sep='\t',
                       storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})
    print(run_metrics_uri)
    return run_metrics

def write_pipeline_params_json():
    # set pipeline submission parameter values
    run_name = 'v6_test_1000'
    batch_size = 25      # The number of samples to be processed per Cromwell workflow
    num_concurrent = 4000 # The max number of samples that can be processed concurrently (i.e. to stay within GCP quotas)
    min_fail_num = 5      # Allow this many failures before worrying about the percentage failure
    max_fail_pct = 15.0   # Halt the submission loop if more than this percentage of batches fail
    submission_wait_sec = 30  # Wait seconds in between workflow submissions to smooth out the GCP resource bottleneck
    sub_check_wait_min = 2 # Number of minutes to wait in between checking for whether jobs need to be submitted
    pipeline_status_path = './v6_test_1000.pipeline_status.tsv'
    pipeline_status_uri = os.path.join(os.getenv('WORKSPACE_BUCKET'), run_name, 'v6_test_1000.pipeline_status.tsv')

    #make cloud FS instance to write to GCS
    cloud_fs = gcsfs.GCSFileSystem(project=os.getenv("GOOGLE_PROJECT"), requester_pays=True)

    # save variables to a json file so that they can be changed mid-run if necessary
    pipe_params = {'RUN_NAME':run_name,
                   'BATCH_SIZE':batch_size,
                   'NUM_CONCURRENT':num_concurrent,
                   'MIN_FAIL_NUM':min_fail_num,
                   'MAX_FAIL_PCT':max_fail_pct,
                   'SUBMISSION_WAIT':submission_wait_sec,
                   'SUB_CHECK_WAIT':sub_check_wait_min,
                   'PIPELINE_STATUS_PATH':pipeline_status_path,
                   'PIPELINE_STATUS_URI':pipeline_status_uri}
    #save to local fs
    params_path = './pipeline_submission_params.json'
    with open(params_path, 'w') as out:
        json.dump(pipe_params, out)
    #save to cloud
    params_uri = os.path.join(os.getenv('WORKSPACE_BUCKET'), run_name, 'pipeline_submission_params.json')
    with cloud_fs.open(params_uri, 'w') as out:
        json.dump(pipe_params, out)
    print(f'To change parameter values for the running submission process, edit the file here:\n{params_uri}')
    return

def run_pipeline():
    #Now, submit all workflows until 40k samples have terminated
    pipeline_status = {'timestamp':[],
                       'Submission pending':[],
                       'Submitted':[],
                       'Running':[],
                       'Succeeded':[],
                       'Failed':[]}

    # submit the initial set of jobs and then periodically check to submit more as they complete, up to 40k samples
    try:
        while numpy.sum(test_run['run_status'] != 'Submission pending') < 1000:
            #reload submission parameters each iteration to allow them to be tuned over the course of the run
            with cloud_fs.open(params_uri, 'r') as params_in:
                params_json = json.load(params_in)
            run_name = params_json['RUN_NAME']
            batch_size = params_json['BATCH_SIZE']
            num_concurrent = params_json['NUM_CONCURRENT']
            min_fail_num = params_json['MIN_FAIL_NUM']
            max_fail_pct = params_json['MAX_FAIL_PCT']
            submission_wait_sec = params_json['SUBMISSION_WAIT']
            sub_check_wait_min = params_json['SUB_CHECK_WAIT']
            pipeline_status_path = params_json['PIPELINE_STATUS_PATH']
            pipeline_status_uri = params_json['PIPELINE_STATUS_URI']

            #Now, do the submission
            test_run_uri = pss.submit_cromwell_workflows(test_run, run_name=run_name, batch_size=batch_size, 
                                                         mtSwirl_root='/home/jupyter/mtSwirl_fork/mtSwirl/',
                                                         num_concurrent_crams=num_concurrent,
                                                         addl_sub_interval_sec=submission_wait_sec)

            time_now = datetime.now(tz=ZoneInfo('America/New_York')).isoformat()
            print(time_now)

            # update and record the summary of the pipeline status
            status_counts = test_run['run_status'].value_counts()
            status_counts_dict = status_counts.to_dict()
            pipeline_status['timestamp'].append(time_now)
            for k in sorted(set(pipeline_status.keys()) | set(status_counts_dict.keys())):
                if k == 'timestamp':
                    continue
                try:
                    pipeline_status[k].append(status_counts_dict.get(k, 0))
                except KeyError:
                    #any less-frequent pipeline status types will be added to the dataframe on the fly
                    pipeline_status[k] = ([0]*(len(pipeline_status['timestamp'])-1)) + [status_counts_dict[k]]
            pipeline_status_df = pandas.DataFrame(pipeline_status)
            pipeline_status_df.to_csv(pipeline_status_path, sep='\t', index=False)
            pipeline_status_df.to_csv(pipeline_status_uri, sep='\t', index=False,
                                      storage_options={'project':os.getenv('GOOGLE_PROJECT'), 'requester_pays':True})

            # wait before updating the job status and attempting to submit more jobs
            time.sleep(60*sub_check_wait_min)

            # update the status results
            test_run = pss.update_cromwell_status(test_run_uri, verbose=False)
            test_run['person_id'] = test_run['person_id'].astype(str)

            # check whether we are having excessive failures (reload status_counts to get updated numbers)
            status_counts = test_run['run_status'].value_counts()
            status_counts_dict = status_counts.to_dict()
            success_count = status_counts_dict.get('Succeeded', 0)
            failure_count = status_counts_dict.get('Failed', 0)
            pct_failed = (failure_count/(failure_count+success_count))*100

            # if so, stop the loop
            if (failure_count > min_fail_num) and (pct_failed > max_fail_pct):
                msg = (f'{pct_failed:0.1f}% of terminated jobs have ended in failure, which is greater\n'
                       f'than the threshold setting of {max_fail_pct}%. Halting the submission loop.\n'
                       'Please check to see if something is wrong.')
                with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
                    out.write(msg + '\n')
                print(msg)
                break
        else:
            print('All samples submitted. Submission loop ending.')
    except Exception:
        excpt_str = traceback.format_exc()
        print(excpt_str)
        with open('./PIPELINE_SUBMISSION_ERROR', 'w') as out:
            out.write(excpt_str + '\n')
        raise

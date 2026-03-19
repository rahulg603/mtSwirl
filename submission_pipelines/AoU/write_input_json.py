#! /usr/bin/env python

from copy import deepcopy
import gcsfs
import json
import os
import pandas as pd
import subprocess

#### PARAMETERS
n_test = 1 # number of samples per iteration
n_iter = 1 # number of iterations
JOBLIMIT = 5000 # how many jobs to allow simulataneously
path_completed_samples = "mt_pipeline_single_2_5_stats.tsv" # expects this to be in the current bucket
path_failed_samples = "mt_pipeline_single_2_5_failures.tsv" # samples that have failed along with log files
tf_rerun_fail = False
tf_force_dl = False

# assumes the structure of the repo is .../mtSwirl/generate_mtdna_call_mt/AoU/write_input_json.py
package_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

bucket = os.getenv("WORKSPACE_BUCKET")
project = os.getenv("GOOGLE_PROJECT")
cloud_fs = gcsfs.GCSFileSystem(project=project, requester_pays=True)
# set up cromwell directories
path_indiv_save = 'test_output_new_printreads'
output_bucket = os.path.join(bucket, path_indiv_save)  # This is where the output of the WDL will be.

options_filename = "options.json"

if __name__ == "__main__":
    #get template of pipeline JSON input
    json_template_path = os.path.join(package_root, 'WDL/files/prepopulated_inputs.json')
    with open(json_template_path) as json_in:
        json_template = json.load(json_in)

    json_template['MitochondriaPipelineWrapper.requester_pays_project'] = project

    #upload to GCS and add the auxiliary workflow script locations to the json_template
    aux_script_params = {
        "MitochondriaPipelineWrapper.nuc_interval_list": os.path.join(package_root, 'WDL/files/NUMTv3_all385.hg38.interval_list'),
        "MitochondriaPipelineWrapper.CheckVariantBoundsScript": os.path.join(package_root, 'WDL/scripts/check_variant_bounds.R'),
        "MitochondriaPipelineWrapper.FaRenamingScript": os.path.join(package_root, 'WDL/scripts/compatibilify_fa_intervals_consensus.R'),
        "MitochondriaPipelineWrapper.CheckHomOverlapScript": os.path.join(package_root, 'WDL/scripts/check_overlapping_homoplasmies.R'),
        "MitochondriaPipelineWrapper.HailLiftover": os.path.join(package_root, 'WDL/scripts/fix_liftover.py'),
        "MitochondriaPipelineWrapper.JsonTools": os.path.join(package_root, 'WDL/scripts/jsontools.py'),
        "MitochondriaPipelineWrapper.MergePerBatch": os.path.join(package_root, 'WDL/scripts/merge_per_batch.py')}
    for k, aux_path in aux_script_params.items():
        aux_uri = os.path.join(output_bucket, '/'.join(aux_path.split('/')[-2:]))
        subprocess.run(['gsutil', 'cp', aux_path, aux_uri], check=True)
        aux_script_params[k] = aux_uri
    json_template.update(aux_script_params)

    # grab list of completed samples
    comp_samps_uri = os.path.join(bucket, path_completed_samples)
    try:
        with cloud_fs.open(comp_samps_uri) as comp_samps_in:
            df_stats = pd.read_csv(comp_samps_in, sep='\t')
            sample_list_completed = list(df_stats.s)
    except FileNotFoundError:
        sample_list_completed = []

    # grab list of failed samples
    fail_samps_uri = os.path.join(bucket, path_failed_samples)
    try:
        with cloud_fs.open(fail_samps_uri) as fail_samps_in:
            df_failed = pd.read_csv(fail_samps_in, sep='\t')
        sample_list_failed = list(df_failed.s)
    except FileNotFoundError:
        sample_list_failed = []

    # get all cram paths (AoU manifest)
    manifest_uri = os.path.join(os.getenv('CDR_STORAGE_PATH'), 'wgs/cram/manifest.csv')
    with cloud_fs.open(manifest_uri) as manifest_in:
        manifest_df = pd.read_csv(manifest_in)

    #filter manifest for samples that have already been processed or failed
    original_size = manifest_df.shape[0]
    manifest_df = manifest_df[~manifest_df.person_id.isin(sample_list_completed)]
    if tf_rerun_fail:
        manifest_df = manifest_df[manifest_df.person_id.isin(sample_list_failed)]
        new_size = manifest_df.shape[0]
        if not all(x in [x in list(manifest_df.person_id) for x in sample_list_failed]):
            raise ValueError("ERROR: All failed samples should be found in the sample list.")
        print(f'{str(new_size)} previously failed samples added to list for rerunning.')
    else:
        manifest_df = manifest_df[~manifest_df.person_id.isin(sample_list_failed)]
        new_size = manifest_df.shape[0]
        if original_size - new_size != len(sample_list_completed) + len(sample_list_failed):
            raise ValueError("ERROR: the number of samples removed from manifest must be the same as the number of processed samples + number of prior failures.")
        print(f'{str(len(sample_list_completed) + len(sample_list_failed))} samples already processed or failed and removed from manifest.')
        print(f'{str(new_size)} samples remain in the manifest (of {str(original_size)} originally present).')
    manifest_df = manifest_df.reset_index(drop=True)

    
    #Generate the JSON input files for the batches of samples that should be processed together
    json_collection = []
    max_rows = manifest_df.shape[0]

    for idx in range(0, n_iter):
        if (idx + 1) * n_test >= max_rows:
            this_max = max_rows
            break_here = True
        else:
            this_max = (idx + 1) * n_test
            break_here = False

        manifest_df_sub = manifest_df.iloc[(idx * n_test):this_max]
        s = list(manifest_df_sub.person_id.astype(str))
        cram_paths = list(manifest_df_sub.cram_uri)
        crai_paths = list(manifest_df_sub.cram_index_uri)

        sample_list_uri = os.path.join(output_bucket, f'sample_list{idx!s}.txt')
        with cloud_fs.open(sample_list_uri, "w") as f:
            f.write('\n'.join(s) + '\n')

        cram_file_list_uri = os.path.join(output_bucket, f'cram_file_list{idx!s}.txt')
        with cloud_fs.open(cram_file_list_uri, "w") as f:
            f.write('\n'.join(cram_paths) + '\n')

        crai_file_list_uri = os.path.join(output_bucket, f'crai_file_list{idx!s}.txt')
        with cloud_fs.open(crai_file_list_uri, "w") as f:
            f.write('\n'.join(crai_paths) + '\n')

        dct_update = {'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_list': cram_file_list_uri,
                      'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_index_list': crai_file_list_uri,
                      'MitochondriaPipelineWrapper.sample_name_list': sample_list_uri,
                      'MitochondriaPipelineWrapper.force_manual_download': tf_force_dl}

        this_json = deepcopy(json_template)
        this_json.update(dct_update)
        json_collection.append(this_json)

        json_opts_uri = os.path.join(output_bucket, f"input_allofus{idx!s}.json")
        with cloud_fs.open(json_opts_uri, 'w') as f:
            json.dump(this_json, f)

        if break_here:
            break

    #write batches to output file
    batch_json_uri = os.path.join(output_bucket, 'batch_input_allofus.json')
    with cloud_fs.open(batch_json_uri, 'w') as f:
        json.dump(json_collection, f)

    ct_count_uri = os.path.join(output_bucket, 'ct_submissions.txt')
    with cloud_fs.open(ct_count_uri, 'w') as f:
        f.write(str(len(json_collection)))

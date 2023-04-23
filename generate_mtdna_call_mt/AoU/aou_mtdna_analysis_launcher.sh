#!/bin/bash

#### PARAMETERS
export numTest=30 # number of samples per iteration
export numIter=1 # number of iterations
export JOBLIMIT=5000 # how many jobs to allow simulataneously
export outputFold=220920_30_36
export PORTID=8094
export USE_MEM=80
export SQL_DB_NAME="local_cromwell_run.db" # name of local SQL database
export FORCEDOWNLOAD='False' # if enabled, will force download for CRAM via gsutil -u {} cp
export RERUN_FAIL='False' # if enabled, will try to rerun failures. If disabled, will filter out failures.
export FILE_DONE="mt_pipeline_single_2_5_stats.tsv" # expects this to be in the current bucket
export FILE_FAIL="mt_pipeline_single_2_5_failures.tsv" # samples that have failed along with log files


#### INSTALL DEPENDENCIES
pip install pyhocon
git clone git clone https://github.com/rahulg603/mtSwirl.git


#### DOWNLOAD DATA
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/scatterWrapper_MitoPipeline_v2_5.wdl -o scatterWrapper_MitoPipeline_v2_5.wdl
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/check_variant_bounds.R -o check_variant_bounds.R
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/compatibilify_fa_intervals_consensus.R -o compatibilify_fa_intervals_consensus.R
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/jsontools.py -o jsontools.py
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/merge_per_batch.py -o merge_per_batch.py
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/check_overlapping_homoplasmies.R -o check_overlapping_homoplasmies.R
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/scripts/fix_liftover.py -o fix_liftover.py
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/files/NUMTv3_all385.hg38.interval_list -o NUMTv3_all385.hg38.interval_list
curl https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/files/input_allofus.json -o input_allofus.json
curl https://github.com/broadinstitute/cromwell/releases/download/77/cromwell-77.jar -o cromwell-77.jar -L
curl https://github.com/broadinstitute/cromwell/releases/download/77/womtool-77.jar -o womtool-77.jar -L
curl -s "https://get.sdkman.io" -o install_sdkman.sh
bash install_sdkman.sh


#### CONFIGURE AND UPLOAD
sed -i 's|WORKSPACEDIR|'"$WORKSPACE_BUCKET|" input_allofus.json
sed -i 's|REQPAYS|'"$GOOGLE_PROJECT|" input_allofus.json
gsutil cp check_variant_bounds.R "$WORKSPACE_BUCKET"/code/
gsutil cp compatibilify_fa_intervals_consensus.R "$WORKSPACE_BUCKET"/code/
gsutil cp check_overlapping_homoplasmies.R "$WORKSPACE_BUCKET"/code/
gsutil cp fix_liftover.py "$WORKSPACE_BUCKET"/code/
gsutil cp merge_per_batch.py "$WORKSPACE_BUCKET"/code/
gsutil cp jsontools.py "$WORKSPACE_BUCKET"/code/
gsutil cp NUMTv3_all385.hg38.interval_list "$WORKSPACE_BUCKET"/intervals/
gsutil -u $GOOGLE_PROJECT cp $CDR_STORAGE_PATH/wgs/cram/manifest.csv .
gcloud auth list --format json | jq -r .[0].account


#### PREPARE FILESYSTEM
mkdir "${outputFold}"


#### MUNGE INPUTS AND PREPARE COMMANDS
python <<'CODE'
import os, re
import json
import pandas as pd
from copy import deepcopy
from pyhocon import ConfigFactory, ConfigTree, HOCONConverter
from google.cloud import storage


# get globals
mem = int(os.getenv('USE_MEM'))
joblimit = int(os.getenv('JOBLIMIT'))
port = int(os.getenv('PORTID'))
n_test = int(os.getenv("numTest"))
n_iter = int(os.getenv("numIter"))
sql_db = os.getenv("SQL_DB_NAME")
bucket = os.getenv("WORKSPACE_BUCKET")
project = os.getenv("GOOGLE_PROJECT")
tf_force_dl = os.getenv("FORCEDOWNLOAD") == 'True'
path_completed_samples = os.getenv("FILE_DONE")
path_failed_samples = os.getenv("FILE_FAIL")
tf_rerun_fail = os.getenv("RERUN_FAIL") == 'True'
path_indiv_save = os.getenv("outputFold") + '/'
 
# set up cromwell directories
cromwell_test_workdir = bucket + "/"  # Later, "cromwell-executions" will be appended to this for cromwell-workflow storage.
output_bucket = bucket + "/" + os.getenv("outputFold")  # This is where the output of the WDL will be.
 
print(f'Workspace bucket: {bucket}')
print(f'Workspace project: {project}')
print(f'Workspace cromwell working bucket: {cromwell_test_workdir}')
print(f'Workspace output bucket: {output_bucket}')

# grab list of completed and failed samples
bucket_name_str = re.sub('^gs://', '', bucket)
storage_client = storage.Client()
bucket_obj = storage_client.bucket(bucket_name_str)
if storage.Blob(bucket=bucket_obj, name=path_completed_samples).exists(storage_client):
    df_stats = pd.read_csv(f"gs://{bucket_name_str}/{path_completed_samples}", sep='\t')
    sample_list = list(df_stats.s)
else:
    sample_list = []
if storage.Blob(bucket=bucket_obj, name=path_failed_samples).exists(storage_client):
    df_failed = pd.read_csv(f"gs://{bucket_name_str}/{path_failed_samples}", sep='\t')
    sample_list_failed = list(df_failed.s)
else:
    sample_list_failed = []

# set up filenames
options_filename = "options.json"
wdl_filename = "scatterWrapper_MitoPipeline_v2_5.wdl"
json_filename = "input_allofus.json"

# create options file
# options_content = f'{{\n  "jes_gcs_root": "{output_bucket}",\n  "workflow_failure_mode": "NoNewCalls"\n}}\n'
options_content = f'{{\n  "jes_gcs_root": "{output_bucket}"\n}}\n'

fp = open(options_filename, 'w')
fp.write(options_content)
fp.close()
print(options_content)

# create cromwell configuration, lifting the job limit and exporting to local SQL database
with open('/home/jupyter/cromwell.conf', 'r') as f:
    input_conf = f.read()
include_str_repl = 'include required(classpath("application"))\n\n'
input_conf_rm = input_conf.replace(include_str_repl, '')
cromwell_config_file = ConfigFactory.parse_string(input_conf_rm)
cromwell_config_file['system'] = ConfigTree({'new-workflow-poll-rate': 1,
                                             'max-concurrent-workflows': joblimit,
                                             'max-workflow-launch-count': 400,
                                             'job-rate-control': ConfigTree({'jobs': 50,
                                                                             'per': '3 seconds'})})
cromwell_config_file['backend']['providers']['PAPIv2-beta']['config']['concurrent-job-limit'] = joblimit
cromwell_config_file['backend']['providers']['PAPIv2-beta']['config']['genomics']['enable-fuse'] = True
cromwell_config_file['database'] = ConfigTree({'profile': "slick.jdbc.HsqldbProfile$",
                                               'insert-batch-size': 6000,
                                               'db': ConfigTree({'driver':"org.hsqldb.jdbcDriver", 
                                                                 'url':f'jdbc:hsqldb:file:{sql_db};shutdown=false;hsqldb.default_table_type=cached;hsqldb.tx=mvcc;hsqldb.large_data=true;hsqldb.lob_compressed=true;hsqldb.script_format=3;hsqldb.result_max_memory_rows=20000',
                                                                 'connectionTimeout': 300000})})
with open('/home/jupyter/cromwell.new.conf', 'w') as f:
    f.write(include_str_repl + HOCONConverter.to_hocon(cromwell_config_file))

# load cromwell configuration for modification
with open(json_filename, 'r') as json_file:
    jf = json.load(json_file)

# generate workflows and cromwell commands
df = pd.read_csv('manifest.csv')
original_size = df.shape[0]
df = df[~df.person_id.isin(sample_list)]
if tf_rerun_fail:
    df = df[df.person_id.isin(sample_list_failed)]
    new_size = df.shape[0]
    if not all(x in [x in list(df.person_id) for x in sample_list_failed]):
        raise ValueError("ERROR: All failed samples should be found in the sample list.")
    print(f'{str(new_size)} previously failed samples added to list for rerunning.')
else:
    df = df[~df.person_id.isin(sample_list_failed)]
    new_size = df.shape[0]
    if original_size-new_size != len(sample_list)+len(sample_list_failed):
        raise ValueError("ERROR: the number of samples removed from manifest must be the same as the number of processed samples + number of prior failures.")
    print(f'{str(len(sample_list) + len(sample_list_failed))} samples already processed or failed and removed from manifest.')
    print(f'{str(new_size)} samples remain in the manifest (of {str(original_size)} originally present).')
df = df.reset_index(drop=True)

cromwell_run_cmd = f'source "/home/jupyter/.sdkman/bin/sdkman-init.sh" &&  sdk install java 11.0.14-tem && echo "Validating WDL..." && java -jar womtool-77.jar validate _WDL_FILE_ && java -Xmx{str(mem)}g -classpath ".:sqlite-jdbc.jar" -Dconfig.file=/home/jupyter/cromwell.new.conf -Dwebservice.port={str(port)} -jar cromwell-77.jar server'
cromwell_run_cmd_final = cromwell_run_cmd.replace("_WDL_FILE_", wdl_filename)

with open(f"cromwell_startup_script.sh", "w") as text_file:
    text_file.write("#!/bin/bash\n")
    text_file.write(cromwell_run_cmd_final + '\n')

with open(f"cromwell_submission_script_individual_jobs.sh", "w") as text_file:
    text_file.write("#!/bin/bash\n")

json_collection = []
max_rows = df.shape[0]

for idx in range(0, n_iter):

    if (idx+1)*n_test >= max_rows:
        this_max = max_rows
        break_here = True
    else:
        this_max = (idx+1)*n_test
        break_here = False
    
    df_sub = df.iloc[idx*n_test: this_max]
    s = list(df_sub.person_id)
    cram_paths = list(df_sub.cram_uri)
    crai_paths = list(df_sub.cram_index_uri)

    f = open(f"{path_indiv_save}sample_list{str(idx)}.txt", "w")
    f.writelines('\n'.join([str(x) for x in s]) + '\n')
    f.close()

    f = open(f"{path_indiv_save}cram_file_list{str(idx)}.txt", "w")
    f.writelines('\n'.join(cram_paths) + '\n')
    f.close()

    f = open(f"{path_indiv_save}crai_file_list{str(idx)}.txt", "w")
    f.writelines('\n'.join(crai_paths) + '\n')
    f.close()

    dct_update = {'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_list': f"{path_indiv_save}cram_file_list{str(idx)}.txt",
                  'MitochondriaPipelineWrapper.wgs_aligned_input_bam_or_cram_index_list': f"{path_indiv_save}crai_file_list{str(idx)}.txt",
                  'MitochondriaPipelineWrapper.sample_name_list': f"{path_indiv_save}sample_list{str(idx)}.txt",
                  'MitochondriaPipelineWrapper.force_manual_download': tf_force_dl}

    this_json = deepcopy(jf)
    this_json.update(dct_update)
    json_collection.append(this_json)

    this_json_filename = f"input_allofus{str(idx)}.json"
    with open(path_indiv_save + this_json_filename, 'w') as f:
        json.dump(this_json, f)

    this_cromwell_cmd = f'curl -X POST "http://localhost:{str(port)}/api/workflows/v1" -H "accept: application/json" -F workflowSource=@_WDL_FILE_ -F workflowInputs=@_INPUTS_ -F workflowOptions=@_OPTIONS_FILE_'
    this_cromwell_cmd = this_cromwell_cmd.replace("_WDL_FILE_", wdl_filename)
    this_cromwell_cmd = this_cromwell_cmd.replace("_INPUTS_", this_json_filename)
    this_cromwell_cmd = this_cromwell_cmd.replace("_OPTIONS_FILE_", options_filename)

    with open(f"{path_indiv_save}cromwell_submission_script_individual_jobs.sh", "a") as text_file:
        text_file.write(this_cromwell_cmd + '\n')

    if break_here:
        break

batch_json_filename = f"batch_input_allofus.json"
with open(batch_json_filename, 'w') as f:
    json.dump(json_collection, f)
    
with open('ct_submissions.txt', 'w') as f:
    f.write(str(len(json_collection)))

batch_cromwell_cmd = f'curl -X POST "http://localhost:{str(port)}/api/workflows/v1/batch" -H "accept: application/json" -F workflowSource=@_WDL_FILE_ -F workflowInputs=@_INPUTS_ -F workflowOptions=@_OPTIONS_FILE_'
batch_cromwell_cmd = batch_cromwell_cmd.replace("_WDL_FILE_", wdl_filename)
batch_cromwell_cmd = batch_cromwell_cmd.replace("_INPUTS_", batch_json_filename)
batch_cromwell_cmd = batch_cromwell_cmd.replace("_OPTIONS_FILE_", options_filename)

count_n_submit = "submission_count=$(grep -o 'Submitted' batch_submission_ids.txt | wc -l)\n"
test_ct = 'if [ "$submission_count" -ne "$(cat ct_submissions.txt)" ]; then echo "ERROR: submission count is incorrect."; exit 1; fi\n'
get_batch_ids = 'cat batch_submission_ids.txt | sed \'s/{"id"://g\' | sed \'s/","status":"Submitted"}//g\' | sed \'s/"//g\' | sed \'s/,/\\n/g\' | sed \'s/\\[//g\' | sed \'s/\\]//g\' > ordered_batch_ids.txt\n'

with open(f"cromwell_submission_script_batch.sh", "w") as text_file:
    text_file.write("#!/bin/bash\n")
    text_file.write(batch_cromwell_cmd + ' | tee batch_submission_ids.txt\n')
    text_file.write(count_n_submit)
    text_file.write(test_ct)
    text_file.write(get_batch_ids)
    text_file.write('echo "" >> ordered_batch_ids.txt\n')

CODE


#### LAUNCH SERVER AS A SUBPROCESS
chmod 777 cromwell_startup_script.sh
setsid ./cromwell_startup_script.sh > cromwell_server_stdout.log 2>cromwell_server_stderr.log &


#### CREATE MONITORING COMMAND
sleep 150
echo "Server started."
echo "Here is the tail of the current stdout.log. Examine this to make sure the server is running:"
tail -n10 cromwell_server_stdout.log
echo ""
echo "Run cromwell_submission_script_batch.sh to submit the desired jobs."
echo ""
echo "Run the following command to track the progress of the various runs:"
echo ""
export success_file_pref="${outputFold}_prog_$(date +'%T' | sed 's|:|.|g')"
echo "python gnomad-mitochondria/gnomad_mitochondria/pipeline/cromwell_run_monitor.py --run-folder ${outputFold} --sub-ids ordered_batch_ids.txt --sample-lists ${outputFold}/sample_list{}.txt --check-success --output ${success_file_pref}" | tee check_workflow_status.sh
echo ""
echo "We have outputted this command in check_workflow_status.sh."
echo ""


#### CREATE UPLOADING COMMAND
echo "Upon completion of the workflow, don't forget to upload data file paths to gs:// !!"
echo "Use compile_paths.sh to do this, which will both merge files from this run and append to the database."
echo '#!/bin/bash' > compile_paths.sh
export tsvPREF="${WORKSPACE_BUCKET}/tsv/${outputFold}"
export htPREF="${WORKSPACE_BUCKET}/ht/${outputFold}"
echo "gsutil cp ${success_file_pref}'*' ${tsvPREF}/" >> compile_paths.sh
echo "python mtSwirl/generate_mtdna_call_mt/AoU/aou_collate_tables.py --pipeline-output-path ${success_file_pref}.success.tsv --file-paths-table-flat-output ${tsvPREF}/tab_batch_file_paths.tsv --per-sample-stats-flat-output ${tsvPREF}/tab_per_sample_stats.tsv --file-paths-table-output ${htPREF}/tab_batch_file_paths.ht --per-sample-stats-output ${htPREF}/tab_per_sample_stats.ht" >> compile_paths.sh
echo "python mtSwirl/generate_mtdna_call_mt/AoU/aou_update_sample_database.py --new-paths tsv/${outputFold}/tab_batch_file_paths.tsv --new-stats tsv/${outputFold}/tab_per_sample_stats.tsv --new-failures tsv/${outputFold}/${success_file_pref}.failure.tsv" >> compile_paths.sh
echo "" >> compile_paths.sh

import dxpy # type: ignore
import os
import json
import subprocess
import datetime

from constants import *
from query_projects import *


# used to get the number of batches ran
def check_progress_nexus(project_id=PROJECTS[PROJECT_NAMES[1]]):
    """Checks the number of folders in a specified DNANexus directory.

    Counts the number of folders within the `M3_500K_FOLDER_NAME` directory in the given project.

    Args:
        project_id (str, optional): The ID of the DNANexus project. Defaults to PROJECTS[PROJECT_NAMES[1]].

    Returns:
        int: The number of folders in the specified directory.
    """

    objects = dxpy.find_data_objects(
        project=project_id,
        classname='*',
        folder=M3_500K_FOLDER_NAME,
        recurse=False
    )
    folders = [obj for obj in objects if obj.get('type') == 'folder']
    return len(folders)

def munge_submission_script(json_name, ind, parallel=False, batch_size=""):
    if parallel:
        p_str = "parallelized_"
    else:
        p_str = ""

    this_nm = f"MitochondriaPipelineSwirl_Multi-{batch_size}_{ind}_m2ignored"
    return f"dx run workflows/{p_str}MitochondriaPipelineSwirl_v2.6_Multi/MitochondriaPipeline -y \
            --input-json-file {JSON_DIRECTORY}/{json_name} \
            --destination {M3_500K_FOLDER_NAME}/v2.6_Multi_batch_{this_nm} \
            --name {this_nm}", this_nm
            
# can queue up jobs in order or if specified
def find_jobs(jobs=None, n_samples=100, ignore_set=set()):
    samples_queued = []
    json_for_sample = []
    jobs_to_samples = {}
    
    # if we dont provide specific jsons to run, load jobs from dir and run them
    if jobs is None:
        # load in the generated job data from JSONs
        job_jsons = read_files_in_directory(JSON_DIRECTORY)
        
        # load in names of completed samples from both projects. M2 should be fixed, M3 changes
        completed_samples = set(read_file_at_path(COMPLETED_SAMPLE_LIST_PATH_03))
        completed_samples = completed_samples | set(read_file_at_path(M2_SAMPLE_LIST))
    else:
        job_jsons = jobs
    
    job_jsons.sort(key=extract_index)
    for job in job_jsons:
        if "Multi" in job:
            ind = extract_index(job)
            if ind not in ignore_set:
                with open(f"{JSON_DIRECTORY}/{job}", 'r') as f:
                    job_params = json.load(f)
                samples_in_curr_job = [x.split("_")[0] for x in job_params["stage-common.sample_name"]]
                samples_t = []
                # go through sample names individually and remove them if they are in the list of completed samples
                for sample in samples_in_curr_job:
                    if sample not in completed_samples:
                        samples_queued.append(sample)
                        samples_t.append(sample)

                    json_for_sample.append(job)
                jobs_to_samples[str(job)] = samples_t
                
                if len(samples_t) == 0:    
                    del jobs_to_samples[str(job)]
                    json_for_sample.remove(job)

                if len(samples_queued) >= n_samples:
                    break
            
    return samples_queued, json_for_sample, jobs_to_samples

def find_jobs_for_comparison(n_samples=200):
    """Finds jobs and samples for comparison up to a specified number of samples.

    Reads job JSON files from a specified directory, extracts sample names, and groups jobs by samples.

    Args:
        n_samples (int, optional): The maximum number of samples to include. Defaults to 200 (pretty arbitrary).

    Returns:
        tuple: A tuple containing three lists:
        - samples_queued: A list of sample names.
        - json_for_sample: A list of JSON file names corresponding to samples.
        - jobs_to_samples: A dictionary mapping job names to lists of sample names.
    """
    m2_sample_set = set(read_file_at_path(M2_SAMPLE_LIST))
    # for comparison, we look through the jsons and find samples that are in m2
    job_jsons = read_files_in_directory(JSON_DIRECTORY)
    
    samples_queued = []
    json_for_sample = []
    jobs_to_samples = {}
    
    # define a simple function to extract sample id from batch string
    name_no_batch = lambda x: x.split("_")[0]
    
    # load in names of completed samples
    completed_samples = set(read_file_at_path(COMPLETED_SAMPLE_LIST_PATH_03))
    completed_samples = completed_samples | set(read_file_at_path(M2_SAMPLE_LIST))

    job_jsons.sort(key=extract_index)
    for job in job_jsons:
        if "Multi" in job:
            with open(f"{JSON_DIRECTORY}/{job}", 'r') as f:
                job_params = json.load(f)
            samples_in_curr_job = job_params["stage-common.sample_name"]
            
            # go through sample names individually and remove them if they are in the list of completed samples
            samples_in_curr_job = [item for item in samples_in_curr_job if item not in completed_samples]
        else:
            continue
        # print(samples_in_curr_job)
        
        # per json, maintain a list of samples that are in M2's list
        sample_is_common = []
        for sample in samples_in_curr_job:
            # if the samplename is found in m2's sample name set,
            # add the sample and the job
            if name_no_batch(sample) in m2_sample_set:
                sample_is_common.append(sample)
                samples_queued.append(sample)
                json_for_sample.append(job)
        
        if len(sample_is_common) >= 0:    
            jobs_to_samples[str(job)] = sample_is_common
        if len(samples_queued) > n_samples:
            break
    
    # make sure that the samples we queued are actually done in Mito-02
    for s in samples_queued:
        assert name_no_batch(s) in m2_sample_set
    return samples_queued, json_for_sample, jobs_to_samples

def checkpoint_completed_jobs(samples_queued, verbose=False):
    def append_list_to_checkpoint(list_to_append):
        with open(PROGRESS_M3_FILE, 'a') as file:
            for item in list_to_append:
                file.write("%s\n" % item)
    if verbose:
        sample_list_curr = read_file_at_path(PROGRESS_M3_FILE)
        append_list_to_checkpoint(samples_queued)
        sample_list_post = read_file_at_path(PROGRESS_M3_FILE)
        print(f"Added {len(sample_list_post) - len(sample_list_curr)} samples to the completed list!")
    else:
        append_list_to_checkpoint(samples_queued)
                

def queue_jobs(refresh_checkpoint_list=False, comparison=False, n_samples=200, sleep=30, parallel=False, batch_size=30, ignore_set=set(), log_dir=LOG_DIRECTORY):
    """Queues jobs, ensuring consistent ordering and logging.

    Reads job configurations, selects a maximum number of samples, ensures consistent
    ordering, and submits jobs. Logs command execution and sleeps for
    a specified time between submissions.

    Args:
        n_samples (int, optional): The maximum number of samples to compare. Defaults to 200.
        
    Returns:
        None
    """
    # ensures consistent ordering
    json_to_ind = index_batches()
    
    # used to ensure that we have an updated list of completed sample names. takes a while though
    if refresh_checkpoint_list:
        generate_completed_sample_names(overwrite=True)
    
    # if we want to queue jobs for comparison, we use a function that calculates common samples' jsons first
    if comparison:
        ignore_set_inds = {extract_index_analysis(analysis_name, batch_size) for analysis_name in ignore_set}
        samples_to_queue, json_for_sample, jobs_to_sample = find_jobs_for_comparison(n_samples=n_samples, ignore_set_inds=ignore_set_inds)
        print(f"Queueing up {len(jobs_to_sample)} jobs with {len(samples_to_queue)} common samples.")
    else:
        ignore_set_inds = {extract_index_analysis(analysis_name, batch_size) for analysis_name in ignore_set}
        samples_to_queue, json_for_sample, jobs_to_sample = find_jobs(n_samples=n_samples, ignore_set=ignore_set_inds)
        print(f"Queueing up {len(jobs_to_sample)} jobs for a total of {len(samples_to_queue)} samples.")
        
    # main loop
    batch_names = []
    batch_ids = []
    for i, (sample_json, sample_names) in enumerate(jobs_to_sample.items()):
        if parallel:
            # sample_json = sample_json.replace("dx_jobs_100_updatedjson20/", "dx_jobs_100_updatedjson20_parallelized/")
            # job = json.loads(open(sample_json, "r").read())
            # job["stage-common.common.gatk_samtools_docker"] = 'dx://mitochondria-03:/Docker/gatk_4_6_0_samtools_1_9_220527_gnu_parallel.tar.gz'
            # with open(sample_json_parallel, "w") as outfile: 
                # json.dump(job, outfile)
            pass
        cmd, this_batch_nm = munge_submission_script(sample_json, json_to_ind[sample_json], parallel=parallel, batch_size=batch_size)
        os.makedirs(LOG_DIRECTORY, exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file = os.path.join(LOG_DIRECTORY, f"command_{timestamp}.log")

        with open(log_file, 'w') as f:
            # note: process.wait() only waits for the command to be ran, not for the job to be finished.
            process = subprocess.Popen(cmd, stdout=f, stderr=f, shell=True)
            process.wait()
            
        with open(log_file, 'r') as f:
            this_batch_id = f.read().split("Analysis ID: ")[1].strip()

        os.system(f"sleep {sleep}")
        print(f"Running {cmd}")
        print(f"[{i+1} / {len(jobs_to_sample)}] : Queued {sample_json}. Contains {len(sample_names)} samples:")
        print(f"{sample_names}\n{'=' * LINE_LENGTH}")
        
        if i == len(jobs_to_sample) - 1:
            print(f"Finished queueing batches! Now we wait...")
            
        batch_names.append(this_batch_nm)
        batch_ids.append(this_batch_id)
    
    checkpoint_completed_jobs(samples_to_queue, verbose=True)
    print(batch_ids)
    return dict(zip(batch_names, batch_ids))
        
    
# queues up some jobs that are in common
# queue_jobs_for_comparison()
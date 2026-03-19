import threading
import time
import json
import os
import dxpy
import datetime

from queue_jobs_m3 import *

def get_curr_time():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def get_job_status(job_id):
    job_desc = dxpy.describe(job_id)
    return job_desc.get('state')

def check_job_success(job_id):
    return get_job_status(job_id) == 'done'

def check_job_completed(job_id):
    return get_job_status(job_id) in ['failed', 'done', 'terminated']

class DNANexusJobManager:
    def __init__(self, queue_size, batch_size, n_machines, poll_interval=None, log_path=None, completed_jobs_path=None, failure_path=None, 
                 check_history=None, refresh_completed_jobs=False, reload_sample_list=False, log_dir=None):
        assert(queue_size % batch_size == 0)
        assert(check_history > n_machines * 1.2)
        self.STATUS_SET = {"done", "failed", "in_progress"}
        self.BATCH_SEARCH_PATTERN = "batch_.*"
        self.NUM_BATCH_FILES_GENERATED = 5
        self.ALL_DIR_SEARCH_PATTERN = f"*/{SEARCH_PATTERN}"
        
        self.batch_size = batch_size
        
        if poll_interval is None:
            self.poll_interval_monitor = int(10 * 60)
            # self.poll_interval_monitor = int(2 * 60)
            # self.poll_interval_check_completed = int(2.2 * 60)
        else:
            self.poll_interval_monitor = poll_interval["monitor"]
            # self.poll_interval_check_completed = poll_interval["update"]
            
        if failure_path is None:
            self.failed_jobs_path = "failed_jobs.json"
        if completed_jobs_path is None:
            self.completed_jobs_path = "completed_jobs.json"
        if log_dir is None:
            self.log_dir = f"./logs_queue_manager_{datetime.date.today()}/"
            os.makedirs(self.log_dir, exist_ok=True)
        if log_path is None:
            self.log_path = f"{self.log_dir}mito_03_batch_40_runs.log"
            
        self.check_history = check_history
        self.queue_size = queue_size
        
        if refresh_completed_jobs:
            new_jobs = self.get_status_jobs("done")
        else:
            new_jobs = None
        
        # root jobs is a dict of [json job description name] : [analysis ID]
        if os.path.exists(self.completed_jobs_path):
            with open(self.completed_jobs_path, "r") as f:
                self.completed_jobs = json.load(f)
            self.completed_set = set(self.completed_jobs.keys())
            self.num_completed_jobs = len(self.completed_jobs)
            print(f"Loaded {self.num_completed_jobs} completed jobs.")
        else:
            open(self.completed_jobs_path, 'a').close()
            self.completed_jobs = {}
            self.completed_set = set()
            self.num_completed_jobs = 0
                     
        self.active_jobs = self.get_status_jobs("in_progress")
        print(f'Found {str(len(self.active_jobs))} active jobs.')
        
        self.n_machines = n_machines

        self.failed_jobs = {}
        self.lock = threading.Lock()

        # load failed jobs from file if exists
        if os.path.exists(self.failed_jobs_path):
            with open(self.failed_jobs_path, "r") as f:
                self.failed_jobs = json.load(f)
                self.failing_set = set(self.failed_jobs.keys())
                print(f"Loaded {len(self.failed_jobs)} jobs to ignore.")
        else:
             open(self.failed_jobs_path, 'a').close()
             self.failed_jobs = set()

        if new_jobs is not None:
            print(f'Searched {str(self.check_history)} recently completed jobs and found {str(len(new_jobs))}.')
            new_jobs_f = {k: v for k, v in new_jobs.items() if k not in self.completed_set}
            print(f'Of these, {str(len(new_jobs_f))} were not previously recorded.')
            self.active_jobs.update(new_jobs_f)
            print(f'Tracking {str(len(self.active_jobs))} jobs, a subset of which may be completed.')
            self.update_completed_jobs()
                     
        self.write_to_log(f'INIT MESSAGE: Starting runs with n_machines={self.n_machines}, batch_size={self.batch_size}. Current time is {get_curr_time()}.')
        if reload_sample_list:
            self.completed_samples = generate_completed_sample_names(M3_500K_FOLDER_NAME, search_pattern=self.ALL_DIR_SEARCH_PATTERN, overwrite=True)
            print(f"Finished generating successful sample names. Current progress: {len(self.completed_samples)} samples")

    def write_to_log(self, message):
        if not os.path.exists(self.log_path):
            with open(self.log_path, 'w') as f:
                f.write(message)
                f.write("\n")
        else:
            with open(self.log_path, 'a') as f:
                f.write(message)
                f.write("\n")

    def submit_jobs(self, n_active_jobs=0):
        print(f"active jobs: {n_active_jobs}\tn machines: {self.n_machines}")
        # relies on monitor jobs to submit jobs. uses n_active_jobs to determine how many to submit
        if n_active_jobs < self.n_machines:
            available_slots = self.n_machines - n_active_jobs
            samples_to_submit = available_slots * self.batch_size
            self.write_to_log(f"SUBMIT_LOG\t{get_curr_time()}\tn samples submitting: {samples_to_submit}")
            submitted_jobs = self.queue_job_n_samples(samples_to_submit)
            
            for job, job_id in submitted_jobs.items():
                self.active_jobs[job] = job_id

            # dequeue jobs
            self.queue_size -= samples_to_submit

    def monitor_jobs(self):
        # uses dnanexus api to get list of active jobs so we know what to submit
        while True:
            active_jobs_t = self.get_status_jobs("in_progress")
            for job, job_id in active_jobs_t.items():
                self.active_jobs[job] = job_id
            n_active_jobs = len(self.active_jobs)

            self.update_completed_jobs()
            if n_active_jobs < self.n_machines and self.queue_size > 0:
                self.submit_jobs(n_active_jobs)

            time.sleep(self.poll_interval_monitor)

            print(f"MONITOR_LOG\t{get_curr_time()}\tn active jobs: {n_active_jobs}\tcurr queue size: {self.queue_size}")
            if self.queue_size == 0:
                return True

    def track_failures(self, job_name, job_id):
        # add to failing set so we can keep track of jobs to ignore
        self.failed_jobs[job_name] = job_id
        self.failing_set = set(self.failed_jobs.keys())    
        self.write_to_log(f"FAIL_LOG\t{get_curr_time()}\t{job_name} {job_id} failed. Ignoring in future.")
        with open(self.failed_jobs_path, "w") as f:
            json.dump(self.failed_jobs, f)
                
    def track_successes(self, job_name, job_id):
        self.completed_jobs[job_name] = job_id
        self.completed_set = set(self.completed_jobs.keys())  
        self.write_to_log(f"UPDATE_LOG\t{get_curr_time()}\t{job_name} : {job_id} completed.")
        with open(self.completed_jobs_path, "w") as f:
            json.dump(self.completed_jobs, f)

    def update_completed_jobs(self):
        currently_active_jobs = self.active_jobs.copy()
        
        for job_name, job_id in list(currently_active_jobs.items()): 
            # if job is not in progress, either done successfully or failed for some reason
            job_status = get_job_status(job_id)
            if job_status in ['failed', 'done', 'terminated']:
                self.num_completed_jobs += 1
                is_success = (job_status == 'done') and self.confirm_batch_files(job_name)
                self.write_to_log(f"UPDATE_LOG\t{get_curr_time()}\t{job_name} : {job_id} finished with success {is_success}")
                if is_success:
                    self.track_successes(job_name, job_id)
                    samples_completed_in_job = generate_completed_sample_names(f"{M3_500K_FOLDER_NAME}/v2.6_Multi_batch_{job_name}", save=True, overwrite=True)
                    with open(COMPLETED_SAMPLE_LIST_PATH_03, "a") as f:
                        for sample in samples_completed_in_job:
                            f.write(f"{sample}\n")
                else:
                    # if job failed, track failure
                    self.track_failures(job_name, job_id)

                # remove from active batches
                print(f"{get_curr_time()} Removing {job_name} from queue.")
                del self.active_jobs[job_name]
                # upload to gs://ukb_500k/queue_status
                

        self.completed_samples = read_file_at_path(COMPLETED_SAMPLE_LIST_PATH_03)
        self.write_to_log(f"UPDATE_LOG\t{get_curr_time()}\tnum completed samples: {len(self.completed_samples)}\tqueue size: {self.queue_size}\tcompleted {self.num_completed_jobs} batches this session")
            

    def confirm_batch_files(self, job_name):
        files = dxpy.find_data_objects(
            project=PROJECTS["mitochondria-03"],
            recurse=True,
            folder=f"{M3_500K_FOLDER_NAME}/v2.6_Multi_batch_{job_name}",
            name=self.BATCH_SEARCH_PATTERN,
            name_mode='regexp',
            describe=True
        )
        return len(list(files)) == self.NUM_BATCH_FILES_GENERATED


    def run(self):
        stop_event = threading.Event()
        try:
            job_monitor_thread = threading.Thread(target=self.monitor_jobs, daemon=True)
            # update_thread = threading.Thread(target=self.update_completed_jobs, daemon=True)
            
            job_monitor_thread.start()
            # update_thread.start()
            
            job_monitor_thread.join()
            # update_thread.join()
            
        except KeyboardInterrupt:
            print("Received interrupt signal. Stopping threads...")
            stop_event.set()
            print("All threads stopped successfully.")
            self.write_to_log(f"TERMIN_LOG\tAll threads stopped successfully.")
            
    def get_status_jobs(self, status):
        assert(status in self.STATUS_SET)
        # might need to add a small buffer here
        print(f"Searching most recent {self.check_history} jobs to find '{status}' analyses.")
        arr = ['dx', 'find', 'analyses', 
               '--project', PROJECT_NAMES[1], 
               '--origin-jobs',
               '--name', 'MitochondriaPipeline*',
               '--num-results', f'{self.check_history}',
               '--state', f'{status}']
        dx_api_list_analyses = subprocess.run(arr, capture_output=True, text=True).stdout.split('\n')
        
        # returns a map of job name -> analysis id
        return dict(
            (line.split(" (MitochondriaPipeline) ")[0].strip(' *'), line.split(f"({status}) ")[1])
            for line in dx_api_list_analyses 
            if f'({status}) analysis-' in line and 'm2ignored' in line and f'MitochondriaPipelineSwirl_Multi-{self.batch_size}' in line
        )
        
    def queue_job_n_samples(self, n_samples_to_queue):
        # when looking for batches to queue, ignore currently queued jobs and jobs known to fail and successful jobs
        ignore_set_t = set(self.active_jobs.keys()) | self.failing_set | self.completed_set
        batch_names = queue_jobs(refresh_checkpoint_list=False, 
                comparison=False, 
                n_samples=n_samples_to_queue, 
                sleep=60, 
                parallel=True, 
                batch_size=self.batch_size,
                ignore_set=ignore_set_t,
                log_dir=self.log_dir)
        return batch_names

        
if __name__ == "__main__":
    print("Starting.")
    job_manager = DNANexusJobManager(queue_size=360, batch_size=40, n_machines=3, 
                                     check_history=120, 
                                     refresh_completed_jobs=False,
                                     reload_sample_list=False)
    job_manager.run()
    

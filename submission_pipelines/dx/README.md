# DNANexus Job Manager

This toolset provides a job management system for DNANexus, allowing for efficiently monitoring and submission of batch jobs.

Developed by @gnchau.

## Features

- **Job Submission**: Automatically submits jobs based on available resources.
- **Job Monitoring**: Continuously checks the status of jobs and tracks completed and failed jobs.
- **Logging**: Records actions and job statuses to log files for future reference.
- **Sample Management**: Generates and tracks completed samples from successful jobs.

## Class: `DNANexusJobManager`

### Constructor Parameters
- `queue_size` (int): Total number of jobs to manage.
- `batch_size` (int): Number of samples to process per job.
- `n_machines` (int): Number of machines available for processing.
- `poll_interval` (dict): Optional intervals for monitoring.
- `log_path` (str): Optional path for log file.
- `completed_jobs_path` (str): Optional path for completed jobs log.
- `failure_path` (str): Optional path for failed jobs log.
- `check_history` (int): Number of past jobs to check for status.
- `refresh_completed_jobs` (bool): Refresh the list of completed jobs at initialization.
- `reload_sample_list` (bool): Reload the list of samples at initialization.
- `log_dir` (str): Directory for storing logs.

### Methods

- `run()`: Starts the job manager, initiating monitoring.
- `monitor_jobs()`: Continuously checks active jobs and submits new jobs as needed.
- `submit_jobs(n_active_jobs)`: Submits jobs based on the number of currently active jobs.
- `update_completed_jobs()`: Updates the status of completed jobs and tracks successes or failures.
- `track_successes(job_name, job_id)`: Records successful jobs in the log.
- `track_failures(job_name, job_id)`: Records failed jobs in the log.
- `confirm_batch_files(job_name)`: Checks if all batch files for a job are present.
- `get_status_jobs(status)`: Retrieves the current status of jobs.
- `queue_job_n_samples(n_samples_to_queue)`: Queues a specific number of samples for processing.

to use:
```
screen -r
python3 manage_queue.py
```

from util_funcs import *
import re

DESCRIBE_PROPERTIES = ['id', 'project', 'class', 'sponsored', 'name', 'types', 'state', 
                        'hidden', 'links', 'folder', 'tags', 'created', 'modified', 'createdBy', 'media', 
                        'archivalState', 'size', 'cloudAccount']

PROJECT_NAMES = ["mitochrondrial-02", "mitochondria-03"]
PROJECTS = {PROJECT_NAMES[0] : "project-G92qjkjJKqV4240b3jzyxb92", 
            PROJECT_NAMES[1] : "project-GbqfBQ8Jg52pv4zyq1fZJB9G"}
JSON_STAT_M2 = "../mito_02_stats.json"
JSON_STAT_M3 = "../mito_03_stats.json"

COMPLETED_SAMPLE_LIST_PATH_03 = "./completed_sample_list.txt"
SEARCH_DIRS_MITO_02 = ["/220618_MitochondriaPipelineSwirl/next8500_upto110k",
                        "/220618_MitochondriaPipelineSwirl/next8500_upto110k_p2",
                        "/220618_MitochondriaPipelineSwirl/next8500_upto110k_p3",
                        "/220618_MitochondriaPipelineSwirl/next28k_upto137k",
                        "/220618_MitochondriaPipelineSwirl/next9k_upto146k",
                        "/220618_MitochondriaPipelineSwirl/next27k_upto173k",
                        "/220618_MitochondriaPipelineSwirl/final27k_upto200k"]

STATISTICS_NAMES = {"coverage" : "batch_merged_mt_coverage.tsv.bgz",
                    "analysis_stats" : "batch_analysis_statistics.tsv",
                    "index_metrics" : "batch_idxstats_metrics.tsv.gz",
                    "yield_metrics" : "batch_yield_metrics.tsv.gz",
                "vcf" : "batch_merged_mt_calls.vcf.bgz"}

STATISTIC_FILE_SPECS = {STATISTICS_NAMES["coverage"] : ["per_base_coverage"],
                        STATISTICS_NAMES["analysis_stats"] : ["n_liftover_r2_pass", "mean_coverage"],
                        STATISTICS_NAMES["index_metrics"] : ["mapped_reads", "len", "chr"],
                        STATISTICS_NAMES["yield_metrics"] : ["singletons", "mate_diff_chr", "READ_LENGTH", "duplicates"],
                        STATISTICS_NAMES["vcf"] : []}

STATISTIC_FILE_DIR = "../statistics_both_projects"

M3_500K_FOLDER_NAME = "/finalized_parallel_MitochondriaPipelineSwirl/500k"
#                      /first_try_MitochondriaPipelineSwirl/500k
# M3_500K_FOLDER_NAME = "/comparison_MitochondriaPipelineSwirl/500k"


M2_COMPLETED_FOLDER = "./mito_02_completed/"
M2_FILE_LIST = f"{M2_COMPLETED_FOLDER}all_finished_sample_names.txt"
M2_SAMPLE_LIST = "./samples_completed_02.txt"
STATISTIC_FILE_DOWNLOADED_LIST = "./statistics_downloaded_file_ids.txt"
# JSON_DIRECTORY = "../dx_jobs_100_updatedjson20"
JSON_DIRECTORY = "../dx_jobs_40_batch_test_json"
PROGRESS_M3_FILE = "../m3_progress_mutable.txt"
LOG_DIRECTORY = "../logs/"
LINE_LENGTH = 20

COMPARISON_GRAPHS_OUT = "../comparison_plots"

def assert_correct_project():
    project_name, project_id = "mitochondria-03", PROJECTS["mitochondria-03"]
    assert project_name == PROJECT_NAMES[1]
    return project_name, project_id


# get list of json jobs in order to queue
def extract_index(filename) -> int:
    match = re.search(r'MitochondriaPipelineSwirl_Multi_(\d+)\.json', filename)
    return int(match.group(1)) if match else 0

def extract_index_analysis(filename, batch_size) -> int:
    match = re.search(rf'{batch_size}_(\d+)', filename)
    return int(match.group(1)) if match else 0

def index_batches():
    jsons = read_files_in_directory(JSON_DIRECTORY)
    jsons.sort(key=extract_index)
    return dict(zip(jsons, [i for i in range(len(jsons))])) 

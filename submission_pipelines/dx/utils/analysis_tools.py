import dxpy # type: ignore
import time
import subprocess
from collections import defaultdict
import json

from constants import *
from query_projects import *

import matplotlib.pyplot as plt # type: ignore


def read_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def save_plot(x_data, y_data, x_name, y_name, title, filename, out_dir, graph_type):
    if graph_type == "scatter":
        plt.scatter(x_data, y_data, s=3)
        plt.title(title)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.savefig(os.path.join(out_dir, filename), dpi=500)
        plt.close() 
    else:
        plt.bar(x_data, y_data)
        plt.title(title)
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.savefig(os.path.join(out_dir, filename), dpi=500)
        plt.close() 

    
def compare_stats(m2_stats_file=None, m2_stats=None, m3_stats_file=None, m3_stats=None):
    """Compares statistics from two datasets.

    Args:
        m2_stats_file (str, optional): Path to the M2 stats file. Defaults to None.
        m2_stats (pd.DataFrame, optional): M2 stats as a Pandas DataFrame. Defaults to None.
        m3_stats_file (str, optional): Path to the M3 stats file. Defaults to None.
        m3_stats (pd.DataFrame, optional): M3 stats as a Pandas DataFrame. Defaults to None.

    Raises:
        ValueError: If both stats and file are not provided for either M2 or M3.
    """

    # use stats if theyre given. otherwise, read from files
    if m2_stats is None and m2_stats_file is None:
        raise ValueError("Either m2_stats or m2_stats_file must be provided.")
    if m3_stats is None and m3_stats_file is None:
        raise ValueError("Either m3_stats or m3_stats_file must be provided.")

    if m2_stats is None:
        m2_stats = read_json(m2_stats_file)
    if m3_stats is None:
        m3_stats = read_json(m3_stats_file)

    common_samples = set([sname.split("_")[0] for sname in m2_stats.keys()]) & set([sname.split("_")[0] for sname in m3_stats.keys()])
    
    # the analysis here is the graphs but we can change this later
    os.makedirs(COMPARISON_GRAPHS_OUT, exist_ok=True)
    for sample in common_samples:
        for suff in list(STATISTICS_NAMES.values()):
            handle_stat_analysis(suff, sample, m2_stats, m3_stats)

    return COMPARISON_GRAPHS_OUT

    # * The VCF file may take some munging to convert into a format that can be easily compared.
    # * A few specific displays (x is old run, y is new run):
    # * mtDNA per-base coverage metrics in "batch_merged_mt_coverage".
    # * Variant counts in the column n_liftover_r2_pass in "batch_analysis_statistics.tsv"
    # * Nuclear coverage by adding up mapped_reads for each sample in "batch_idxstats_metrics"
    # * singletons  in "batch_yield_metrics"
    # * mate_diff_chr  in "batch_yield_metrics"

def handle_stat_analysis(stat_type, sample, m2_stats, m3_stats, graphs_out=COMPARISON_GRAPHS_OUT):
    # doesnt work well yet with multiple stat types per dataset
    assert stat_type in list(STATISTICS_NAMES.values())
    if stat_type == "batch_merged_mt_coverage.tsv.bgz":
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            m2_points = m2_stats[sample][stat_type][spec]
            m3_points = m3_stats[sample][stat_type][spec]
            
        filename = f"{sample}/{spec}_comparison.png"
        title = f"Comparison of {spec} from {stat_type}"
        x_name, y_name = f"{PROJECT_NAMES[0]} {spec}", f"{PROJECT_NAMES[1]} {spec}"
        save_plot(m2_points, m3_points, x_name, y_name, 
                        title, filename, graphs_out, "scatter")
            
    elif stat_type == "batch_analysis_statistics.tsv":
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            m2_points = m2_stats[sample][stat_type][spec]
            m3_points = m3_stats[sample][stat_type][spec]
        
        filename = f"{sample}/{spec}_comparison.png"
        title = f"Comparison of {spec} from {stat_type}"
        x_name, y_name = f"{PROJECT_NAMES[0]} {spec}", f"{PROJECT_NAMES[1]} {spec}"
        save_plot(m2_points, m3_points, x_name, y_name, 
                        title, filename, graphs_out, "scatter")
    else:
        # handles an aggregate case
        m2_points = m2_stats[sample][STATISTICS_NAMES["index_metrics"]]["mapped_reads"] + \
                    m2_stats[sample][STATISTICS_NAMES["yield_metrics"]]["singletons"] + \
                    m2_stats[sample][STATISTICS_NAMES["yield_metrics"]]["mate_diff_chr"]
        m3_points = m2_stats[sample][STATISTICS_NAMES["index_metrics"]]["mapped_reads"] + \
            m2_stats[sample][STATISTICS_NAMES["yield_metrics"]]["singletons"] + \
            m2_stats[sample][STATISTICS_NAMES["yield_metrics"]]["mate_diff_chr"]
            
        filename = f"{sample}/{spec}_comparison.png"
        title = f"Comparison of {spec} from {stat_type}"
        x_name, y_name = f"{PROJECT_NAMES[0]} {spec}", f"{PROJECT_NAMES[1]} {spec}"
        save_plot(m2_points, m3_points, x_name, y_name, 
                        title, filename, graphs_out)



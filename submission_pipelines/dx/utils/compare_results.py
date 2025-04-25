import dxpy # type: ignore
import os
import pandas as pd # type: ignore
import json
import pysam # type: ignore

from constants import *
from analysis_tools import *

# i want to go through the list of completed samples then go to the other project and download and compare the results
# per filename, we generate the graphs we want and save them in a directory
def open_file(file_path,
              local_directory="./comparisons/"):
    project_name, project_id = assert_correct_project()
    
    original_filename = os.path.basename(file_path)
    local_file_path = os.path.join(local_directory, original_filename)

    original_filename = os.path.basename(file_path)
    dxpy.download_dxfile(
        dxpy.dxlink({"project": project_id, "path": file_path}),
        local_file_path
    )
    # os.system("gunzip -c batch_analysis_coverage_statistics.tsv.bgz > batch_analysis_coverage_statistics.tsv")
    return pd.read_csv("batch_analysis_coverage_statistics.tsv", sep='\t', header=0)

    
def compare_sample_lists(path_m3_samples=COMPLETED_SAMPLE_LIST_PATH_03, path_m2_samples=M2_SAMPLE_LIST):    
    m3_sample_list = read_file_at_path(path_m3_samples)
    m2_sample_list = read_file_at_path(path_m2_samples)
    
    # one caveat: I stripped the batch suffix from m2 lsit so we do the same for M3
    m3_sample_list = list(map(lambda x: x.split("_")[0], m3_sample_list))
    
    # sort so we get a deterministic order
    return sorted(set(m3_sample_list) & set(m2_sample_list))
        
        
def get_common_sample_locations(common_samples):
    """Searches projects for the location of common samples. Doesn't break if it can't find a sample in a project.

    Retrieves sample locations for specified common samples from two DNANexus projects.
    
    Args:
        common_samples (list): A list of common sample names to find.

    Returns:
        tuple: A tuple of two tuples:
            - (samples_m3, sample_locations_m3): A tuple containing lists of sample names and locations for mito-03.
            - (samples_m2, sample_locations_m2): A tuple containing lists of sample names and locations for mito-02.
    """

    # get a list of the common samples in each projects and their respective locations
    samples_m3, sample_locations_m3 = find_parent_dirs_sample(common_samples, 
                                                              project_id=PROJECTS[PROJECT_NAMES[1]], 
                                                              search_dir=M3_500K_FOLDER_NAME,
                                                              ignore=set())
    
    # the m2 project is structured differently
    samples_m2, sample_locations_m2 = [], []
    for dna_dir in SEARCH_DIRS_MITO_02:
        samples_m2_t, sample_locations_m2_t = find_parent_dirs_sample(common_samples, 
                                                                        project_id=PROJECTS[PROJECT_NAMES[0]], 
                                                                        search_dir=dna_dir,
                                                                        ignore=set(samples_m2))
        samples_m2.extend(samples_m2_t)
        sample_locations_m2.extend(sample_locations_m2_t)
        # stop looking if all the samples are found. could be more optimal if in the find_parent_dirs_sample function
        if len(samples_m2) == len(common_samples):
            break
    return (samples_m3, sample_locations_m3), (samples_m2, sample_locations_m2)

def find_stats_file(stat, folder, project_name):
    results = dxpy.find_data_objects(project=PROJECTS[project_name],
                                name=f"*{stat}", 
                                folder=folder,
                                name_mode='glob', 
                                describe=True)
    
    # might break if it finds more than one for whatever reason, but it should be okay
    # for res in results:
    #     print(res['describe']['id'])
    for res in results:
        if res is not None:
            return res['describe']['id']
    # return list(results)[0]['describe']['id']

def collect_statistics(sample_names, sample_locations, project_name=None, downloaded_set=None, save=True):   
    print(STATISTIC_FILE_SPECS)  
    if project_name is None:
        raise ValueError("Provide a project name.")
    # apparently it costs to download files locally so lets keep track of them
    if downloaded_set is None:
        # if the file doesnt exist, we start with empty sets
        if not os.path.exists(STATISTIC_FILE_DOWNLOADED_LIST):
            with open(STATISTIC_FILE_DOWNLOADED_LIST, 'w') as _:
                pass
            downloaded_set = set()
            new_fname_dict = {}
        else:
            # if the file exists, we can read it and avoid redownloadind
            downloaded_map = list(map(lambda x: x.split("\t"), 
            read_file_at_path(STATISTIC_FILE_DOWNLOADED_LIST)))
            downloaded_set = set()
            new_fname_dict = {}
            for fid, new_fname in downloaded_map:
                new_fname_dict[fid] = new_fname
                downloaded_set.add(fid)
    else:
        downloaded_set = set()
        new_fname_dict = {}

    sample_stats = {}
    if project_name == PROJECT_NAMES[0]:
        munge_sample_name = lambda x: str(x) + "_23193_0_0"
    elif project_name == PROJECT_NAMES[1]:
        munge_sample_name = lambda x: str(x) + "_23372_0_0"
    else:
        raise ValueError(f"Please provide a project that's either {PROJECT_NAMES}.")
    
    for (sample, sample_location) in zip(sample_names, sample_locations):
        sample = munge_sample_name(sample)
        sample_stats[sample] = {}
        for suf in list(STATISTICS_NAMES.values()):
            statistic_file_id = find_stats_file(suf, sample_location, project_name)
            # it looks like some cases generated some but not all stat files, so next two lines handles it
            if statistic_file_id is None:
                continue
            print(f"Parsing {statistic_file_id}.")
            if statistic_file_id in downloaded_set:
                # if the statistic file has already been downloaded, get the local url
                new_fname = new_fname_dict[statistic_file_id]
                print(f"Found stats file {new_fname} for sample {sample}.")
            else:
                # if the statistic file has not been downloaded, generate a local url and download it
                new_dir = f"{STATISTIC_FILE_DIR}/{project_name}/{sample_location[1:].replace('/', '_').split('_')[-1]}"
                new_fname = f"{new_dir}/{suf}"
                make_dir(STATISTIC_FILE_DIR)
                make_dir(new_dir)
                dxpy.download_dxfile(statistic_file_id, filename=new_fname)
                print(f"Saved {statistic_file_id} to {new_fname}.\n(processing sample {sample})!")
                downloaded_set.add(statistic_file_id)
                new_fname_dict[statistic_file_id] = new_fname
            sample_stats[sample][suf] = handle_statistic_case(suf, new_fname, sample)
            
    with open(STATISTIC_FILE_DOWNLOADED_LIST, "w") as f:
        for fid, new_fname in new_fname_dict.items():
            f.write(f"{fid}\t{new_fname}\n")
            
    print(sample_stats[sample]['batch_yield_metrics.tsv.gz'])
    if save:
        path = JSON_STAT_M2 if project_name == "mitochrondrial-02" else JSON_STAT_M3
        if os.path.exists(path):
            with open(path, 'r') as f:
                existing_data = json.load(f)
                existing_data.update(sample_stats)
            with open(path, 'w') as f:
                json.dump(existing_data, f, indent=4)
        else:
            with open(path, 'w') as f:
                json.dump(sample_stats, f, indent=4)

    return path, sample_stats
    
    
def handle_statistic_case(stat_type, stat_local_url, sample):
    # kind of messy, will work for stats that we wanted for the comparison but can be adjusted pretty easily if we want more
    def read_zipped_to_dataframe(file_path):
        return pd.read_csv(file_path, header=0, sep='\t', compression='gzip')

    assert stat_type in list(STATISTICS_NAMES.values())
    res = {}
    if stat_type == STATISTICS_NAMES["coverage"]:
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            res[spec] = list(read_zipped_to_dataframe(stat_local_url)[sample])
        return res
    elif stat_type == STATISTICS_NAMES["analysis_stats"]:
        df = pd.read_csv(stat_local_url, sep='\t', header=0)
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            res[spec] = list(df[df['s']==sample][spec])
        return res
    elif stat_type == STATISTICS_NAMES["index_metrics"]:
        df = read_zipped_to_dataframe(stat_local_url)
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            res[spec] = list(df[df["s"]==sample][spec])
        return res
    elif stat_type == STATISTICS_NAMES["yield_metrics"]:
        df = read_zipped_to_dataframe(stat_local_url)
        for spec in STATISTIC_FILE_SPECS[stat_type]:
            res[spec] = list(df[df["s"]==sample][spec])
        return res
    elif stat_type == STATISTICS_NAMES["vcf"]:
        return analyze_heteroplasmy(sample, stat_local_url)
    else:
        raise ValueError("Not a valid statistic file suffix.")

def download_and_analyze_stats(refresh_m3=False, compare=True):
    # 5 mains steps here:
    # step 1: get a list of samples from completed runs in mito 3
    # step 2: get those results from runs in mito 2
    # step 3: find the common ones
    # step 4: collect the associated statistics
    # step 5: generate comparisons
    def load_batch_prefixes():
        return read_file_at_path(COMPLETED_SAMPLE_LIST_PATH_03)
    
    #  step 1
    if refresh_m3:
        completed_samples_m3 = generate_completed_sample_names(overwrite=True)
    else:
        completed_samples_m3 = load_batch_prefixes()
    print(f"Loaded a total of {len(completed_samples_m3)} samples from {PROJECT_NAMES[1]}!")
    
    # step 2
    # if sample list exists we shouldnt compute it again since it takes a while
    if os.path.isfile(M2_FILE_LIST):
        with open(M2_FILE_LIST, 'r') as f:
            return f.read()
    else:
        m2_fnames, m2_fids, m2_sample_names = get_completed_sample_names_m2(save=True)
    
    # step 3
    # find the common ones
    common_samples = compare_sample_lists()
    (samples_m3, sample_locations_m3), (samples_m2, sample_locations_m2) = \
        get_common_sample_locations(common_samples)
    
    # step 4
    m2_stats_file, m2_stats = collect_statistics(samples_m2, sample_locations_m2, project_name=PROJECT_NAMES[0])
    print(f"Saved {PROJECT_NAMES[0]} to {m2_stats_file}!")
    
    m3_stats_file, m3_stats = collect_statistics(samples_m3, sample_locations_m3, project_name=PROJECT_NAMES[1])
    print(f"Saved {PROJECT_NAMES[1]} to {m3_stats_file}!")
    
    m2_heteroplasmy = analyze_heteroplasmy(samples_m2, sample_locations_m2, project_name=PROJECT_NAMES[0])
    m3_heteroplasmy = analyze_heteroplasmy(samples_m3, sample_locations_m3, project_name=PROJECT_NAMES[1])

    # step 5: compare the results
    if compare:
        # i return the graph directory here but we can compute more advanced stats here later if we want

        # heteroplasmy is handled much differently
        save_plot(m2_heteroplasmy, m3_heteroplasmy, "scatter")
        out = compare_stats(m2_stats=m2_stats, m3_stats=m3_stats)
        if out:
            print(f"Success! Comparison graphs can be found at {out}.")
            
    # we can compute some more interesting stats here without comparison, for example more aggregate stats like nuclear cov.
    else:
        pass
    
def analyze_heteroplasmy(sample, sample_location):
    sample_heteroplasmy = []
    vcf_reader = pysam.VariantFile(sample_location)
    
    for record in vcf_reader:
        if record.chrom == "chrM" and record.pos == 302:
            if record.alleles[0] == "A" and record.alleles[1] == "AC":
                sample_heteroplasmy.append(record.samples[sample]["HL"])

    return sample_heteroplasmy

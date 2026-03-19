import dxpy # type: ignore
import os.path
from tqdm import tqdm # type: ignore
from tqdm.contrib import tzip # type: ignore
import pandas as pd # type: ignore

from constants import *

SEARCH_PATTERN = "[0-9]+_[0-9]+_[0-9]+_[0-9]+.*.final.vcf"
SEARCH_PATTERN = "batch_analysis_statistics.tsv"

# SEARCH_PATTERN = "*"

def generate_completed_sample_names(curr_directory,
                                    save=True,
                                    overwrite=False,
                                    search_pattern=SEARCH_PATTERN):
    """Generates a set of completed sample names based on file patterns.
    
    Searches for files matching a specific pattern (defined by `SEARCH_PATTERN`)
    within the specified `curr_directory` in the current DNANexus project.
    Determines completion based on the file name structure.
    Args:
        curr_directory (str, optional): The directory to search for files. Defaults to "/first_try_MitochondriaPipelineSwirl/500k/".
        save (bool, optional): If True, saves the results to a file. Defaults to True.
        overwrite (bool, optional): If True, overwrites the existing file. Defaults to False.
    Returns:
        set: A set of completed sample names.
    Raises:
        AssertionError: If `overwrite` is True but the target file doesn't exist.
    """
    project_name, project_id = assert_correct_project()
    print(f"Searching for completed batches in project [{project_name}] ({project_id}).")
    
    # recurse=true makes this very slow if there are many directories in the curr_dir, so make sure it's mostly flat files
    files = dxpy.find_data_objects(
        project=project_id,
        recurse=True,
        folder=curr_directory,
        name=search_pattern,
        name_mode='regexp',
        describe=True
    )
    
    tmp_batch_dir = "./tmp/verify_complete_samples/"
    os.makedirs(tmp_batch_dir, exist_ok=True)
    
    results = set()


    files = list(files)
    for f in tqdm(files):
        curr_fname = f['describe']['id']
        
        dxpy.download_dxfile(curr_fname, filename=f"{tmp_batch_dir}tmpfile.txt")
        samples = list(pd.read_csv(f"{tmp_batch_dir}tmpfile.txt", sep='\t', header=0)['s'])
        for sample in samples:
            sample_name_t = sample.split("_")
            if len(sample_name_t) == 4:
                results.add(sample_name_t[0])

    if save:
        if not os.path.isfile(COMPLETED_SAMPLE_LIST_PATH_03):
            os.system(f"touch {COMPLETED_SAMPLE_LIST_PATH_03}")
        if overwrite:
            print(f"Writing completed sample list to {COMPLETED_SAMPLE_LIST_PATH_03}.")
            assert os.path.isfile(COMPLETED_SAMPLE_LIST_PATH_03)
            lst = list(results)
            lst.sort()
            with open(COMPLETED_SAMPLE_LIST_PATH_03, 'a') as f:
                f.writelines("%s\n" % l for l in lst)

    return results
    

def generate_similar_substrings(substring=None):
    if not substring:
        raise("so substring provided")
    
    # all this does for now is search for a substring in both projects
    for project_name, project_id in PROJECTS.items():
        print(f"Searching in project: [{project_name}] ({project_id})")
        files = dxpy.find_data_objects(
            project=project_id,
            # folder='directory',
            name=f"*{search_term}*",
            name_mode='glob',
            describe=True
        )
        for file in files:
            print(f"[{project_name}]\tFound: {file['describe']['name']}")
    
# gets the names of completed samples from the previous runs and generates a list of these. 
# don't need to run it again if you don't want to      
def get_completed_sample_names_m2(project_id=PROJECTS["mitochrondrial-02"],
                                  save=False):
    """Finds completed sample names in a specified DNANexus project.
    
    Searches for files with specific patterns in the given project to identify completed samples.
    Args:
        project_id (str, optional): The ID of the DNANexus project. Defaults to PROJECTS["mitochrondrial-02"].
        save (bool, optional): If True, saves the results to a file. Defaults to False.
    Returns:
        list: A list of completed sample names.
    """
    def find_file(project_id, directory_path, file_name):
        files = dxpy.find_data_objects(
            project=project_id,
            recurse=True,
            folder=directory_path,
            name=file_name,
            name_mode='regexp',
            describe=True
        )
        file_ids = []
        file_names = []
        file_paths = []
        for f in files:
            file_names.append(f['describe']['name'])
            file_ids.append(f['describe']['id'])
            file_paths.append(f['describe']['folder'])
        return file_names, file_ids, file_paths
    
    def get_safe_fname(s):
        return s.replace("/", "_")

    sample_names = []
    
    # creates the checkpoint folder for reading in the old sample names
    if not os.path.exists(M2_COMPLETED_FOLDER):
        os.mkdir(M2_COMPLETED_FOLDER)

    
    for m2_folder in SEARCH_DIRS_MITO_02:
        print(f"Searching in {m2_folder} for batch analysis file.")
        file_names, file_ids, file_paths = find_file(project_id, m2_folder, "batch_analysis_statistics.tsv")
        
        # print(file['describe']['name'])
        for file_name, file_id, file_path in tzip(file_names, file_ids, file_paths):
            new_fname = f"{M2_COMPLETED_FOLDER}{get_safe_fname(m2_folder[1:])}_{get_safe_fname(file_path[1:])}_{file_name}"
            # print("new fname    " + new_fname)
            print(f"Saving {file_name} to {new_fname}!")
            dxpy.download_dxfile(file_id, filename=new_fname)
            print(f"Downloaded {new_fname}!")
            sample_names.extend(list(pd.read_csv(new_fname, sep='\t', header=0)['s']))
        
        print(f"Finished searching in {m2_folder}!")
        
    if save:
        with open(M2_FILE_LIST, "w") as f:
            for fn, fi, res in zip(file_names, file_ids, sample_names):
                f.write(f"{res}\t{fn}\t{fi}\n")

    return file_names, file_ids, sample_names
    
def find_parent_dirs_sample(substrings, project_id=PROJECTS[PROJECT_NAMES[1]], 
                            search_dir=None,
                            ignore=set()):
    """Finds parent directories for files containing specified substrings.
    
    Searches for files containing any of the provided substrings within a 
    specific directory in a DNANexus project. Returns a list of tuples where 
    each tuple contains the substring and its corresponding parent directory.
    Args:
        substrings (list): A list of substrings to search for in filenames.
        project_id (str, optional): The ID of the DNANexus project. Defaults to 
                                   PROJECTS[PROJECT_NAMES[1]]. (Assumes PROJECTS 
                                   and PROJECT_NAMES are defined elsewhere).
        search_dir (str, optional): The directory within the project to search. 
                                   Defaults to None to ensure proper behavior.
    Returns:
        list: A list containing substrings.
        list: A list containing the substrings' parent directories.
    """
    # working backwards a bit here because I forgot to record this 
    # when generating filenames but should be fast
    parent_dirs = set()

    # with open(substring_file, 'r') as f:
    #     substrings = f.read().splitlines()

    if search_dir is None:
        raise ValueError("Please provide a search directory and its associated project.")
    print(f"[In {project_id}]:")
    for i, substring in enumerate(substrings):
        if substring in ignore:
            print(f"[{i+1}/{len(substrings)}] {substring} already found in {search_dir}. Ignoring.")
            continue
        # Use dxpy.find_data_objects with appropriate filters to find files
        print(f"[{i+1}/{len(substrings)}] Searching for {substring} in {search_dir}.")
        results = dxpy.find_data_objects(project=project_id,
                                        name=f"*{substring}*", 
                                        folder=search_dir,
                                        name_mode='glob', 
                                        describe=True)

        for data_object in results:
            parent_dir = data_object['describe']['folder']
            parent_dirs.add((substring, parent_dir))

    res = list(parent_dirs)
    res.sort(key=lambda x: x[0])
    substrings = [pair[0] for pair in res]
    parent_dirs = [pair[1] for pair in res]


    return substrings, parent_dirs

def fix_mito_02_ids():
    sample_name_list = []
    fnames = read_files_in_directory(M2_COMPLETED_FOLDER)
    for file_name in fnames:
        file_path = os.path.join(M2_COMPLETED_FOLDER, file_name)
        if os.path.isfile(file_path) and file_path.endswith(".tsv"):
            print(file_path)
            sample_name_list.extend(list(pd.read_csv(file_path, sep='\t', header=0)['s']))

    with open(M2_SAMPLE_LIST, "w") as f:
        for sn in sample_name_list:
            name_no_batch = lambda x: x.split("_")[0]
            f.write(f"{name_no_batch(sn)}\n")
    return sample_name_list
    
               
# snames = generate_completed_sample_names(save=True, overwrite=True) 
# generate_similar_substrings("2307006")
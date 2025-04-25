import os

def read_file_at_path(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.read().splitlines()
        return lines
    except FileNotFoundError:
        print(f"File at path '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []
    
def read_files_in_directory(directory_path):
    try:
        # List only files in the specified directory
        files = [file for file in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, file))]
        return files
    except FileNotFoundError:
        print(f"Directory at path '{directory_path}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []
    
def make_dir(directory_path):
    try:
        os.makedirs(directory_path, exist_ok=True)
        print(f"Directory '{directory_path}' created successfully.")
    except Exception as e:
        print(f"An error occurred while creating the directory: {e}")
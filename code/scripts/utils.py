import os
import json
import pickle

def load_config(dirname=None):
    if dirname is None:
        dirname = os.getcwd()
    config_path = os.path.join(dirname, "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


def save_list(list, path):
    with open(path, "wb") as file:
        pickle.dump(list, file)


def load_list(path):
    with open(path, "rb") as file:
        new_list = pickle.load(file)
    return new_list


def first_file_is_newer(file1:str, file2:str):
    """Return True if file1 was created after file2 and False otherwise."""
    if os.path.exists(file2): 
        return os.path.getctime(file1) > os.path.getctime(file2)
    else:
        return False
"""directory_handlers.py: module for handling directory operations."""
import os


def file_directory_handler(file_name, folder=None, path="normal"):
    path_working = path_handler(path)
    try:
        new_file = path_working + folder + file_name
        directory = path_working + folder
    except TypeError:
        new_file = path_working + file_name
        directory = path_working
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        raise
    check_if_and_makedir(directory)
    return new_file


def check_if_and_makedir(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
        return


def check_if_and_makefile(filename, path=None, file_type=".csv", path_type="normal"):
    if path is None:
        file_directory = file_directory_handler(filename, path=path_type)
    else:
        file_directory = path + filename
    for i in range(1, 1000):
        file_directory_updated = file_directory + "-" + str(i) + file_type
        if os.path.exists(file_directory_updated):
            continue
        else:
            break
    return file_directory_updated


def get_path(path_type):
    if path_type == "parent":
        path_new = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    elif path_type == "normal":
        path_new = os.getcwd()
    elif path_type == "parent_parent":
        path_new = os.path.dirname(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    else:
        raise Exception("path name not recognized!")
    return path_new


def path_handler(path_type):
    """ Gets working path in your OS """
    if path_type == "parent":
        path_working = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    elif path_type == "normal":
        path_working = os.getcwd()
    else:
        raise Exception("path name not recognized!")
    return path_working
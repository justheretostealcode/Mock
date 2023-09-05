import os


def filter_files_by_type(files, file_type=".py"):
    splitted_files = [(file, os.path.splitext(file)) for file in files]
    rel_files = [(file, (splitted_file[0].split(".")[0], splitted_file[1])) for file, splitted_file in splitted_files if splitted_file[1] == file_type]
    return rel_files


def delete_cython_files_in_dir(dir):

    files_in_dir = os.listdir(dir)

    python_files = filter_files_by_type(files_in_dir, ".py")
    c_files = filter_files_by_type(files_in_dir, ".c")
    shared_lib_files = filter_files_by_type(files_in_dir, ".so")

    python_files_prefixes = [file[1][0] for file in python_files]

    for file in c_files + shared_lib_files:
        if file[1][0] in python_files_prefixes:
            joined_path = os.path.join(dir, file[0])
            os.remove(joined_path)
            print("Removed:", joined_path)

    pass



if __name__ == '__main__':

    # ToDo Call setup from here
    # Move files to intended locations

    # Delete .so and .c files being equally named as .py files

    directories = ["simulator/", "models/"]

    for dir in directories:
        delete_cython_files_in_dir(dir)


    pass
import sys
import os
import glob
import shutil
from pathlib import Path

# python 3.5-3.8 has a problem with __file__ being relative to the initial cwd, (which isnt saved anywhere AFAIK) so when the user does a chdir() all the relative __file__ paths are meaningless
# this is an attempt at a workaround, although it does fail if the user changed BOTH the PWD env variable and did os.chdir() before importing this module
inital_pwd = None
if os.name == 'nt': # "if windows"
    inital_pwd = os.getenv("CD") # this is the windows pwd var
else:
    inital_pwd = os.getenv("PWD")

intial_cwd = os.getcwd()
if inital_pwd != intial_cwd:
    initial_cwd = inital_pwd

def write(data, *, to=None, path=None, force=True):
    path = to or path
    # make sure the path exists
    if force: ensure_is_folder(parent_path(path))
    with open(path, 'w') as the_file:
        the_file.write(str(data))

def read(filepath):
    try:
        with open(filepath,'r') as f:
            output = f.read()
    except:
        output = None
    return output    
    

def remove(path):
    if not Path(path).is_symlink() and os.path.isdir(path):
        shutil.rmtree(path)
    else:
        try:
            os.remove(path)
        except:
            pass

def ensure_is_folder(path, *, force=True):
    """
    Note: very forcefull
    """
    parent_path = parent_folder(path)
    # root is always a folder
    if parent_path == path:
        return 
    
    # make sure the parent actually is a folder
    if not is_folder(parent_path):
        # recurse up
        ensure_is_folder(parent_path, force=force)
    
    # delete anything in the way
    if force and exists(path) and not is_folder(path):
        remove(path)
    
    # now make the folder
    os.makedirs(path, exist_ok=True)
    return path

def ensure_is_file(path, *, force=True):
    """
    Note: very forcefull
    """
    ensure_is_folder(parent_path(path), force=force)
    # delete a folder if its in the way
    if force and is_folder(path):
        remove(path)
    if not exists(path):
        write("", to=path)
    return path

def move_out_of_the_way(path, extension=".old"):
    if exists(path):
        new_path = path+extension
        move_out_of_the_way(new_path, extension)
        move(item=path, to=os.path.dirname(new_path), new_name=os.path.basename(new_path))
    
def clear_a_path_for(path, overwrite=False, extension=".old"):
    original_path = path
    paths = []
    while os.path.dirname(path) != path:
        paths.append(path)
        path = os.path.dirname(path)
    
    paths.reverse()
    for each_path in paths:
        if not exists(each_path):
            break
        if is_file(each_path):
            if overwrite:
                remove(each_path)
            else:
                move_out_of_the_way(each_path, extension)
    ensure_is_folder(os.path.dirname(original_path))
    if overwrite:
        remove(original_path)
    else:
        move_out_of_the_way(original_path, extension)
    
    return original_path

def copy(item, *, to, new_name="", force=True):
    if new_name == "":
        raise Exception('copy() needs a new_name= argument:\n    copy(item="location", to="directory", new_name="")\nif you want the name to be the same as before do new_name=None')
    elif new_name is None:
        new_name = os.path.basename(item)
    
    # get the full path
    to = os.path.join(to, new_name)
    # if theres a file in the target, remove it
    if force and exists(to):
        remove(to)
    # make sure the containing folder exists
    ensure_is_folder(os.path.dirname(to))
    if os.path.isdir(item):
        shutil.copytree(item, to)
    else:
        return shutil.copy(item, to)

def move(item, *, to=None, new_name="", force=True):
    if new_name == "":
        raise Exception('move() needs a new_name= argument:\n    move(item="location", to="directory", new_name="")\nif you want the name to be the same as before do new_name=None')
    elif new_name is None:
        new_name = os.path.basename(item)
    
    # get the full path
    to = os.path.join(to, new_name)
    # make sure the containing folder exists
    ensure_is_folder(os.path.dirname(to))
    shutil.move(item, to)

def exists(path):
    return os.path.exists(path)

def final_target_of(path):
    # resolve symlinks
    if os.path.islink(path):
        have_seen = set()
        while os.path.islink(path):
            path = os.readlink(path)
            if path in have_seen:
                return None # circular broken link
            have_seen.add(path)
    return path

def is_folder(path):
    # resolve symlinks
    final_target = final_target_of(path)
    if not final_target:
        return False
    return os.path.isdir(final_target)
# aliases    
is_dir = is_directory = is_folder
    
def is_file(path):
    return os.path.isfile(path)

def ls(path="."):
    if is_folder(path):
        return os.listdir(path)
    else:
        return [path]

def list_paths_in(path):
    if is_folder(path):
        return [ join(path, each) for each in os.listdir(path) ]
    else:
        return []

def list_basenames_in(path):
    if is_folder(path):
        return os.listdir(path)
    else:
        return []

def list_file_paths_in(path):
    if is_folder(path):
        return [ join(path, each) for each in os.listdir(path) if is_file(join(path, each)) ]
    else:
        return []

def list_folder_paths_in(path):
    if is_folder(path):
        return [ join(path, each) for each in os.listdir(path) if is_folder(join(path, each)) ]
    else:
        return []

# iterate is MUCH faster than listing for large folders
def iterate_paths_in(path, recursively=False):
    if is_folder(path):
        with os.scandir(path) as iterator:
            for entry in iterator:
                yield os.path.join(path, entry.name)
        if recursively:
            for each_sub_folder in iterate_folder_paths_in(path, recursively=True):
                yield from iterate_paths_in(each_sub_folder, recursively=True)

def iterate_basenames_in(path):
    if is_folder(path):
        with os.scandir(path) as iterator:
            for entry in iterator:
                yield entry.name

def iterate_file_paths_in(path, recursively=False):
    if is_folder(path):
        with os.scandir(path) as iterator:
            for entry in iterator:
                if entry.is_file():
                    yield os.path.join(path, entry.name)
        if recursively:
            for each_sub_folder in iterate_folder_paths_in(path, recursively=True):
                yield from iterate_file_paths_in(each_sub_folder, recursively=False)

def iterate_folder_paths_in(path, recursively=False, have_seen=None):
    if recursively and have_seen is None: have_seen = set()
    if is_folder(path):
        with os.scandir(path) as iterator:
            for entry in iterator:
                entry_is_folder = False
                if entry.is_dir():
                    entry_is_folder = True
                    each_subpath_final_target = os.path.join(path, entry.name)
                # check if symlink to dir
                elif entry.is_symlink():
                    each_subpath = os.path.join(path, entry.name)
                    each_subpath_final_target = final_target_of(each_subpath)
                    if is_folder(each_subpath_final_target):
                        entry_is_folder = True
                
                if entry_is_folder:
                    yield os.path.join(path, entry.name)
                    if recursively:
                        have_not_already_seen_this = each_subpath_final_target not in have_seen
                        have_seen.add(each_subpath_final_target) # need to add it before recursing
                        if have_not_already_seen_this:
                            yield from iterate_folder_paths_in(each_subpath_final_target, recursively, have_seen)

def glob(path):
    return glob.glob(path)

def touch(path):
    ensure_is_folder(dirname(path))
    if not exists(path):
        write("", to=path)

def touch_dir(path):
    ensure_is_folder(path)

def parent_folder(path):
    return os.path.dirname(path)

parent_path = parent_folder
dirname = parent_folder

def basename(path):
    return os.path.basename(path)

def name(path):
    filename, file_extension = os.path.splitext(os.path.basename(path))
    return filename

def extname(path):
    filename, file_extension = os.path.splitext(path)
    return file_extension

def without_ext(path):
    started_with_dot_slash = path.startswith("./")
    parent_folders = os.path.dirname(path)
    filename, file_extension = os.path.splitext(os.path.basename(path))
    output = os.path.join(parent_folders, filename)
    if not started_with_dot_slash and output.startswith("./"):
        output = output[2:]
    return output

def without_any_ext(path):
    started_with_dot_slash = path.startswith("./")
    parent_folders = os.path.dirname(path)
    filename = os.path.basename(path).split('.')[0]
    output = os.path.join(parent_folders, filename)
    if not started_with_dot_slash and output.startswith("./"):
        output = output[2:]
    return output

without_extension = without_ext # alias
without_any_extension = without_any_ext # alias

def path_pieces(path):
    """
    example:
        *folders, file_name, file_extension = path_pieces("/this/is/a/filepath.txt")
    """
    folders = []
    while 1:
        path, folder = os.path.split(path)

        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)

            break
    folders.reverse()
    *folders, file = folders
    filename, file_extension = os.path.splitext(file)
    return [ *folders, filename, file_extension ]

def join(*paths):
    return os.path.join(*paths)

def is_absolute_path(path):
    return os.path.isabs(path)

def is_relative_path(path):
    return not os.path.isabs(path)

def make_absolute_path(to, coming_from=None):
    # if coming from cwd, its easy
    if coming_from is None:
        return os.path.abspath(to)
    
    # source needs to be absolute
    coming_from_absolute = os.path.abspath(coming_from)
    # if other path is  absolute, make it relative to coming_from
    relative_path = to
    if os.path.isabs(to):
        relative_path = os.path.relpath(to, coming_from_absolute)
    return os.path.join(coming_from_absolute, relative_path)

def make_relative_path(*, to, coming_from=None):
    if coming_from is None:
        coming_from = get_cwd()
    return os.path.relpath(to, coming_from)

def get_cwd():
    return os.getcwd()

def walk_up_until(file_to_find, start_path=None):
    here = start_path or os.getcwd()
    if not os.path.isabs(here):
        here = os.path.join(os.getcwd(), file_to_find)
    
    while 1:
        check_path = os.path.join(here, file_to_find)
        if os.path.exists(check_path):
            return check_path
        
        # reached the top
        if here == os.path.dirname(here):
            return None
        else:
            # go up a folder
            here = os.path.dirname(here)

def local_path(*paths):
    import os
    import inspect
    
    cwd = os.getcwd()
    # https://stackoverflow.com/questions/28021472/get-relative-path-of-caller-in-python
    try:
        frame = inspect.stack()[1]
        module = inspect.getmodule(frame[0])
        directory = os.path.dirname(module.__file__)
    # if inside a repl (error =>) assume that the working directory is the path
    except (AttributeError, IndexError) as error:
        directory = cwd
    
    if is_absolute_path(directory):
        return join(directory, *paths)
    else:
        # See note at the top
        return join(intial_cwd, directory, *paths)
        

def line_count_of(file_path):
    # from stack overflow "how to get a line count of a large file cheaply"
    def _make_gen(reader):
        while 1:
            b = reader(2**16)
            if not b: break
            yield b
    with open(file_path, "rb") as file:
        count = sum(buf.count(b"\n") for buf in _make_gen(file.raw.read))
    
    return count

def get_home():
    return str(Path.home())

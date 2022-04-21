import os
import glob
import re

def find_files(fd, ext_list, pat_list=None, pat_v_list=None, abs_pat_list=None, abs_pat_v_list=None,
               N=None, flag_ext_add_dot=True):
    """
    finds a list of files given extensions/patterns

    :param fd: root folder to start search
    :param ext_list: list of extensions (e.g. .trc or .trc.gz for gromos trajectory files)
    :param pat_list: list of patterns to be included in the file name
    :param pat_v_list: list of patterns to be excluded in the file name
    :param abs_pat_list: list of patterns to be included in the abs file name
    :param abs_pat_v_list: list of patterns to be excluded in the abs file name
    :param flag_ext_add_dot: make sure that extensions in ext_list start with a "."
    :return: list of paths
    """
    for pat_variable in (pat_list, pat_v_list, abs_pat_list, abs_pat_v_list):
        if pat_variable is None:
            pat_variable = []
    if flag_ext_add_dot:
        ext_list_dot = []
        for ext in ext_list:
            if not ext.startswith('.'):
                ext = '.' + ext
            ext_list_dot.append(ext)
    else:
        ext_list_dot = ext_list
    for fd, FDs, Fs in os.walk(fd):
        if N is not None and N == 0:
            break
        for f_name in Fs:
            flag_ext = False
            for ext in ext_list:
                if f_name.endswith(ext):
                    flag_ext = True
                    break
            if flag_ext:
                flag = False
                for pat in pat_list:
                    if pat not in f_name:
                        flag = True
                        break
                if flag:continue
                for patv in pat_v_list:
                    if patv in f_name:
                        flag = True
                        break
                if flag:continue
                f_path = os.path.abspath(os.path.join(fd, f_name))
                for pat in abs_pat_list:
                    if pat not in f_path:
                        flag = True
                        break
                if flag:continue
                for patv in abs_pat_v_list:
                    if patv in f_path:
                        flag = True
                        break
                if flag:continue
                yield(f_path)
                if N is not None:
                    N -= 1
                if N == 0:
                    break

def find_additional_files(file_folder, file_pattern, flag_file=True):
    """
    finds additional files in a folder (based on a pattern using glob)
    :param file_folder: root folder (or file based on which the root folder is defined - in case flag_file==True)
    :param file_pattern: pattern used to find additional files (root folder + pattern)
    :param flag_file: if True (by default), file_folder argument assumed to be a file in the root folder
    """
    if flag_file:
        fd, fname = os.path.split(file_folder)
    else:
        fd = file_folder
    files = glob.glob(os.path.join(fd, file_pattern))
    return files

def get_ext(fname):
    try:
        return re.match(r'\S*\.(\S*)', fname).group(1)
    except:
        return ''
"""
Author: Drazen Petrov
"""

"""
to import SMArt package add the path to SMArt to your path

to convert back to GROMOS format:
from SMArt.md.data_st import gr_EnegyTrajectory
tre = gr_EnegyTrajectory()
tre.load_tre_npz(path_to_ene_trajectory_file) # *.tre.npz
tre.write_tre(path_to_new_tre_file) # *.tre
run gzip path_to_new_tre_file in terminal to get *.tre.gz file

trg = gr_EnegyTrajectory()
trg.load_trg_npz(path_to_freene_trajectory_file) # *.trg.npz
trg.write_trg(path_to_freene_trajectory_file) # *.trg
run gzip path_to_freene_trajectory_file in terminal to get *.trg.gz file
"""

import os
import glob
import numpy as np
from SMArt.md.data_st import gr_EnegyTrajectory
from SMArt.incl import test_time

def check_ext(f_name, extensions=('tre', 'trg')):
    for ext in extensions:
        if f_name.endswith(ext) or f_name.endswith(ext + '.gz'):
            return ext
    return False

def find_trj_files(fd, pat_list, pat_v_list, abs_pat_list, abs_pat_v_list, N=None):
    """
    :param fd: root folder to start search
    :param pat_list: list of patterns to be included in the file name
    :param pat_v_list: list of patterns to be excluded in the file name
    :param abs_pat_list: list of patterns to be included in the abs file name
    :param abs_pat_v_list: list of patterns to be excluded in the abs file name
    :return: list of paths
    """
    for fd, FDs, Fs in os.walk(fd):
        if N is not None and N == 0:
            break
        for f_name in Fs:
            if check_ext(f_name):
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

def reduce_tre_trg(lib_path, f_path, data_type, flag_rm, dt=None, skip=None, fname_prefix=None, fname_sufix=None):
    """
    reduce the trajectory file size - reads gromos tre or trg trajectory file and stores it as numpy binary
    :param lib_path: ene_ana library file
    :param f_path: trajectory file
    :param data_type: data type for values (e.g. float, np.float64, np.float32, np.float16)
    :param flag_rm: remove the original file (True/False)
    :param skip: write every `skip` steps in 
    :param dt: used to calculate `skip` param; write every dt steps in the numpy trajectory
    :param fname_prefix: prefix to add to the output trajectory file name
    :param fname_sufix: sufix to add to the output trajectory file name
    """
    ext = check_ext(f_path)
    if ext=='tre':
        treg = gr_EnegyTrajectory(lib_path, f_path, real_num_dtype=data_type)
        np_trj = treg.tre
    elif ext=='trg':
        treg = gr_EnegyTrajectory(lib_path, fene_trj_path=f_path, real_num_dtype=data_type)
        np_trj = treg.trg
    if skip is None and dt:
        np_trj['TIMESTEP']['TIME'][1, 1, 0] - np_trj['TIMESTEP']['TIME'][0, 1, 0]
        skip = int(dt / np_trj)
    if skip:
        setattr(treg, ext, np_trj[::skip])
    if f_path.endswith('.gz'):
        new_f_path = f_path[:-3]
    else:
        new_f_path = f_path
    if fname_prefix or fname_sufix:
        new_f_path = os.path.abspath(new_f_path)
        new_f_fd, new_f_name = os.path.split(new_f_path)
        if fname_prefix:
            new_f_name = fname_prefix + new_f_name
        if fname_sufix:
            assert new_f_path.endswith(ext)
            new_f_name = new_f_name[:-4] + fname_sufix + '.' + ext
        new_f_path = os.path.join(new_f_fd, new_f_name)
    print('\tnew trj file', new_f_path)
    if ext=='tre':
        treg.write_tre_npz(new_f_path)
    elif ext=='trg':
        treg.write_trg_npz(new_f_path)
    if flag_rm:
        os.remove(f_path)
    return treg


if __name__ == '__main__':
    #------------------------------------------------------
    #import argparse
    #parser = argparse.ArgumentParser()
    from SMArt.incl import ArgParser
    parser = ArgParser(fromfile_prefix_chars='@')
    parser.add_argument('-fd', type=str, help='folder to search')
    parser.add_argument('-trj_files', type=str, nargs='+', help='list of trajectory files')
    parser.add_argument('-lib', type=str, help='ene ana library file')
    parser.add_argument('-data_t', type=str, default='s', choices=['s', 'd', 'h'], help = 'data type to use: s-single (default), d-double, h-half (ISSUE WITH inf VALUES!)')
    parser.add_argument('-pat', type=str, nargs='+', help='pattern to be included in the file name')
    parser.add_argument('-pat_v', type=str, nargs='+', help='pattern to be excluded in the file name')
    parser.add_argument('-abs_pat', type=str, nargs='+', help='pattern to be included in the abs path of the file')
    parser.add_argument('-abs_pat_v', type=str, nargs='+', help='pattern to be excluded in the abs path of the file')

    parser.add_argument('-dt', type=float, help='skip every dt frames (in ps)')
    parser.add_argument('-skip', type=int, help='skip frames')
    parser.add_argument('-f_pref', type=str, help='prefix to add to the output trajectory file name')
    parser.add_argument('-f_suf', type=str, help='sufix to add to the output trajectory file name')

    parser.add_argument('-N', type=int, help = 'stop after N *trc(.gz) files (does not stop when -1)')
    parser.add_argument('-run', default = False, action = 'store_true', help = 'run, otherwise just print *trc(.gz) files')
    parser.add_argument('-remove', default = False, action = 'store_true', help = 'remove the original files, otherwise just convert and keep')
    parser.add_argument('-v', '--verbose', default=False, action = 'store_true', help = 'print more info')
    args = parser.parse_args()

    assert args.fd or args.trj_files

    for pat_name in ('pat', 'pat_v', 'abs_pat', 'abs_pat_v'):
        if getattr(args, pat_name) is None:
            setattr(args, pat_name, [])

    data_type = {'s':np.float32, 'd':np.float64, 'h':np.float16}[args.data_t]
    #------------------------------------------------------
    if args.trj_files:
        TRJs_gen = args.trj_files[:args.N]
    else:
        TRJs_gen = find_trj_files(args.fd, args.pat, args.pat_v, args.abs_pat, args.abs_pat_v, args.N)
    if args.verbose:
        print(args)
    for f_path in TRJs_gen:
        if args.verbose:
            print(f_path)
        if args.run:
            assert args.lib
            treg = reduce_tre_trg(args.lib, f_path, data_type, args.remove, args.dt, args.skip, args.f_pref, args.f_suf)
        else:
            print(f_path)
  

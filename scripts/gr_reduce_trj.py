"""
Author: Drazen Petrov
"""

"""
to import SMArt package to your path

to convert back to GROMOS format:
from SMArt.md.data_st import Trajectory
trc = Trajectory()
trc.load_trc_npz(path_to_trajectory_file.trc.npz)
trc.write_trc(path_to_new_trc_file.trc)
run gzip path_to_new_trc_file.trc in terminal to get *.trc.gz file
"""

import os
import glob
import numpy as np
from SMArt.md.data_st import Trajectory
from SMArt.md import parse_cnf
from SMArt.incl import test_time
#from SMArt.md import parse_trc

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
            if f_name.endswith('.trc') or f_name.endswith('.trc.gz'):
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

def find_cnf_files(trj_file, cnf_pattern):
    trj_fd, trj_fname = os.path.split(trj_file)
    cnf_files = glob.glob(os.path.join(trj_fd, cnf_pattern))
    return cnf_files

def find_last_solute(cnf_files, file_i=0, flag_check=True):
    cnf = parse_cnf(cnf_files[file_i])
    last_at = cnf.find_last_solute()
    if flag_check:
        for cnf_i, cnf_f in enumerate(cnf_files):
            if cnf_i != file_i:
                cnf = parse_cnf(cnf_f)
                temp_last_at = cnf.find_last_solute()
                assert last_at == temp_last_at
    return last_at


def reduce_trc(f_path, data_type, flag_rm, dt=None, skip=None, N_atoms=None, atom_sel=None,
               fname_prefix=None, fname_sufix=None):
    """
    reduce the trajectory file size - reads gromos trc trajectory file and stores it as numpy binary
    :param f_path: trajectory file
    :param data_type: data type for coordinates (e.g. float, np.float64, np.float32, np.float16)
    :param flag_rm: remove the original file (True/False)
    :param skip: write every `skip` steps in 
    :param dt: used to calculate `skip` param; write every dt steps in the numpy trajectory
    :param N_atoms: write only first N_atoms in the output trajectory
    :param atom_sel: atom selection to write in the output trajectory
    :param fname_prefix: prefix to add to the output trajectory file name
    :param fname_sufix: sufix to add to the output trajectory file name
    """

    trc = Trajectory(f_path, real_num_dtype = data_type)
    if N_atoms or atom_sel:
        trc_array = trc.reduce_N_atoms(N_atoms, atom_sel)
        trc.trj = trc_array
        trc._trj = trc_array
    if skip is None and dt:
        trc_dt = float(trc.trj['time'][1] - trc.trj['time'][0])
        skip = int(dt / trc_dt)
    if skip:
        trc.trj = trc.trj[::skip]
        trc._trj = trc.trj
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
            assert new_f_path.endswith('.trc')
            new_f_name = new_f_name[:-4] + fname_sufix + '.trc'
        new_f_path = os.path.join(new_f_fd, new_f_name)
    print('\tnew trj file', new_f_path)
    trc.write_trc_npz(new_f_path)
    if flag_rm:
        os.remove(f_path)


if __name__ == '__main__':
    #------------------------------------------------------
    #import argparse
    #parser = argparse.ArgumentParser()
    from SMArt.incl import ArgParser
    parser = ArgParser(fromfile_prefix_chars='@')
    parser.add_argument('-fd', type=str, help='folder to search')
    parser.add_argument('-trj_files', type=str, nargs='+', help='list of trajectory files')
    parser.add_argument('-data_t', type=str, default='s', choices=['s', 'd', 'h'], help = 'data type to use: s-single (default), d-double, h-half')
    parser.add_argument('-pat', type=str, nargs='+', help='pattern to be included in the file name')
    parser.add_argument('-pat_v', type=str, nargs='+', help='pattern to be excluded in the file name')
    parser.add_argument('-abs_pat', type=str, nargs='+', help='pattern to be included in the abs path of the file')
    parser.add_argument('-abs_pat_v', type=str, nargs='+', help='pattern to be excluded in the abs path of the file')

    parser.add_argument('-dt', type=float, help='skip every dt frames (in ps)')
    parser.add_argument('-skip', type=int, help='skip frames')
    parser.add_argument('-N_atoms', type=int, help='only a subselection of atoms (first N)')
    parser.add_argument('-atom_sel', type=int, nargs='+', help='used to define the subselection of atoms (starts with 0)')
    parser.add_argument('-cnf_pat', type=str, help='patter to find a cnf file (to get the number of atoms)')
    parser.add_argument('-flag_check_cnf', default=False, action='store_true', help='check if all cnf files give the same number of atoms')
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
            if args.cnf_pat:
                cnf_files = find_cnf_files(f_path, args.cnf_pat)
                last_solute_atom = find_last_solute(cnf_files, flag_check=args.flag_check_cnf)
                if args.verbose:
                    print('reduce to first N atoms:', last_solute_atom)
                N_atoms, atom_sel = last_solute_atom, None
            else:
                N_atoms, atom_sel = args.N_atoms, args.atom_sel
            reduce_trc(f_path, data_type, args.remove, args.dt, args.skip, N_atoms, atom_sel,
                       args.f_pref, args.f_suf)
        else:
            print(f_path)
  

"""
Author: Drazen Petrov
"""

"""
to import SMArt package to your path

to convert back to GROMOS format:
trc = Trajectory()
trc.load_trc_npz(path_to_trajectory_file.trc.npz)
trc.write_trc(path_to_new_trc_file.trc)
run gzip path_to_new_trc_file.trc in terminal to get *.trc.gz file
"""

import os
import numpy as np
from SMArt.md.data_st import Trajectory
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

def reduce_trc(f_path, data_type, flag_rm, dt=None, skip=None, N_atoms=None, atom_sel=None):
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
    print(new_f_path)
    trc.write_trc_npz(new_f_path)
    if flag_rm:
        os.remove(f_path)



if __name__ == '__main__':
    #------------------------------------------------------
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fd', type=str, required=True, help='folder to search')
    parser.add_argument('-data_t', type=str, default='s', choices=['s', 'd', 'h'], help = 'data type to use: s-single (default), d-double, h-half')
    parser.add_argument('-pat', type=str, nargs='+', help='pattern to be included in the file name')
    parser.add_argument('-pat_v', type=str, nargs='+', help='pattern to be excluded in the file name')
    parser.add_argument('-abs_pat', type=str, nargs='+', help='pattern to be included in the abs path of the file')
    parser.add_argument('-abs_pat_v', type=str, nargs='+', help='pattern to be excluded in the abs path of the file')

    parser.add_argument('-dt', type=float, help='skip every dt frames (in ps)')
    parser.add_argument('-skip', type=int, help='skip frames')
    parser.add_argument('-N_atoms', type=int, help='only a subselection of atoms (first N)')
    parser.add_argument('-atom_sel', type=int, nargs='+', help='used to define the subselection of atoms (starts with 0)')

    parser.add_argument('-N', type=int, help = 'stop after N *trc(.gz) files (does not stop when -1)')
    parser.add_argument('-run', default = False, action = 'store_true', help = 'run, otherwise just print *trc(.gz) files')
    parser.add_argument('-remove', default = False, action = 'store_true', help = 'remove the original files, otherwise just convert and keep')
    args = parser.parse_args()
    for pat_name in ('pat', 'pat_v', 'abs_pat', 'abs_pat_v'):
        if getattr(args, pat_name) is None:
            setattr(args, pat_name, [])

    data_type = {'s':np.float32, 'd':np.float64, 'h':np.float16}[args.data_t]
    #------------------------------------------------------
    TRJs_gen = find_trj_files(args.fd, args.pat, args.pat_v, args.abs_pat, args.abs_pat_v, args.N)
    for f_path in TRJs_gen:
        if args.run:
            red_kwargs = {}
            red_kwargs['dt'] = args.dt
            red_kwargs['skip'] = args.skip
            reduce_trc(f_path, data_type, args.remove, args.dt, args.skip, args.N_atoms, args.atom_sel)
        else:
            print(f_path)
  

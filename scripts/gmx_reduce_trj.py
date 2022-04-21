"""
Author: Drazen Petrov

uses gmx trjconv to reduce trajetory files is size
"""

import os
from helper_fnc import find_files, find_additional_files, get_ext
from SMArt.md.wrappers import run_gm_prog

def reduce_gmx_trc(f_path, flag_rm=False, dt=None, skip=None, atom_sel=None, tpr_file=None,
               fname_prefix='new_', fname_sufix=None, gmx_kwargs=None, **kwargs):
    """
    reduce the trajectory file size - reads gromos trc trajectory file and stores it as numpy binary
    :param f_path: trajectory file
    :param flag_rm: remove the original file (True/False)
    :param skip: write every `skip` steps in 
    :param dt: write every dt steps in the numpy trajectory
    :param atom_sel: atom selection to write in the output trajectory (0 by default - System)
    :param fname_prefix: prefix to add to the output trajectory file name
    :param fname_sufix: sufix to add to the output trajectory file name
    :param gmx_kwargs: extra kwargs to pass to gmx trjconf (flags)
    :kwargs:
        out_f_path: path to a new file that will be created (if not given, it is derived from the input f_path)
        out_ext: extension for the out file
        flag_add_ext: add ext to the end of the file (False by default)
        flag_replace_orig_file: move the new file to the original f_path (replace it)
    """
    if gmx_kwargs is None:
        gmx_kwargs = {}
    gmx_kwargs = gmx_kwargs.copy()
    gmx_kwargs['f'] = f_path
    if dt:
        gmx_kwargs['dt'] = dt
    if skip:
        gmx_kwargs['skip'] = skip
    out_trj_file = kwargs.get('out_f_path')
    if out_trj_file is None:
        abs_f_path = os.path.abspath(f_path)
        fd, fname = os.path.split(abs_f_path)
        ext = get_ext(fname)
        if ext:
            fname = fname[:-(len(ext)+1)]
        if fname_prefix is None:fname_prefix=''
        if fname_sufix is None:fname_sufix=''
        new_fname = fname_prefix + fname + fname_sufix
        out_ext = kwargs.get('out_ext')
        if out_ext:
            ext = out_ext
        flag_add_ext = kwargs.get('flag_add_ext')
        print(False==flag_add_ext)
        if kwargs.get('flag_add_ext'):
            print('HELLO!')
            print(kwargs.get('flag_add_ext'))
            print(flag_add_ext)
            if not ext.startswith('.'):
                ext = '.' + ext
            new_fname += ext
        out_trj_file = os.path.join(fd, new_fname)
    gmx_kwargs['o'] = out_trj_file
    if tpr_file:
        gmx_kwargs['s'] = tpr_file
    if atom_sel:
        pre_comm = 'echo "{:}" | gmx '.format(atom_sel)
    else:
        pre_comm = 'gmx '
    comm = run_gm_prog('trjconv', pre_comm=pre_comm, **gmx_kwargs)
    res = os.system(comm)
    if res==0:
        if flag_rm:
            os.remove(f_path)
            if kwargs.get('flag_replace_orig_file'):
                os.rename(out_trj_file, f_path)
            replace_orig_file_script = kwargs.get('replace_orig_file_script')
            if replace_orig_file_script:
                replace_orig_file_script.write('mv {:} {:}\n'.format(out_trj_file, abs_f_path))
    else:
        if kwargs.get('flag_stop_if_conv_fails', True):
            raise Exception('command failed:\n'+comm)


if __name__ == '__main__':
    #------------------------------------------------------
    #import argparse
    #parser = argparse.ArgumentParser()
    from SMArt.incl import ArgParser
    parser = ArgParser(fromfile_prefix_chars='@')
    parser.add_argument('-trj_files', type=str, nargs='+', help='list of trajectory files')
    parser.add_argument('-fd', type=str, help='folder to search')
    parser.add_argument('-in_file_ext', type=str, nargs='+', default=('trr',), help='list of file extensions - ("trr",) by default)')
    parser.add_argument('-pat', type=str, nargs='+', help='pattern to be included in the file name')
    parser.add_argument('-pat_v', type=str, nargs='+', help='pattern to be excluded in the file name')
    parser.add_argument('-abs_pat', type=str, nargs='+', help='pattern to be included in the abs path of the file')
    parser.add_argument('-abs_pat_v', type=str, nargs='+', help='pattern to be excluded in the abs path of the file')

    parser.add_argument('-dt', type=float, help='skip every dt frames (in ps)')
    parser.add_argument('-skip', type=int, help='skip frames')
    parser.add_argument('-atom_sel', type=str, help='used to define the subselection of atoms (e.g. 0 or System)')
    parser.add_argument('-tpr_pat', type=str, help='patter to find a tpr file')
    parser.add_argument('-f_pref', type=str, default='new_', help='prefix to add to the output trajectory file name')
    parser.add_argument('-f_suf', type=str, help='sufix to add to the output trajectory file name')
    parser.add_argument('-out_file_ext', type=str, help='out file extension (taken over from the in file if not defined and flag_out_ext set)')
    parser.add_argument('-flag_out_ext', default=False, action='store_true', help='use extension from the input file')
    parser.add_argument('-flag_replace_orig_file', default=False, action='store_true', help='replaces the original file with the out file')
    parser.add_argument('-replace_orig_file_script', help='writes commands that replace the original file with the out file')
    parser.add_argument('-flag_stop_if_conv_fails', default=True, action='store_false', help='stop execution if one of the commands fail')

    parser.add_argument('-N', type=int, help = 'stop after N *trc(.gz) files (does not stop when -1)')
    parser.add_argument('-run', default = False, action = 'store_true', help = 'run, otherwise just print *trc(.gz) files')
    parser.add_argument('-remove', default = False, action = 'store_true', help = 'remove the original files, otherwise just convert and keep')
    parser.add_argument('-v', '--verbose', default=False, action = 'store_true', help = 'print more info')
    args = parser.parse_args()

    assert args.fd or args.trj_files
    if args.flag_replace_orig_file:
        assert not args.replace_orig_file_script
        assert args.flag_out_ext and args.remove

    if args.replace_orig_file_script and args.run:
        replace_orig_file_script = open(args.replace_orig_file_script, 'w')
        replace_orig_file_script_fpath = os.path.abspath(args.replace_orig_file_script)
    else:
        replace_orig_file_script = None

    find_file_kwargs = dict()
    for pat_name in ('pat', 'pat_v', 'abs_pat', 'abs_pat_v'):
        if getattr(args, pat_name) is None:
            setattr(args, pat_name, [])
        find_file_kwargs[pat_name + '_list'] = getattr(args, pat_name)

    if args.trj_files:
        TRJs_gen = args.trj_files[:args.N]
    else:
        find_file_kwargs['N'] = args.N
        TRJs_gen = find_files(args.fd, args.in_file_ext, **find_file_kwargs)
    if args.verbose:
        print(args)
    
    reduce_gmx_kwargs = dict(flag_rm=args.remove, dt=args.dt, skip=args.skip, atom_sel=args.atom_sel,
        fname_prefix=args.f_pref, fname_sufix=args.f_suf, out_ext=args.out_file_ext,
        flag_add_ext=args.flag_out_ext, flag_replace_orig_file=args.flag_replace_orig_file,
        flag_stop_if_conv_fails=args.flag_stop_if_conv_fails, replace_orig_file_script=replace_orig_file_script)
    for f_path in TRJs_gen:
        if args.verbose:
            print(f_path)
        if args.run:
            tpr_file = None
            if args.tpr_pat or args.atom_sel:
                if args.tpr_pat is None:
                    args.tpr_pat = '*tpr'
                tpr_files = find_additional_files(os.path.abspath(f_path), args.tpr_pat)
                tpr_file = tpr_files[0]
            reduce_gmx_kwargs['tpr_file'] = tpr_file
            reduce_gmx_trc(f_path, **reduce_gmx_kwargs)
        else:
            print(f_path)
    
    if args.replace_orig_file_script and args.run:
        replace_orig_file_script.write('rm '+replace_orig_file_script_fpath)
        replace_orig_file_script.close()
        os.system('chmod u+x '+replace_orig_file_script_fpath)


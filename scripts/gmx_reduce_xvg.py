import os
from helper_fnc import find_files, find_additional_files, get_ext
from SMArt.md.gromacs.io.ana import read_xvg_data
import numpy as np
import pandas as pd


def reduce_gmx_xvg(f_path, data_type=np.float32, remove_ext=False, out_f_path=None, **kwargs):
    """
    reduce the xvg file size - reads gromos xvg file and stores it as numpy binary
    :param f_path: xvg file
    :param data_type: data type for values (e.g. float, np.float64, np.float32, np.float16)
    :param remove_ext: removes the original file extension in the output file name before adding .npz
    :param out_f_path: path to a new file that will be created (if not given, it is derived from the input f_path)
    """
    if out_f_path is None:
        abs_f_path = os.path.abspath(f_path)
        fd, fname = os.path.split(abs_f_path)
        if remove_ext:
            ext = get_ext(fname)
            if ext:
                fname = fname[:-(len(ext)+1)]
        out_f_path = os.path.join(fd, fname + '.npz')
    else:
        if not out_f_path.endswith('.npz'):
            out_f_path += '.npz'
    df = read_xvg_data(f_path, flag_pdDF=True)
    np.savez_compressed(out_f_path, cols=df.columns, data=df.values.astype(data_type))

def read_reduced_gmx_xvg(f_path, flag_pdDF=False):
    """
    reads reduced gmx xvg file
    :param f_path: path to the reduced gmx xvg file (.npz)
    :param flag_pdDF: return data as pandas DataFrame
    :return: cols - list of column names, data - numpy array with the data
    """
    loaded = np.load(f_path, allow_pickle=True)
    cols = loaded['cols'].tolist()
    data = loaded['data']
    if flag_pdDF:
        return pd.DataFrame(data, columns=cols)
    else:
        return (cols, data)


if __name__ == '__main__':
    #------------------------------------------------------
    import argparse
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-trj_files', type=str, nargs='+', help='list of trajectory files')
    parser.add_argument('-fd', type=str, help='folder to search')
    parser.add_argument('-in_file_ext', type=str, nargs='+', default=('xvg',), help='list of file extensions')
    parser.add_argument('-pat', type=str, nargs='+', help='pattern to be included in the file name')
    parser.add_argument('-pat_v', type=str, nargs='+', help='pattern to be excluded in the file name')
    parser.add_argument('-abs_pat', type=str, nargs='+', help='pattern to be included in the abs path of the file')
    parser.add_argument('-abs_pat_v', type=str, nargs='+', help='pattern to be excluded in the abs path of the file')

    parser.add_argument('-data_t', type=str, default='s', choices=['s', 'd', 'h'], help = 'data type to use: s-single (default), d-double, h-half (NOT RECOMMENDED FOR ENERGY!)')
    parser.add_argument('-remove_ext', default=False, action='store_true', help = 'removes the original file extension in the output file name before adding .npz')

    parser.add_argument('-N', type=int, help = 'stop after N files (does not stop when -1 or None)')
    parser.add_argument('-run', default = False, action = 'store_true', help = 'run, otherwise just print the files')
    parser.add_argument('-remove', default = False, action = 'store_true', help = 'remove the original files, otherwise just convert and keep')
    parser.add_argument('-v', '--verbose', default=False, action = 'store_true', help = 'print more info')

    args = parser.parse_args()

    if args.verbose:print(args)

    assert args.fd or args.trj_files
    data_type = {'s':np.float32, 'd':np.float64, 'h':np.float16}[args.data_t]

    if args.trj_files:
        TRJs_gen = args.trj_files[:args.N]
    else:
        find_file_kwargs = dict()
        for pat_name in ('pat', 'pat_v', 'abs_pat', 'abs_pat_v'):
            find_file_kwargs[pat_name + '_list'] = getattr(args, pat_name)
        find_file_kwargs['N'] = args.N
        TRJs_gen = find_files(args.fd, args.in_file_ext, **find_file_kwargs)

    for f_path in TRJs_gen:
        if args.verbose or not args.run:
            print(f_path)
        if args.run:
            reduce_gmx_xvg(f_path, data_type=data_type, remove_ext=args.remove_ext)
            if args.remove:
                os.remove(f_path)



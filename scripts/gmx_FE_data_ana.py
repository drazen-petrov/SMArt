import glob, os
from SMArt.md.ana import pert_FE
from SMArt.md.gromacs.io.ana import read_bar_data, read_xvg_data

def read_data(FD_list, data_file_pattern='*/FE_MD*.xvg', skip=None, stride=None, dl_max=None):
    bar_data = {}
    for sim_index, fd in enumerate(FD_list):
        data_files = glob.glob(fd + '/' + data_file_pattern)
        for f in data_files:
            dl_max_kw = {}
            if dl_max:
                dl_max_kw['dl_max'] = dl_max
            data, sim_l = read_bar_data(f, skip=skip, stride=stride, **dl_max_kw)
            pert_FE.combine_bar_dhdl(bar_data, data, sim_l, append_index=sim_index)
    return bar_data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.convert_arg_line_to_args = lambda arg_line:arg_line.split()
    parser.add_argument('-d', '--fd', dest='fd', help='target folder', type=str, required=True, nargs='+')
    parser.add_argument('-p', '--pattern', help='pattern to find data files in each of the folders', type=str, default='*/FE_MD*.xvg')
    parser.add_argument('-stride', type=int, help='stride step for reading the data files (read every stride steps)')
    parser.add_argument('-skip', type=int, help='skip first `skip` steps when reading the data files')
    parser.add_argument('-dl_max', type=float, help='max delta lambda for reading (how many neighbouring LPs to take into account)')
    parser.add_argument('-BS_steps', type=int, help='number of samples for bootstrapping', default=0)
    parser.add_argument('-n', '--n_cpu', type=int, help='number of CPUs for bootstrapping', default=8)
    parser.add_argument('-t', '--LPs_times', type=float, help='simulation time for each lambda point', nargs='+')
    parser.add_argument('-max_LPs_t', type=float, help='max total time for the next iteration', default=11.)
    parser.add_argument('-max_LP_t', type=float, help='max time for each LP for the next iteration', default=1.)
    parser.add_argument('-max_tot_LP_t', type=float, help='max total time for each LP', default=5.)
    parser.add_argument('-max_tot_LPs_t', type=float, help='max total time for all LPs for all iterations', default=100.)

    args = parser.parse_args()


    print('reading...')
    bar_data = read_data(args.fd, args.pattern, skip=args.skip, stride=args.stride, dl_max=args.dl_max)
    print('DONE reading\n\nanalysis...')

    pert_FE.dG_err_tols._dG_err_tols__BS_steps = args.BS_steps
    seg_width_slide_kw = {}
    if args.dl_max and args.dl_max>=1.:
        seg_width_slide_kw['seg_width_slide'] = [(1.,1.)]
    result = pert_FE.update_LPs_times(bar_data, ncpu=args.n_cpu, **seg_width_slide_kw)
    seg_score_flag, converged_segments, new_LPs, seg_data_dG_err, segs2calc_dG = result
    dG = pert_FE.get_full_dG_from_segs(seg_data_dG_err, segs2calc_dG[0])
    if args.BS_steps:
        err_methods = dict(full=['mbar_err'], BS={args.BS_steps:['mbar']})
    else:
        err_methods = dict(full=['mbar_err'])
    ddG = pert_FE.get_full_dG_err_from_segs(seg_data_dG_err, segs2calc_dG[0], err_methods)

    print('dG = ',dG, '+-', ddG)
    print('\nconverged segments')
    print(converged_segments)
    
    # LPs and times for next iteration
    if args.LPs_times:
        if len(args.LPs_times)==1:
            LP_t = args.LPs_times[0]
            LPs_times = dict((lp, LP_t) for lp in sorted(bar_data))
        else:
            assert len(args.LPs_times) == len(bar_data)
            LPs_times = dict(zip(sorted(bar_data), args.LPs_times))
        new_LPs_times_kwargs = dict(max_iter_LPs_t=args.max_LPs_t, max_iter_t_LP=args.max_LP_t, max_t_LP=args.max_tot_LP_t, max_total_LPs_t=args.max_tot_LPs_t)
        new_LPs_times = pert_FE.get_LPs_times(new_LPs, LPs_times, **new_LPs_times_kwargs)
        print('\nnew_LPs_times')
        print(new_LPs_times)



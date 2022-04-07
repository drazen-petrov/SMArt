import os
import SMArt
from SMArt.md import pipeline
from SMArt.md.data_st import MD_Parameters
from SMArt.md.wrappers import GMX_FE_sim_set_processor

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.convert_arg_line_to_args = lambda arg_line:arg_line.split()
    parser.add_argument('-d', '--fd', dest='fd', help='target folder', type=str, required=True)
    parser.add_argument('-t', '--top', dest='top', type=str, help='topology', required=True)
    parser.add_argument('-m', '--mdp', dest='mdp', type=str, help='input parameter file', required=True)
    parser.add_argument('-c', '--conf', dest='conf', type=str, help='configuration file', required=True)
    parser.add_argument('-r', '--rconf', dest='rconf', help='reference configuration file (for posres)', type=str)
    parser.add_argument('-w', '--maxwarn', help='max warn for gromacs', type=int, default=1)
    parser.add_argument('-eq_t', type=float, help='equlibration step time in ns', default=0.2)
    parser.add_argument('-sim_t', type=float, help='produciton simulation time in ns', default=1.)
    parser.add_argument('-max_sim_t', type=float, help='max simulation time in ns (e.g. due to slurm time limits)', default=None)
    parser.add_argument('-flag_cont', action='store_true', help='continuation of FE simualtions (using *cpt and appending into trajectories)', default=False)
    parser.add_argument('-N_LPs', type=int, help='number of lambda points', default=11)
    
    args = parser.parse_args()

    top = args.top.strip()
    top = os.path.abspath(top)

    gro = args.conf.strip()
    gro = os.path.abspath(gro)

    if args.rconf:
        ref = args.rconf.strip()
        ref = os.path.abspath(ref)
    else:
        ref = None

    mdp = args.mdp.strip()
    mdp = os.path.abspath(mdp)

    fd = args.fd.strip()
    fd = os.path.abspath(fd) + '/'

    flag_continue_cpt=args.flag_cont

    sim_set = pipeline.Initial_FE(args.eq_t, args.sim_t, N_lam=args.N_LPs, max_t_sim=args.max_sim_t)
    for sim in sim_set.sim_set:
        print(sim)

    sim_set.submit_cmd = 'sbatch --ntasks 1 --cpus-per-task 4 --mem 1000 --gres=mps:33 --partition NGN '
    for sim in sim_set.sim_set:
        if not sim['eq']:
            if flag_continue_cpt:
                flag_cp_all = (True, False)
            else:
                flag_cp_all = False
            sim['slurm_job_kwargs'] = dict(temp_fd='/scratch/${SLURM_JOBID}/', flag_cp_all=flag_cp_all)
            sim['flag_cp_before_new_sub'] = True

    GMX_sim_process = GMX_FE_sim_set_processor(gro, top, mdp)    

    print('\n\n\n')

    grompp_kwargs = {'maxwarn':args.maxwarn}
    if ref:
        grompp_kwargs['r'] = ref
    sim_set.generate_sim_files_jobs(fd, GMX_sim_process.fnc2process, grompp_kwargs=grompp_kwargs, flag_continue_cpt=flag_continue_cpt)


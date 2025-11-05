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
    parser.add_argument('-N_LPs', type=int, help='number of lambda points', default=11)
    parser.add_argument('-cmd_prefix', type=str, help='command prefix to be used for running/submitting the run files, e.g. sbatch', default="./")
    parser.add_argument('-gmx', type=str, help='path to gmx', default="gmx")
    parser.add_argument('-temp_fd', type=str, help='temp folder on a local disk to write during simulation and copy to main storate after simulation done (e.g. /scratch/\${SLURM_JOBID}/)')
    parser.add_argument('-nproc', type=int, help='if partial usage of a node is required, use this. also adjust the submit command (e.g. --cpus-per-task for slurm')
    
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

    sim_set = pipeline.Initial_FE(args.eq_t, args.sim_t, N_lam=args.N_LPs, max_t_sim=args.max_sim_t)
    for sim in sim_set.sim_set:
        print(sim)

    sim_set.submit_cmd = args.cmd_prefix
    for sim in sim_set.sim_set:
        if not sim['eq']:
            if args.temp_fd:
                sim['job_kwargs'] = dict(temp_fd=args.temp_fd, flag_cp_all=False)
                sim['flag_cp_before_new_sub'] = True

    GMX_FE_sim_set_processor.gmx_path = args.gmx
    GMX_sim_process = GMX_FE_sim_set_processor(gro, top, mdp)
    GMX_sim_process.additional_mdrun_kwargs = {}
    if args.nproc:
        GMX_sim_process.additional_mdrun_kwargs['nt'] = args.nproc

    print('\n\n\n')

    grompp_kwargs = {'maxwarn':args.maxwarn}
    if ref:
        grompp_kwargs['r'] = ref
    sim_set.generate_sim_files_jobs(fd, GMX_sim_process.fnc2process, grompp_kwargs=grompp_kwargs)


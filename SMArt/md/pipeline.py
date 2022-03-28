from SMArt.incl import os, np
from SMArt.md.ana.incl import get_lset, fix_float, Real

def get_slurm_job_settings(**slurm_kw):
    """
    :param slurm_kw: key value dictionary that will add the following lines: '#SBATCH -{:} {:}\n'
    :return: header string
    """
    s = '#!/bin/bash\n'
    for item in slurm_kw.items():
        s+= '#SBATCH -{:} {:}\n'.format(*item)
    return s + '\n'

def gen_slurm_job(f_path, commands, **slurm_kw):
    """
    :param f_path: run file
    :param commands: commands to be executed in the run file
    :param slurm_kw: additional slurm settings (see  get_slurm_job_settings function)
    :return: None
    """
    f = open(f_path, 'w')
    f.write(get_slurm_job_settings(**slurm_kw))
    for comm in commands:
        if not comm.endswith('\n'):
            comm+='\n'
        f.write(comm)
    f.close()
    os.system('chmod a+x '+ f_path)

def prep_pdb(pdb):
    pass

def eq(top, conf, md_in, eq_jobs):
    pass

def run_intial_FE(t_eq, t_sim, LPs = None, flag_run_eq=True,  **kwargs):
    """
    gets a simulation set for FE calculations
    :param t_eq: equilibation time [ns]
    :param t_sim: production run time [ns]
    :param LPs: list of lambda points
    :param flag_run_eq: run all eq simulations in one job
    :param kwargs:
        N_lam: number of equidistant lambda points (for LPs)
    :return:
        set of simulations (dictionary)
    """
    if LPs is None:
        N_lam = kwargs.get('N_lam', 11)
        LPs = get_lset(np.linspace(0, 1., N_lam))
    sim_set = []
    c = 0
    current_eq_sim_id = 0
    for i,l in enumerate(LPs):
        if i != len(LPs) - 1:
            if flag_run_eq:
                new_sim = dict(eq=True, l=l, t=t_eq, sub=(c+1, ), run=(c+2, ))
            else:
                new_sim = dict(eq=True, l=l, t=t_eq, sub=(c+1, c+2), run=())
            sim_set.append(new_sim)
        else:
            new_sim = dict(eq=True, l=l, t=t_eq, sub=(c+1, ), run=())
            sim_set.append(new_sim)
        c+=1
        new_sim = dict(eq=False, l=l, t=t_sim, sub=(), run=())
        sim_set.append(new_sim)
        c+=1
    return sim_set




    
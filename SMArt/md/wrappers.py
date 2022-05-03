from SMArt.incl import os
from SMArt.md.data_st import MD_Parameters

def _run_prog(prog, flag_character, pre_comm='', pipe=(None, None), **flags):
    comm = pre_comm + prog
    for flag, val in flags.items():
        comm += ' {:}{:} {:}'.format(flag_character, flag, val)
    if pipe[0]:
        comm += ' > ' + pipe[0]
    if pipe[1]:
        comm += ' 2> ' + pipe[1]
    return comm

def run_gr_prog(prog, pre_comm='', pipe=(None, None), **flags):
    """
        run GROMOS program
    :param prog: program to run (e.g. make_top)
    :param pipe: pipe the output to a file
    :param flags: flags to be passed to the program (e.g. dict(topo=path_to_topo_file))
    :return: command
    """
    return _run_prog(prog, '@', pre_comm=pre_comm, pipe=pipe, **flags)

def run_gm_prog(prog, pre_comm='gmx ', pipe=(None, None), **flags):
    """
        run GROMACS program
    :param prog: program to run (e.g. grompp)
    :param pipe: pipe the output to a file
    :param flags: flags to be passed to the program (e.g. dict(f=path_to_mdp_file))
    :return: command
    """
    return _run_prog(prog, '-', pre_comm=pre_comm, pipe=pipe, **flags)

def run_prog(prog, pipe=(None, None), format_type='gr', **flags):
    """
        run GROMOS/GROMACS program
    :param prog: program to run (e.g. make_top / grompp)
    :param pipe: pipe the output to a file
    :param flags: flags to be passed to the program (e.g. dict(f=path_to_mdp_file) or dict(topo=path_to_topo_file))
    :return: command
    """
    flag_character = {'gr':'@', 'gm':'-'}[format_type]
    pre_comm = dict(gr='', gm='gmx ')
    return _run_prog(prog, flag_character, pre_comm=pre_comm, pipe=pipe, **flags)

class GMX_FE_sim_set_processor:

    additional_mdrun_kwargs = dict(v='', nt=4)

    def __init__(self, conf, top, mdp, **kwargs):
        self.conf = os.path.abspath(conf)
        self.current_conf = os.path.abspath(conf)
        self.top = os.path.abspath(top)
        self.mdp = MD_Parameters(mdp)
        self.kwargs = kwargs
        self.dt = self.mdp.get_dt()
        self.LPs_pred = self.mdp.get_LPs_pred()

    def fnc2process(self, sim, sim_set, job_f, sim_fd, name, **kwargs):
        flag_continue_cpt = kwargs.get('flag_continue_cpt')
        tot_sim_time = sim['t']
        sys_name = name
        if flag_continue_cpt:
            tot_sim_time += sim.get('prev_FE_sim_t', 0)
            sys_name = sim.get('sys_name', name)
        mdp_kw = dict(nsteps = int(1000 * tot_sim_time / self.dt))
        mdp_kw['init-lambda-state'] = self.LPs_pred.index(sim['l'])
        mdp_out = os.path.join(sim_fd, name + '.mdp')
        self.mdp.change_mdp(mdp_kw, mdp_out = mdp_out)
        commands = []
        #grompp
        grompp_kwargs = kwargs.get('grompp_kwargs', self.kwargs.get('grompp_kwargs'))
        in_conf = sim.get('input_conf', self.current_conf)
        if not flag_continue_cpt:
            in_conf = sim.get('prev_FE_sim_final_conf', in_conf)
        comm = run_gm_prog('grompp', f=mdp_out, c=in_conf, p=self.top, o=sys_name, **grompp_kwargs)
        commands.append(comm)
        if sim['eq']:
            self.current_conf = os.path.join(sim_fd, sys_name + '.gro')
        #mdrun
        mdrun_kwargs = kwargs.get('mdrun_kwargs', self.kwargs.get('mdrun_kwargs', dict()))
        mdrun_kwargs = mdrun_kwargs.copy()
        mdrun_kwargs.update(self.additional_mdrun_kwargs)
        if flag_continue_cpt:
            mdrun_kwargs['cpi'] = sys_name + '.cpt'
        comm = run_gm_prog('mdrun', deffnm=sys_name, **mdrun_kwargs)
        commands.append(comm)
        # update sim setup dictionary for next FE sims
        if not sim['eq']:
            for new_sim_i in sim.get('sub'):
                new_sim = sim_set[new_sim_i]
                if 'input_conf' in sim:
                    new_sim['input_conf'] = sim['input_conf']
                if flag_continue_cpt:
                    new_sim['prev_FE_sim_t'] = tot_sim_time
                else:
                    new_sim['prev_FE_sim_final_conf'] = os.path.join(sim_fd, sys_name + '.gro')
        else:
            for new_sim_i in sim.get('sub'):
                new_sim = sim_set[new_sim_i]
                new_sim['input_conf'] = os.path.join(sim_fd, sys_name + '.gro')
        return commands

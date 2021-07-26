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


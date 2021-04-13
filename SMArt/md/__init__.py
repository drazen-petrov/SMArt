from . import ana
from .data_st import FF, BBdb, Topology, Configuration, Trajectory, MolSystem


def parse_ff(parse_from, parse_from_file = True, format_type = 'gr', **kwargs):
    """
    :param parse_from: this can also be gromos or gromacs stream
    :param parse_from_file:
    :param format_type: gr or gm for gromos and gromacs
    :param kwargs:
    :return: force field object
    """
    return FF(parse_from, parse_from_file = parse_from_file, format_type = format_type, **kwargs)

def parse_top(parse_from, parse_from_file = True, format_type = 'gr', **kwargs):
    """
    :param parse_from: this can also be gromos or gromacs stream
    :param parse_from_file:
    :param format_type: gr or gm for gromos and gromacs
    :param kwargs:
    :return: toplogy object
    """
    if kwargs.get('flag_molsys',False):
        top = MolSystem()
        top.parse_top(parse_from, parse_from_file=parse_from_file, format_type = format_type, **kwargs)
    else:
        top = Topology(parse_from, parse_from_file=parse_from_file, format_type = format_type, **kwargs)
    return top

def parse_mtb(parse_from, parse_from_file = True, **kwargs):
    """
    :param parse_from: this can also be gromos or gromacs stream
    :param parse_from_file:
    :param kwargs:
    :return: mtb
    """
    mtb = BBdb()
    ifp_file = kwargs.get('ifp_file')
    if ifp_file:
        mtb.parse_ifp(ifp_file)
    mtb.parse_mtb(parse_from, parse_from_file=parse_from_file, **kwargs)
    return mtb

def parse_trc(f_path, **kwargs):
    """
    :param parse_from: this can also be gromos stream
    :param parse_from_file:
    :param kwargs:
    :return: trajectory
    """
    return Trajectory(f_path, **kwargs)


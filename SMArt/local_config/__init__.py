from SMArt.incl import subprocess, os, do_warn, Defaults
from .__local_config_data import __bin_paths, _gm_FDs, _gr_FDs

if __bin_paths['gmx'] is None:
    __bin_paths['gmx'] = 'gmx'

try:
    p = subprocess.Popen([__bin_paths['gmx']], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    temp_liness = str(err).split('\\n')
    for temp_line in temp_liness:
        temp_lines = temp_line.split('\n')
        for l in temp_lines:
            if l.startswith('Executable'):
                temp = l.split()
                __bin_paths['_gmx_bin'] = temp[1]
            if l.startswith('Data prefix:'):
                temp = l.split()
                gmx_data = temp[2]
                _gm_FDs.append(os.path.join(gmx_data, 'share', 'gromacs', 'top'))
except:
    do_warn('gmx not found')
    __bin_paths['gmx'] = None

if __bin_paths['gpp'] is None:
    try:
        make_top_path = os.popen('which make_top').readlines()
        __bin_paths['gpp'] = os.path.split(make_top_path)[0]
        gr_fd =  os.path.split(make_top_path)[0]
        for i in range(3):
            gr_fd =  os.path.split(gr_fd)[0]
        print(gr_fd)
    except:
        do_warn('gpp not found')
        __bin_paths['gpp'] = None

if __bin_paths['gxx'] is None:
    try:
        md_path = os.popen('which md').readlines()
        __bin_paths['gxx'] = os.path.split(md_path)[0]
    except:
        do_warn('gxx not found')
        __bin_paths['gxx'] = None

if __bin_paths['pdb2pqr'] is None:
    try:
        pdb2pqr_path = os.popen('which pdb2pqr').readlines()
        __bin_paths['pdb2pqr'] = os.path.split(pdb2pqr_path)[0]
    except:
        do_warn('pdb2pqr not found')
        __bin_paths['pdb2pqr'] = None


class paths(Defaults):
    @classmethod
    def gmx(cls):
        return cls.__gmx

    @classmethod
    def gpp(cls):
        return cls.__gpp

    @classmethod
    def gxx(cls):
        return cls.__gxx

    @classmethod
    def pdb2pqr(cls):
        return cls.__pdb2pqr

    @classmethod
    def gr_FDs(cls):
        return cls.__gr_FDs

    @classmethod
    def gm_FDs(cls):
        return cls.__gm_FDs

paths._add_class_defaults(__bin_paths, flag_set=1)

from SMArt.md.gromos import FF_DB as gr_FF_DB

gr_FF_DB = os.path.abspath(gr_FF_DB.__path__[0])
_gr_FDs.append(gr_FF_DB)
_gr_FDs = tuple(_gr_FDs)

# gm_FDs are set in gm_io_defaults (make sure it's not defined/used from 2 places!!! ####FIX

__other_paths = dict(gr_FDs = _gr_FDs)
paths._add_class_defaults(__other_paths, flag_set=1)

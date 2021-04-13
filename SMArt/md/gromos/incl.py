from SMArt.incl import Defaults, _copy_file
from .__local_defaults import _gromos_block_names, _gromos_dvalues


class GromosBlockNames(Defaults):
    pass

GromosBlockNames._add_defaults(_gromos_block_names, flag_set=1)

class GromosDefaults(Defaults):
    _gr_block_names = GromosBlockNames()

GromosDefaults._add_defaults(_gromos_dvalues, flag_set=1)


def _add_flag_text(flag, text, at=True, flag_cp=False):
    if at and not flag.startswith('@'):
        flag = '@' + flag
    if not at and flag.startswith('@'):
        flag = flag[1:]
    if flag_cp:
        text = _copy_file(text, flag_cp)
    return ' ' + flag + ' ' + text + ' '


def _add_flags_kwargs(comm, flags, in_kwargs):
    for flag in flags:
        temp_v = in_kwargs.get(flag, None)
        if not (temp_v is None):
            comm += _add_flag_text(flag, temp_v)
    return comm


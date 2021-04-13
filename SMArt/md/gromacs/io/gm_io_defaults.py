from SMArt.incl import os, Defaults
from SMArt.local_config import paths, _gm_FDs

from SMArt.md.gromacs import FF_DB as gm_FF_DB
gm_FF_DB = os.path.abspath(gm_FF_DB.__path__[0])
_gm_FDs.append(gm_FF_DB)
_gm_FDs = tuple(_gm_FDs)
paths._add_class_defaults(dict(gm_FDs=_gm_FDs), flag_set=1)

##### put these defauls in a separate file... ####FIX
_include_directive = '#include'
_define_directive = '#define'
_gm_keywords = (_include_directive, _define_directive, '#undef', '#ifdef', '#ifndef', '#else', '#endif')
_gm_seg_types = ('include', 'gm_define', 'undef', 'if', 'if', None, None) # gm_directive missing
_gm_kw_segtype_map = dict(zip(_gm_keywords, _gm_seg_types))
_gm_comm_indicators = (';',)
_gm_eol = '\\\n'


class gm_FDs_Defaults(Defaults):
    @property
    def _gm_path_FDs(self):
        return self.__path_FDs

    @staticmethod
    def __find_file_gm_path(f_path, path_FDs):
        if os.path.isfile(f_path):
            return f_path
        for temp_fd in path_FDs:
            temp_path = os.path.join(temp_fd, f_path)
            if os.path.isfile(temp_path):
                return temp_path

    def _find_file_gm_path(self, f_path):
        return self.__find_file_gm_path(f_path, self.__path_FDs)

    @classmethod
    def _find_file_gm_path_cl(cls, f_path):
        return cls.__find_file_gm_path(f_path, cls.__path_FDs)

_def_dict = {'path_FDs': paths.gm_FDs()}
gm_FDs_Defaults._add_class_defaults(_def_dict, flag_set=True)


class gm_io_Defaults(Defaults):
    @property
    def _gm_eol(self):
        return self.__eol

    @property
    def _include_directive(self):
        return self._gm_kw[0]

    @property
    def _define_directive(self):
        return self._gm_kw[1]

    def _add_gm_eol(self, l):
        return l + self.__eol

    def _remove_gm_eol(self, l):
        if l[-2:] == self.__eol:
            return l[:-2] + "\t", True
        else:
            return l, False

    def _get_gm_comm_indicators(self, flag_cl = False):
        if flag_cl:
            return self.__class__.__comm_indicators
        return self.__comm_indicators

    @property
    def _gm_kw(self):
        return self.__gm_keywords

    @property
    def _gm_kw_map(self):
        return self.__gm_kw_segtype_map

    def _get_segtype(self, gm_kw, default_value = None):
        return self._gm_kw_map.get(gm_kw, default_value)

_def_dict = {'comm_indicators':(';',), 'eol':'\\\n', 'gm_keywords':_gm_keywords, 'gm_kw_segtype_map':_gm_kw_segtype_map}
gm_io_Defaults._add_class_defaults(_def_dict, flag_set=True)


class gm_FF_Defaults(list):
    def add_defs(self, gm_defaults):
        self.extend(gm_defaults)

    def write_defs(self):
        txt_format = ' {:>9}' * len(self) + '\n'
        return txt_format.format(*self)

from SMArt.incl import Defaults, get_id_string
from .gm_io_defaults import gm_io_Defaults
from .additional_classes import Define

class gmInteractionTypeWriter(Defaults):
    def __just_for_readability(self):
        self.atoms = None
        self.fnc_type = None
        self.p_type = None

    _gm_int_type_atom_format = '{:<6} '
    _gm_int_type_fnc_format = '{:<1}'
    _gm_int_type_p_format = {str: ' {:12s}', float: ' {:12.7e}', int: ' {:12d}'}
    _gm_int_type_gmdef_format = ' {:>6}'

    def _write_gm_atoms(self):
        temp_v = []
        for temp_at in self.atoms:
            temp_v.append(get_id_string(temp_at, flag_id='gm_id'))
        temp_txt = self._gm_int_type_atom_format * len(self.atoms)
        return temp_txt.format(*temp_v)

    def _write_gm_fnc(self):
        if self.fnc_type=='gr_fnc':
            return self.__class__.get_int_type_gr2gm()
        return self._gm_int_type_fnc_format.format(self.fnc_type)

    def __check_define_params(self, temp_def):
        temp_int_type = self.__class__()
        temp_int_type.p_type = self.p_type
        temp_int_type.add_params(self.p)
        return self.check_eq_params(temp_int_type.p)

    def _write_gm_params(self, defines = None, **kwargs):
        if defines and hasattr(self, 'gmdef'):
            if self.__check_define_params(defines[self.gmdef]):
                return self._gm_int_type_gmdef_format.format(self.gmdef)
        temp_txt = ''
        for i in self.p_type:
            temp_txt+= self._gm_int_type_p_format[i]
        return temp_txt.format(*self.p)

    def _write_gm_int_type(self, defines = None, **kwargs):
        if self.fnc_type=='gr_fnc':
            gm_int_type = self.gr2gm()
            return gm_int_type._write_gm_atoms() + gm_int_type._write_gm_fnc() + \
                   gm_int_type._write_gm_params(defines = defines, **kwargs)
        return self._write_gm_atoms() + self._write_gm_fnc() + self._write_gm_params(defines = defines, **kwargs)

    def gm_params2define(self, define_id = None, **kwargs):
        if define_id is None:
            define_id = self.get_gr2gm_def_pref() + str(self.id)
        split_line = [define_id]
        split_line.extend(self._write_gm_params().split())
        temp_def = Define(split_line)
        self.gmdef = define_id
        return temp_def


class gmCMAPWriter(gm_io_Defaults, gmInteractionTypeWriter):
    def __get_cmap_p_txt(self, n_cmap_perline = 10, **kwargs):
        n = n_cmap_perline
        params = list(self.p)
        temp_txt = ''
        while len(params) > n:
            ttxt = n*self._gm_int_type_p_format[float]
            temp_txt += self._add_gm_eol(ttxt.format(*params[:n]))
            params = params[n:]
        n = len(params)
        temp_txt += n * self._gm_int_type_p_format[float].format(*params[:n])
        return temp_txt

    def _write_gm_int_type(self, defines = None, **kwargs):
        """
        :param defines:
        :param kwargs:
            n_cmap_perline - how many cmap parameters per line, default is 10
        :return:
        """
        temp_txt = self._write_gm_atoms()
        temp_v = []
        temp_v.extend(self.grid_ind)
        temp_txt += self._add_gm_eol(' {:<6}' * 2)
        temp_txt = temp_txt.format(*temp_v)
        temp_txt += self.__get_cmap_p_txt(**kwargs)
        return temp_txt


from SMArt.incl import os, math, np, __VersionCompatibility, Defaults, do_warn
from SMArt.md.gromacs.io.incl_md_only import gmInteractionTypeWriter, gmCMAPWriter

kb = 0.00831451

def d_periodic(p1, p2, l):
    temp_diff = (p2 -p1) % l
    return min(temp_diff, l - temp_diff)

_gr_title_gm_system = 'sys_title'
_file_type_ext_map = dict(ifp='FF', mtb='BB', cnf='conf')

class DescriptionPart:
    """in gromos for title and in gromacs for system/first line of gro file"""
    container2write = _gr_title_gm_system

    def __init__(self, lines=None, file_type = None, form = None, **kwargs):
        if lines is None:
            lines = list()
        self.lines = lines
        self.file_type = file_type
        self.form = form

    def add_line(self, line):
        self.lines.append(line)

    def _add_source(self, f_path = None, additional_txt = '', **kwargs):
        if f_path:
            self.add_line('\tfrom ' + f_path + ' ({:})\n'.format(additional_txt))
            _, file_extension = os.path.splitext(f_path)
            file_extension = file_extension[1:]
            self.file_type = _file_type_ext_map.get(file_extension, file_extension)
        else:
            self.add_line('\tfrom string input ({:})\n'.format(additional_txt))

def _add_DescriptionPart_from_fpath(temp_description, f_path, additional_txt = '', **kwargs):
    temp_description.add_line('\tfrom ' + f_path + ' ({:})\n'.format(additional_txt))
    _, file_extension = os.path.splitext(f_path)
    file_extension = file_extension[1:]
    temp_description.file_type = _file_type_ext_map.get(file_extension, file_extension)

def _add_DescriptionPart(temp_stream, additional_txt = '', **kwargs):
    temp_description = DescriptionPart(**kwargs)
    if temp_stream.f_path is None:
        temp_description.add_line('\tfrom string input ({:})\n'.format(additional_txt))
    else:
        _add_DescriptionPart_from_fpath(temp_description, temp_stream.f_path, additional_txt, **kwargs)
    return temp_description

def check_ind_f_params(param1, param2, cutoff=0.00001, **kwargs):
    """does the checking for individual parameters"""
    params = [float(param1), float(param2)]
    for i in range(2):
        if params[i]!=0.:
            if math.log10(abs(params[i])) < -15:
                params[i]=0.
    if params[1]==0.:
        if params[0]!=0.:
            return False
    else:
        rel_diff = abs((params[1] - params[0]) / params[1])
        if rel_diff > cutoff:
            return False
    return True

def __check_if_eq(param1, param2, fnc = np.isclose, check_type = True, **kwargs):
    """checking params equality"""
    if not isinstance(param1, (list, tuple)):
        param1 = (param1,)
    if not isinstance(param2, (list, tuple)):
        param2 = (param2,)
    n = len(param1)
    if n != len(param2):
        return False
    else:
        for i in range(n):
            if check_type:
                if type(param1[i]) != type(param2[i]):
                    return False
                if not isinstance(param1[i], float) and param1[i]!=param2[i]:
                    return False
            if not fnc(param1[i], param2[i], **kwargs):
                return False
        return True

def check_if_eq_np(param1, param2, check_type = True, **kwargs):
    return __check_if_eq(param1, param2, fnc=np.isclose, check_type=check_type, **kwargs)

def check_if_eq(param1, param2, cutoff=0.00001, check_type = True, flag_np = False, **kwargs):
    if flag_np:
        return __check_if_eq(param1, param2, fnc=np.isclose, check_type=check_type, **kwargs)
    else:
        return __check_if_eq(param1, param2, fnc=check_ind_f_params, check_type=check_type, cutoff=cutoff)

_p_type_code = {"f": float, "i": int}

def _atoms_tuple(self):
    if len(self.atoms) == self.na:
        self.atoms = tuple(self.atoms)

def _add_atom(self, *atoms):
    for at in atoms:
        self.atoms.append(at)
        self._atoms_tuple()


class Interaction(__VersionCompatibility):
    """base class for each interaction type (bonds, angles...)
    each interaction type should be a new class inheriting from this class and defining na
    where na stands for number of atoms - 2 for bonds, 3 for angles, etc.
    attrib
        atoms - real atoms (e.g. from topology or building block
        na - number of atoms
        _sol - A, B - but could have more as well
        int_type - interaction type (e.g. BondType)
    """
    #__slots__ = ('atoms', '_sol', 'na', 'int_type', 'container2write')

    add_atom = _add_atom
    _atoms_tuple = _atoms_tuple

    def __init__(self, int_type, atoms=None, states=None):
        self.int_type = int_type
        self.na = self.int_type.na
        self.container2write = int_type.find_int_container2write()
        if atoms:
            self.atoms = atoms
            self._atoms_tuple()
        else:
            self.atoms = []
        if states:
            self.states = states  # interaction parameters [state1, state2,...]; state = (fnc type, type code, topB params)
        else:
            self.states = []

    def add_state(self, state = None, **kwargs):
        if state is None:
            state = self.int_type(**kwargs)
        self.states.append(state) # state = interaction parameters (fnc type, type code, params)

    def check_atoms_int_match(self, atoms, **kwargs):
        return self.int_type._cls_check_atoms_int_match(self.atoms, atoms, **kwargs)


from SMArt.md.gro2gro.g2g import IntType_g2g


class InteractionType(gmInteractionTypeWriter, IntType_g2g):
    """base class for each interaction type (bonds, angles...)
    each interaction type should be a new class inheriting from this class and defining na, p_type and flag
    p_type is variable type for each of the parameters, for dihedrals (float, float, int)
    flag is gromacs-related flag representing if a parameter is allowed in ptp change -  (for dihedrals (1,1,0))"""

    def __init__(self, int_code=None, params=None, fnc_type = 'gr_fnc', atoms = None, **kwargs):
        self.fnc_type = fnc_type
        self.add_p_type_flag(self.fnc[fnc_type])
        if len(self.p_type) != len(self.flag):
            temp_m = "class " + self.__class__.__name__ + ": len(p_type) != len(flag)"
            raise Exception(temp_m)
        self.id = int_code
        self.p = []
        self.atoms = []
        if params:
            self.set_params(params)
        if atoms:
            if hasattr(atoms,'__iter__'):
                self.add_atom(*atoms)
            else:
                self.add_atom(atoms)

    def __just_for_readability(self):
        self.fnc = dict()
        self.na = 0

    def __str__(self):
        return str(self.__class__.__name__) + ' ' +str(self.p)

    def __repr__(self):
        return str(self.__class__.__name__) + ' ' +str(self.p)

    add_atom = _add_atom
    _atoms_tuple = _atoms_tuple

    def check_eq_params(self, params2, cutoff=0.00001, check_type = True, flag_np = False, **kwargs):
        if not flag_np:
            return check_if_eq(self.p, params2, cutoff=cutoff, check_type=check_type)
        else:
            return check_if_eq_np(self.p, params2, check_type = check_type, **kwargs)

    def check_eq_ptp_flag_params(self, params2, cutoff=0.00001, check_type = True, flag_np = False, **kwargs):
        if len(self.p)!=len(params2):
            return False
        flag_params1 = []
        flag_params2 = []
        for i, temp_param in enumerate(self.p):
            if self.flag[i]==False:
                flag_params1.append(temp_param)
                flag_params2.append(params2[i])
        if not flag_np:
            return check_if_eq(flag_params1, flag_params2, cutoff=cutoff, check_type=check_type)
        else:
            return check_if_eq_np(flag_params1, flag_params2, check_type = check_type, **kwargs)

    @staticmethod
    def __calc_diff_params(param1, param2):
        return (param1 - param2)**2

    def calc_ptp_score(self, int_type2, **kwargs):
        if int_type2 is None:
            return [np.inf] * len(self.p)
        params2 = int_type2.p
        score = []
        if kwargs.get('flag_ptp'):# use only allowed flags
            for i, param1 in enumerate(self.p):
                if self.flag[i]:
                    score.append(self.__calc_diff_params(param1, params2[i]))
        elif kwargs.get('params_order'):# order params according to flags
            for i in kwargs.get('params_order'):
                score.append(self.__calc_diff_params(self.p[i], params2[i]))
        elif kwargs.get('flag_ptp_order'):# order params according to flags
            for i, param1 in enumerate(self.p):
                if self.flag[i]:
                    score.append(self.__calc_diff_params(param1, params2[i]))
                else:
                    score.insert(0, self.__calc_diff_params(param1, params2[i]))
        return score

    def check_number_of_params(self):
        if len(self.p) != len(self.p_type):
            exc_m = 'wrong number of parameters given\n' \
                '{0:d} given - {1:d} expected'.format(len(self.p), len(self.p_type))
            raise Exception(exc_m)

    def set_params(self, params, replace=None):
        """sets parameters for an interaction type - len(params) has to be the same as len(p_type)"""
        if self.p:
            if replace:
                self.p = []
            else:
                do_warn("parameters already defined, use replace=True if you want to replace them")
                return
        self.add_params(params)
        self.check_number_of_params()

    def add_params(self, params, flag_check_len=False, flag_len_break=False):
        """adds successively parameters to the interaction type
        len(params) could be one or more, even more than len(p_type) in combination with flag_len"""
        temp_pos = len(self.p)
        end_pos = len(self.p_type) - 1
        if temp_pos == end_pos + 1:
            params = list(params)
            if params:
                do_warn("all parameters already defined")
            return
        if type(params) in (str, int, float):
            params = (params,)
        try:
            if type(params) == unicode: # python2
                params = (params,)
        except:pass
        for i in params:
            try:
                temp = self.p_type[temp_pos](i)
            except ValueError:
                raise Exception(str(self.p_type[temp_pos]) + " expected - interaction parameter: " + str(i))
            self.p.append(temp)
            if temp_pos == end_pos:
                self.p = tuple(self.p)
                if flag_len_break:
                    return
            temp_pos += 1
        if flag_check_len:
            self.check_number_of_params()

    @staticmethod
    def check_param(param, p_type):
        try:
            temp = p_type(param)
        except:
            if p_type == float:
                raise Exception("float expected;", param, "given")
            else:
                raise Exception("integer expected;", param, "given")
        return temp

    def add_atom_types(self, atoms):
        self.atoms = atoms

    def add_p_type_flag(self, code=None):
        self.p_type = []
        self.flag = []
        if code is None:
            code = self.fnc[self.fnc_type]
        h = int(len(code) / 2)
        for i in range(h):
            self.p_type.append(_p_type_code[code[i]])
            self.flag.append(int(code[h + i]))

    @staticmethod
    def check_at_types(at1, at2):
        for i in range(len(at1)):
            if at1[i] != "X":
                if at1[i] != at2[i].a_t.b_id:
                    return False
            else:
                continue
        return True

    def check_atom_types(self, atoms, ratoms):
        if self.check_at_types(self.atoms, atoms):
            return True
        if self.check_at_types(self.atoms, ratoms):
            return True
        return False

    @classmethod
    def find_int_container2write(cls):
        if hasattr(cls, 'int_container2write'):
            return cls.int_container2write
        elif hasattr(cls, 'container2write'):
            return cls.container2write
        else:
            return 'intDB_' + cls.__name__

    @staticmethod
    def _check_atoms_int_match_presence(int_atoms, atoms, **kwargs):
        for at in atoms:
            if at not in int_atoms:
                return False
        return True

    @staticmethod
    def _check_atoms_int_match_ordered(int_atoms, atoms, **kwargs):
        #assert len(int_atoms) == len(atoms)
        allowed_perturbations = kwargs.get('allowed_perturbations')
        if allowed_perturbations is None:
            allowed_perturbations = [list(range(len(int_atoms)))]
            if kwargs.get('atoms_reversed', True):
                rev_pert = allowed_perturbations[0][::-1]
                allowed_perturbations.append(rev_pert)
        for temp_pert in allowed_perturbations:
            flag = True
            for i, ii in enumerate(temp_pert):
                if atoms[ii] != int_atoms[i]:
                    flag = False
                    break
            if flag:
                return True
        return False


class ExclusionType(InteractionType, Defaults):
    """
    exclusion type - but also carries info on 1-4 pair
    # 0 normal wdv
    # 1 excluded
    # 2 pair
    """
#    container2write = 'exclusion_types'
#    int_container2write = 'exclusions'
    container2write = 'exclusions'
    na = 2
    fnc = {None: 'i1'} # check this?
    fnc['gr_fnc'] = 'i1'

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

ExclusionType._add_class_defaults({'gr2gm_fnc': None}, flag_set=True)


class ExclusionPairType(InteractionType, Defaults):
    """
    exclusion type - but also carries info on 1-4 pair
    # 0 normal wdv
    # 1 excluded
    # 2 pair

    # in GROMOS ptp (PERTATOMPAIR block) these would be the codes:
    # 0 excluded
    # 1 normal vdw
    # 2 pair
    """
    container2write = 'excl_pair'
    na = 2
    fnc = {None: 'i1'} # check this?
    fnc['gr_fnc'] = 'i1'
    _code_map = dict(((0, 1), (1, 0), (2, 2))) # this one translates the internal codes to GROMOS PERTATOMPAIR block codes

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

ExclusionPairType._add_class_defaults({'gr2gm_fnc': None}, flag_set=True)


def generate_params_khq(self, khq):
    params = [0] * len(self.flag) # PTP flag
    for i1, i2 in enumerate(self._khq_pos_fnc_map[self.fnc_type]):
        if i2 is not None:
            params[i2] = khq[i1]
    self.add_params(params)

def get_khm(self):
    khq = []
    for i in self._khq_pos_fnc_map[self.fnc_type]:
        if i is None:
            khq.append(None)
        else:
            khq.append(self.p[i])
    return khq

_deg2rad = math.degrees(1)
_deg2rad_2 = _deg2rad**2
_rad2deg = math.radians(1)
_rad2deg_2 = _rad2deg**2


class BondType(InteractionType, Defaults):
    #__slots__ = ()
    container2write = 'bonds'
    na = 2
    fnc = {"1": "ff11", "2": "ff11", "3": "fff111", "4": "fff000", "5": "",
           "6": "ff11", "7": "ff00", "8": "if01", "9": "if01", "10": "ffff1111"}
    fnc['gr_fnc'] = 'fff111'

    _gr_gm_pos_fnc_map = {} # b0, kh, kq
    _gr_gm_pos_fnc_map['1'] = (0, 1, None)
    _gr_gm_pos_fnc_map['2'] = (0, None, 1)
    _gr_gm_pos_fnc_map['3'] = (0, None, None)
    _gr_gm_pos_fnc_map['4'] = (0, None, None)
    _gr_gm_pos_fnc_map['5'] = (None, None, None)
    _gr_gm_pos_fnc_map['6'] = (0, None, None)
    _gr_gm_pos_fnc_map['7'] = (0, None, None)
    _gr_gm_pos_fnc_map['8'] = (None, None, None)
    _gr_gm_pos_fnc_map['9'] = (None, None, None)
    _gr_gm_pos_fnc_map['10'] = (None, None, None)
    _gr_gm_pos_fnc_map['gr_fnc'] = (2, 1, 0)

    generate_params_khq = generate_params_khq
    get_khm = get_khm

    def convert_harm_quad(self, **kwargs):
        b_khq = self.get_gr_gm_params()
        if b_khq[1] is None and b_khq[2] is None:
            b_khq[1] = b_khq[2] = 0
        elif b_khq[1] is None:
            b_khq[1] = b_khq[2] * 2 * b_khq[0] ** 2
        elif b_khq[2] is None:
            b_khq[2] = b_khq[1] / (2 * b_khq[0] ** 2)
        else:
            #assert check_ind_f_params(b_khq[2], b_khq[1] / (2 * b_khq[0] ** 2)), 'k_harm != 2*k_quad * b0**2'
            assert check_if_eq_np(b_khq[2], b_khq[1] / (2 * b_khq[0] ** 2)), 'k_harm != 2*k_quad * b0**2'
        return b_khq

    def gr2gm(self, **kwargs):
        gm_int_type = self.get_gr2gm_int_type(int_code = self.id, **kwargs)
        b_khq = self.convert_harm_quad(**kwargs)
        gm_int_type.generate_params_gr_gm(b_khq)
        return gm_int_type

    def gm2gr(self, **kwargs):
        b_khq = self.convert_harm_quad(**kwargs)
        if b_khq[0] == None:
            b_khq[0] = kwargs.get('b0', 0.1)
        kwargs['int_code'] = kwargs.get('int_code', self.id)
        gr_int_type = self.get_gm2gr_int_type(**kwargs)
        gr_int_type.generate_params_gr_gm(b_khq)
        return gr_int_type

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

BondType._add_class_defaults({'gr2gm_fnc': '2', 'gr2gm_def_pref': 'gb_'}, flag_set=True)


class AngleType(InteractionType, Defaults):
    #__slots__ = ()
    container2write = 'angles'
    na = 3
    fnc = {'1': 'ff11', '2': 'ff11', '3': 'fff000', '4': 'ffff0000', '5': 'ffff1111', '6': 'fffff00000', '8': 'if01'}
    fnc['gr_fnc'] = 'fff111'

    generate_params_khq = generate_params_khq
    get_khm = get_khm

    _gr_gm_pos_fnc_map = {} # a0, kh, kq
    _gr_gm_pos_fnc_map['1'] = (0, 1, None)
    _gr_gm_pos_fnc_map['2'] = (0, None, 1)
    _gr_gm_pos_fnc_map['3'] = (None, None, None)
    _gr_gm_pos_fnc_map['4'] = (None, None, None)
    _gr_gm_pos_fnc_map['5'] = (0, None, None)
    _gr_gm_pos_fnc_map['6'] = (0, None, None)
    _gr_gm_pos_fnc_map['8'] = (None, None, None)
    _gr_gm_pos_fnc_map['10'] = (0, None, None)
    _gr_gm_pos_fnc_map['gr_fnc'] = (2, 1, 0)

    def convert_harm_quad(self, check_cutoff = 0.001, **kwargs):
        def q2h(a_khq, kt):
            def acos_dom_fix(v):
                return (v + 3) % 2 - 1
            temp_fluc_cosa =math.sqrt(kt / a_khq[2])
            temp_cos_a0 = math.cos(math.radians(a_khq[0]))
            temp_a1 = math.degrees(math.acos(acos_dom_fix(temp_cos_a0 + temp_fluc_cosa)))
            temp_fac1 = d_periodic(temp_a1, a_khq[0], 180) ** 2
            #temp_fac1 = (math.degrees(math.acos(acos_dom_fix(temp_cos_a0 + temp_fluc_cosa))) - a_khq[0]) ** 2
            temp_a2 = math.degrees(math.acos(acos_dom_fix(temp_cos_a0 - temp_fluc_cosa)))
            temp_fac2 = d_periodic(temp_a2, a_khq[0], 180) ** 2
            #temp_fac2 = (math.degrees(math.acos(temp_cos_a0 - temp_fluc_cosa)) - a_khq[0]) ** 2
            return 2*kt / (temp_fac1 + temp_fac2)
        def h2q(a_khq, kt):
            temp_fluc_a = math.sqrt(kt / a_khq[1])
            temp_a1 = math.radians(a_khq[0] + temp_fluc_a)
            temp_a2 = math.radians(a_khq[0] - temp_fluc_a)
            temp_cos_a0 = math.cos(math.radians(a_khq[0]))
            return 2*kt / ((math.cos(temp_a1) - temp_cos_a0)**2 + (math.cos(temp_a2) - temp_cos_a0)**2)
        T = kwargs.get('T', 300)
        kt = kb * T
        kwargs['cutoff'] = check_cutoff
        a_khq = self.get_gr_gm_params()
        if a_khq[1] is None and a_khq[2] is None:
            a_khq[1] = a_khq[2] = 0
        elif a_khq[1] is None:
            a_khq[1] = q2h(a_khq, kt)
            #a_khq[1] = a_khq[2] * (math.sin(math.radians(a_khq[0])) * _rad2deg)**2
        elif a_khq[2] is None:
            a_khq[2] = h2q(a_khq, kt)
#            a_khq[2] = a_khq[1] / (math.sin(math.radians(a_khq[0])) * _rad2deg)**2
        else:
            kq = h2q(a_khq, kt)
            assert check_ind_f_params(a_khq[2], kq, **kwargs), self.id + '\tk_quad != h2q(k_harm) with T = ' + str(T)
            #assert check_if_eq_np(a_khq[2], kq, **kwargs), self.id + '\tk_quad != h2q(k_harm) with T = ' + str(T)
        return a_khq

    def gr2gm(self, **kwargs):
        gm_int_type = self.get_gr2gm_int_type(int_code = self.id, **kwargs)
        a_khq = self.convert_harm_quad(**kwargs)
        gm_int_type.generate_params_gr_gm(a_khq)
        return gm_int_type

    def gm2gr(self, **kwargs):
        a_khq = self.convert_harm_quad(**kwargs)
        if a_khq[0] == None:
            a_khq[0] = kwargs.get('a0', 120.)
        kwargs['int_code'] = kwargs.get('int_code', self.id)
        gr_int_type = self.get_gm2gr_int_type(**kwargs)
        gr_int_type.generate_params_gr_gm(a_khq)
        return gr_int_type

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_ordered(int_atoms, atoms, **kwargs)

AngleType._add_class_defaults({'gr2gm_fnc':'2', 'gr2gm_def_pref': 'ga_'}, flag_set=True)


class Dihedral_Check_atoms_int(InteractionType, Defaults):
    @classmethod
    def _check_atoms_int_match1(cls, int_atoms, atoms, **kwargs):
        if len(atoms) == 4:
            atoms = atoms[1:3]
        return cls._check_atoms_int_match_ordered(int_atoms[1:3], atoms)

    @classmethod
    def _check_atoms_int_match2(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_ordered(int_atoms, atoms)

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, v=1, **kwargs):
        return getattr(cls, '_check_atoms_int_match' + str(v))(int_atoms, atoms, **kwargs)


class DihedralType(Dihedral_Check_atoms_int):
    na = 4

    def _gm_check_dih_imp(self):
        if self.__class__ == DihedralType:
            if self.fnc_type in self.__dih_imp:
                self.__class__ = ImproperType
            else:
                self.__class__ = DihedralType

    #"""
    def __init__(self, *args, **kwargs):
        super(DihedralType, self).__init__(*args, **kwargs)
        self._gm_check_dih_imp()
    #"""


fnc_dih = {"1": "ffi110", "2": "ff11", "3": "ffffff111111", "4": "ffi110", "5": "ffff1111", "8": "if01", "9": "ffi110",
           '10':'fffff00000', '11':'ff00'}
fnc_dih['gr_fnc'] = 'ffi111'
dih_gr_gm_pos_fnc_map = {} # da0, k, m, d0_imp, k_imp
dih_gr_gm_pos_fnc_map['1'] = (0, 1, 2, None, None)
dih_gr_gm_pos_fnc_map['2'] = (None, None, None, 0 , 1) # improper dihedral
dih_gr_gm_pos_fnc_map['3'] = (None, None, None, None, None)
dih_gr_gm_pos_fnc_map['4'] = (0, 1, 2, None, None) # same as 1, but for gmx this counts as improper dihedral
dih_gr_gm_pos_fnc_map['5'] = (None, None, None, None, None)
dih_gr_gm_pos_fnc_map['8'] = (None, None, None, None, None)
dih_gr_gm_pos_fnc_map['9'] = (0, 1, 2, None, None)
dih_gr_gm_pos_fnc_map['10'] = (None, None, None, None, None)
dih_gr_gm_pos_fnc_map['11'] = (None, None, None, None, None)
dih_gr_gm_pos_fnc_map['gr_fnc'] = (1, 0, 2, None, None)


class ImproperType(DihedralType, Defaults):
    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

_DihedralType_defs = {'fnc':fnc_dih, 'container2write':'dihedrals', '_gr_gm_pos_fnc_map':dih_gr_gm_pos_fnc_map}
DihedralType._add_defaults(_DihedralType_defs, flag_set=True)
_DihedralType_class_defs = {'gr2gm_fnc':'1', 'gr2gm_def_pref': 'gd_', 'dih_imp':('2', '4')}
#_DihedralType_class_defs = {'gr2gm_fnc':'1', 'gr2gm_def_pref': 'gd_', 'dih_imp':('2')}
_DihedralType_class_defs['gr2gm_def_pref_fnc_type'] = {'2':'gi_'}
_DihedralType_class_defs['gm2gr_int_type_fnc_type'] = {'2':ImproperType, '4':DihedralType}
DihedralType._add_class_defaults(_DihedralType_class_defs, flag_set=True)

fnc_imp = dict(DihedralType.fnc)
fnc_imp['gr_fnc'] = 'ff11'
imp_gr_gm_pos_fnc_map = dict(dih_gr_gm_pos_fnc_map)
imp_gr_gm_pos_fnc_map['gr_fnc'] = (None, None, None, 1, 0)
_ImproperType_defs = {'fnc':fnc_imp, 'container2write':'impropers', '_gr_gm_pos_fnc_map':imp_gr_gm_pos_fnc_map}
_ImproperType_defs['gr2gm_param_trans_fnc'] = (None, None, None, None, lambda x:x * _deg2rad_2)
_ImproperType_defs['gm2gr_param_trans_fnc'] = (None, None, None, None, lambda x:x * _rad2deg_2)
ImproperType._add_defaults(_ImproperType_defs, flag_set=True)
#ImproperType._add_class_defaults({'gr2gm_fnc':'2', 'gr2gm_int_type':DihedralType, 'gr2gm_def_pref': 'gi_'}, flag_set=True)
ImproperType._add_class_defaults({'gr2gm_fnc':'2', 'gr2gm_def_pref': 'gi_'}, flag_set=True)


def convert_se2c612(se):
    s6 = se[0] ** 6
    c6 = s6 * 4 * se[1]
    c12 = c6 * s6
    return c6, c12


def convert_C612_2_se(c612):
    if c612[0]==0:
        return 0,0
    s6 = (c612[1] / c612[0])
    s = s6 ** (1 / 6.)
    e = c612[0] / (4 * s6)
    return s, e


class _vdWTypeBase(InteractionType, Defaults):
    #__slots__ = ('c612', 'se')

    def convert2se(self, c612 = None, replace = None):
        if c612 is None:
            c612 = self.p
        else:
            self.set_params(c612, replace = replace)
        s, e = convert_C612_2_se(c612)
        self.c612 = self.p
        self.se = (s,e)
        return s, e

    def convert2c612(self, se = None, replace = None):
        if se is None:
            se = self.p
        else:
            self.set_params(se, replace = replace)
        self.se = self.p
        c6, c12 = convert_se2c612(se)
        self.c612 = (c6, c12)
        return c6, c12

    def get_radius(self):
        return self.se[0] / 2

    def gm2gr(self, **kwargs):
        kwargs['int_code'] = kwargs.get('int_code', self.id)
        kwargs['atoms'] = kwargs.get('atoms', self.atoms)
        comb_rule = kwargs.get('comb_rule')
        if comb_rule is not None:
            if comb_rule == '1':
                kwargs['params'] = self.p
            else:
                kwargs['params'] = self.convert2c612()
        else:
            kwargs['params'] = self.c612
        gm_int_type = self.get_gm2gr_int_type(**kwargs)
        return gm_int_type


class vdWType(_vdWTypeBase):
    #__slots__ = ()
    container2write = 'vdw_normal'
    fnc = {"1": "ff11", "2": "fff111"}
    fnc['gr_fnc'] = 'ff11'
    na = 2


vdWType._add_class_defaults({'gr2gm_fnc':'1'}, flag_set=True)


class PairType(_vdWTypeBase):
    #__slots__ = ()
    container2write = 'vdw_pairs'
    fnc = {"1": "ff11", "2": "ffff0000"}
    fnc['gr_fnc'] = 'ff11'
    na = 2

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

PairType._add_class_defaults({'gr2gm_fnc':'1'}, flag_set=True)


class PairNBType(InteractionType, Defaults):
    #__slots__ = ()
    fnc = {"1": "ffff0000"}
    na = 2

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class vdWAtomPairType(__VersionCompatibility, tuple):
    container2write = 'vdw'


class MassType(__VersionCompatibility):
    #__slots__ = ('id', 'm', 'at_id')
    container2write = 'm_type'

    def __init__(self, mass_id, mass, atom_name):
        self.id = mass_id
        self.m = float(mass)
        self.at_id = atom_name

    def __str__(self):
        return str(self.id) + ' ' + str(self.m)

    def __repr__(self):
        return str(self.id) + ' ' + str(self.m)



class AtomType(__VersionCompatibility):
    #__slots__ = ('id', 'name', 'vdw', 'c612', 'rules')
    container2write = 'a_type'

    def __init__(self, atom_id, atom_name, c612 = None, rules=None, vdw=None, bond_atom_id = None, format_type = None,
                 **kwargs):  # atom_id = atom number!
        self.id = atom_id
        self.b_id = bond_atom_id # gromacs thing
        self.name = atom_name
        if c612 is None:
            c612 = list()
        if rules is None:
            rules = list()
        if vdw is None:
            vdw = list()
        self.c612 = c612
        if rules is None:
            rules = list()
        self.rules = rules
        self.vdw = vdw
        self.format_type = format_type
        setattr(self, format_type + '_id', atom_id)

    def convert_se_2_c612(self):
        temp_c612 = convert_se2c612(self.vdw)
        self.se_c612 = tuple(self.vdw) + temp_c612

    def convert_c612_2_se(self):
        temp_se = convert_C612_2_se(self.vdw)
        self.se_c612 = temp_se + tuple(self.vdw)

    def add_genborn(self, params):  # defines genborn params
        self.genborn = params

    def __str__(self):
        return str(self.id) + ' ' + str(self.name)

    def __repr__(self):
        return str(self.id) + ' ' + str(self.name)


class BondAtomType(str):
    container2write = 'b_a_type'


class ConstraintType(InteractionType, Defaults):
    container2write = 'constraints'
    na = 2
    fnc = {"1": "f1", "2": "f1"} # 1 - adds exclusions; 2 - doesn't
    fnc['gr_fnc'] = 'f1'

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)

ConstraintType._add_class_defaults({'gr2gm_fnc':'1'}, flag_set=True)


class SettleType(InteractionType):
    na = 1
    container2write = 'settles'
    fnc = {"1": "ff00"}


class VirtualSite2Type(InteractionType):
    na = 3
    container2write = 'VS2'
    fnc = {"1": "f0"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class VirtualSite3Type(InteractionType):
    na = 4
    container2write = 'VS3'
    fnc = {"1": "ff00", "2": "ff00", "3": "ff00", "4": "fff000"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class VirtualSite4Type(InteractionType):
    na = 5
    container2write = 'VS4'
    fnc = {"2": "ff00"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class VirtualSitenType(InteractionType):
    na = 1
    #gmx_parse_format =
    container2write = 'VSn'

    def __init__(self, fnc_type, int_code=None):
        self.f = fnc_type
        self.w = []

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class Position_r(InteractionType):
    na = 1
    container2write = 'posres'
    fnc = {"1": "fff111", "2": "ff00"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class Distance_r(InteractionType):
    """
    distance restraints
    extra gromos attributes: gr_extra_atoms, gr_atom_types; gr_distres_dis_const (in top)

    """
    na = 2
    container2write = 'disres'
    fnc = {"1": "iiiffff0000000"} ####### probably wrong
    fnc['gr_fnc'] = 'ffi000' ############ probably wrong

    def gr2gm(self, **kwargs):
        """
        :param kwargs:
            index - gromacs dist res index - specifies to which group of dist res this one belongs to
            fac - factor with which the force constant will be multiplied (fc_gm = gr_weight / fac)
            r2 - parameter after which the potential energy becomes linear
        :return:
        """
        index_param = kwargs.get('index', 0)
        params = [1, index_param, 1]
        r2 = kwargs.get('r2', 100)
        fac = kwargs.get('fac', 1.)
        if self.p[2] == 0:
            r012 = [self.p[0], self.p[0], r2]
        if self.p[2] == -1:
            r012 = [self.p[0], self.p[0], self.p[0]]
        if self.p[2] == 1:
            r012 = [0, self.p[0], r2]
        params.extend(r012)
        params.append(self.p[1] / fac)
        temp_dr_type = self._get_gr2gm_int_type_class()(fnc_type = '1')
        temp_dr_type.add_params(params)
        return temp_dr_type

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class Dihedral_r(Dihedral_Check_atoms_int):
    na = 4
    container2write = 'dihres'
    fnc = {"1": "ff11"}


class Orientation_r(InteractionType):
    na = 2
    container2write = 'orires'
    fnc = {"1": "ffffff000000"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


class Angle_r(InteractionType):
    na = 4
    container2write = 'angres'
    fnc = {"1": "ffi110"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        allowed_perturbations = []
        for fac1 in (1, -1):
            for fac2 in (1, -1):
                first_pair = [0, 1][::fac1]
                second_pair = [2, 3][::fac2]
                first_pair.extend(second_pair)
                allowed_perturbations.append(first_pair)
        return cls._check_atoms_int_match_ordered(int_atoms, atoms, allowed_perturbations=allowed_perturbations)

class Angle_r_z(InteractionType):
    na = 2
    container2write = 'angresz'
    fnc = {"1": "ffi110"}

    @classmethod
    def _cls_check_atoms_int_match(cls, int_atoms, atoms, **kwargs):
        return cls._check_atoms_int_match_presence(int_atoms, atoms)


# cmap done separately
class cmap(InteractionType, gmCMAPWriter):
    na = 5
    container2write = 'cmap'

    def __init__(self, int_code = None, params = None, fnc_type = '1', atoms = None, grid_ind=False):
        self.fnc_type = fnc_type
        self.grid_ind = grid_ind
        self.p_type = []
        self.flag = []
        if not grid_ind:
            grid_ind = [24, 24]
        for i in range(grid_ind[0]):
            for j in range(grid_ind[0]):
                self.p_type.append(float)
                self.flag.append(0)
        self.int_code = int_code
        self.p_string = []
        self.p = []
        self.atoms = []
        if params:
            self.set_params(params)
        if atoms:
            if hasattr(atoms,'__iter__'):
                self.add_atom(*atoms)
            else:
                self.add_atom(atoms)

def __check_interaction_type_class(item):
    if hasattr(item, '__mro__') and InteractionType in item.__mro__:
        return True


class AvailableInteractionTypes:
    _available_int_types = []

for k in list(locals().keys()):
    if k.endswith('Type') or __check_interaction_type_class(locals()[k]) and not k.startswith('_'):
        if k!= 'InteractionType':
            AvailableInteractionTypes._available_int_types.append(k)
            setattr(AvailableInteractionTypes, k, locals()[k])

__all__ = list(AvailableInteractionTypes._available_int_types)
__all__.extend(('AvailableInteractionTypes', 'InteractionType', 'Interaction', 'DescriptionPart', 'check_if_eq'))

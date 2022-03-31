from SMArt.incl import copy, np, math, combinations, combinations_with_replacement, OrderedDict, \
    defaultdict, DataDumping, Defaults, do_warn, gzip
from SMArt.incl import GeneralContainer
from SMArt.md.incl import *
from SMArt.geometry import rot
from SMArt.graph.incl import GraphDirected
from SMArt.md.gro2gro.g2g import FFg2g, Top_g2g, MolType_g2g, \
                _atomtype_name__element_map, _element__AtomicNum_map, _AtomicNum__element_map
from .gromos.incl import GromosBlockNames, GromosDefaults
from .gromos.io.incl import GromosFile, GromosString, GromosParser, GromosWriter, IFPBlocksParser, IFPBlocksWriter,\
                            BBParsing, grBBAtomWriting, BBWriting, MTBBlocksParser, topBlocksParser, topBlocksWriter, \
                            grTOPAtomWriting, PTP_EDS_BlocksParser, cnfBlocksParser, cnfBlocksWriter,\
                            TrjCnfBlocksParser, TrjCnfBlocksWriter
from .gromacs.io.incl import gmFFParser, gmFFWriter, gmFragmentMoleculeIO, gmTopologyIO, gmConfigurationIO



class FF(AvailableInteractionTypes, DataDumping, IFPBlocksParser, IFPBlocksWriter, gmFFParser, gmFFWriter, FFg2g):
    """Force field class"""
    def __init__(self, parse_from = None, parse_from_file = True, format_type = 'gr', int_db=None, **kwargs):
        if int_db is not None:
            self._int_db = int_db
        self.block_names = GromosBlockNames() #############################################################
        if parse_from:
            self.parse_ff(parse_from, parse_from_file, format_type = format_type, **kwargs)

    def parse_ff(self, parse_from, parse_from_file = True, format_type = 'gr', **kwargs):
        """
        parse force field parameters
        :param parse_from:
        :param parse_from_file:
        :param format_type:
        :param kwargs:
        :return:
        """
        fnc2call = {'gr':self.parse_ifp, 'gm':self.parse_ff_gm}
        assert format_type in fnc2call
        fnc2call[format_type](parse_from, parse_from_file, **kwargs)

    def write_ff(self, f_path = None, format_type = 'gr', **kwargs):
        fnc2call = {'gr':self.write_ifp, 'gm':self.write_ff_itp}
        assert format_type in fnc2call
        fnc2call[format_type](f_path = f_path, **kwargs)

    def get_intDB(self, int_type=None, create=True, **kwargs):
        #gets the interaction database, if int_type given, returns a container for a given interaction type, e.g. BondType
        if not int_type:
            if self._int_db:
                return self.get_container(self._int_db, create=create, db_type=GeneralContainer, **kwargs)
            else:
                return self
        if self._int_db:
            pre_conts = [self._int_db]
        else:
            pre_conts = []
        return self.get_container(int_type, pre_containers=pre_conts, flag_item=True, create=create,
                                  db_type=GeneralContainer, **kwargs)

    @property
    def get_DUM_type(self):
        return self.get_intDB().DUM_type

    def add_a_type(self, atom_type_id, atom_name, vdw=None, rules=None, replace=False, **kwargs):
        """
atom_id - 1,2,3... if False: defined as the next available number
atom_name - CH3,OW...
vdw - 4 LJ parameters + 2 LJpair parameters
rules - [1,2,3] for every atom type
"""
        if vdw is None:vdw = list()
        if rules is None:rules = list()
        self.get_intDB().add2container(self.AtomType(atom_type_id, atom_name, vdw, rules, **kwargs), replace=replace, create=True)
        self.__a_type_missing_rules()

    def get_a_type(self, atom_type_id):
        return self.get_intDB().get_item(atom_type_id, AtomType, allow_not_found=True)

    def __a_type_missing_rules(self):
        # with every new atom type, all atom types need to have one more rule... 1 is added which might be wrong
        n_at = len(self.get_intDB().a_type)
        for at_id, at in self.get_intDB().a_type.items():
            for i in range(n_at - len(at.rules)):
                at.rules.append(1)

    def __generate_index_map(self):
        at_list = list(self.get_intDB().a_type.keys())
        at_index_map = {}
        for i in range(len(at_list)):
            at_index_map[at_list[i]] = i
        return at_index_map

    def _generate_vdw_atpair_gr(self, at_id1, at_id2, replace=False, at_index_map = None):
        if at_index_map is None:
            at_list = list(self.get_intDB().a_type.keys())
            ind1 = at_list.index(at_id1)
            ind2 = at_list.index(at_id2)
        else:
            ind1 = at_index_map[at_id1]
            ind2 = at_index_map[at_id2]
        at1 = self.get_intDB().a_type[at_id1]
        at2 = self.get_intDB().a_type[at_id2]
        c6 = at1.c612[0] * at2.c612[0]
        c12 = at1.c612[at1.rules[ind2]] * at2.c612[at2.rules[ind1]]
        c6_14 = at1.c612[4] * at2.c612[4]
        c12_14 = at1.c612[5] * at2.c612[5]
        vdw_normal = self.vdWType(params=[c6, c12])
        vdw_pair = self.PairType(params=[c6_14, c12_14])
        self.add_vdw(at_id1, at_id2, vdw_normal, vdw_pair, replace=replace)

    def __comb_c612(self, at1, at2):
        c6 = math.sqrt(at1.vdw[0] * at2.vdw[0])
        c12 = math.sqrt(at1.vdw[1] * at2.vdw[1])
        return c6, c12, c6 * self.gm_ff_defaults[3], c12 * self.gm_ff_defaults[3]

    def __comb_se1(self, at1, at2):
        s = (at1.vdw[0] + at2.vdw[0]) / 2.
        e = math.sqrt(at1.vdw[1] * at2.vdw[1])
        return s, e, s, e * self.gm_ff_defaults[3]

    def __comb_se2(self, at1, at2):
        s = math.sqrt(at1.vdw[0] * at2.vdw[0])
        e = math.sqrt(at1.vdw[1] * at2.vdw[1])
        return s, e, s, e * self.gm_ff_defaults[3]

    def _generate_vdw_atpair_gm(self, at_id1, at_id2, replace=False, **kwargs):
        at1 = self.get_intDB().a_type[at_id1]
        at2 = self.get_intDB().a_type[at_id2]
        if self.gm_ff_defaults[1] == '1':
            temp_params = self.__comb_c612(at1, at2)
        elif self.gm_ff_defaults[1] == '2':
            temp_params = self.__comb_se1(at1, at2)
        elif self.gm_ff_defaults[1] == '3':
            temp_params = self.__comb_se2(at1, at2)
        else:
            raise Exception('unknown comb-rule - [1-3] expected, got: ' + str(self.gm_ff_defaults[1]))
        vdw_normal = self.vdWType(fnc_type=self.gm_ff_defaults[0], params=temp_params[:2])
        if self.gm_ff_defaults[2]=='yes':
            vdw_pair = self.PairType(fnc_type = self.__gm_pair_fnc_default, params=temp_params[2:])
        else:
            vdw_pair = None
        self.add_vdw(at_id1, at_id2, vdw_normal, vdw_pair, replace=replace, **kwargs)
        return vdw_normal, vdw_pair

    def __at_pair(self, at1, at2, flag_sort = True):
        if flag_sort:
            if not hasattr(self, 'at_index_map'):
                self.at_index_map = self.__generate_index_map()
            if self.at_index_map[at1.id] > self.at_index_map[at2.id]:
                return tuple((at2, at1))
        return tuple((at1, at2))

    def generate_vdw(self, at_ids=False, replace=False, format_type = None, **kwargs):
        """
        generates vdw interaction parameters and sets them in self.vdw
        :param at_ids:
        :param replace:
        :param format_type:
        :param kwargs:
            flag_generate_vdw_c612_se = True
        :return:
        """
        if not at_ids:
            at_ids = list(self.get_intDB().a_type)
        if format_type is None:
            format_type = self.get_intDB().a_type[at_ids[0]].format_type
        if format_type not in ('gr', 'gm'):
            raise Exception('format_type is None, but it should be either gr or gm')
#        if format_type == 'gr':
        if 1:
            self.__a_type_missing_rules()
        if format_type == 'gm':
            self.gm_ff_defaults[3] = float(self.gm_ff_defaults[3])
        at_index_map = self.__generate_index_map()
        temp_gen_fnc ={'gr':self._generate_vdw_atpair_gr, 'gm':self._generate_vdw_atpair_gm}[format_type]
        for i in combinations_with_replacement(at_ids, 2):
            _ = temp_gen_fnc(*i, replace=replace, at_index_map = at_index_map)
        if format_type == 'gm':
            gm_vdw_normal = self.get_intDB(vdWType, create= False, cont_pref = self._gm_cont_pref, allow_not_found = 1)
            if gm_vdw_normal:
                temp_vdw_list = []
                for vdw_normal in gm_vdw_normal:
                    at_pair = self.__at_pair(*vdw_normal.atoms)
                    if at_pair[0].id in at_ids and at_pair[1].id in at_ids:
                        if at_pair not in temp_vdw_list:
                            temp_vdw_list.append(at_pair)
                        else:
                            do_warn('vdw_normal parameters already defined ' + str(at_pair))
                        self.get_intDB().vdw[at_pair][0].p = vdw_normal.p
            gm_vdw_pairs = self.get_intDB(PairType, create= False, cont_pref = self._gm_cont_pref, allow_not_found = 1)
            if gm_vdw_pairs:
                temp_vdw_list = []
                for vdw_pair in gm_vdw_pairs:
                    at_pair = self.__at_pair(*vdw_pair.atoms)
                    if at_pair[0].id in at_ids and at_pair[1].id in at_ids:
                        if at_pair not in temp_vdw_list:
                            temp_vdw_list.append(at_pair)
                        else:
                            do_warn('vdw_pair parameters already defined ' + str(at_pair))
                        temp_pair_type = self.get_intDB().vdw[at_pair][1]
                        if temp_pair_type:
                            temp_pair_type.p = vdw_pair.p
                        else:
                            vdw_atpair = self.vdWAtomPairType((self.get_intDB().vdw[at_pair][0], vdw_pair))
                            self.get_intDB().vdw[at_pair] = vdw_atpair
        if kwargs.get('flag_generate_vdw_c612_se', True):
            kwargs['format_type'] = format_type
            self.generate_vdw_c612_se(**kwargs)

    def generate_vdw_c612_se(self, format_type = None, **kwargs):
        if format_type is None:
            for at in self.get_intDB().a_type.values():
                break
            format_type = at.format_type
        if format_type not in ('gr', 'gm'):
            raise Exception('format_type is None, but it should be either gr or gm')
        if format_type == 'gm' and self.gm_ff_defaults[1] != '1':
            for vdw_pair in self.get_intDB().vdw.values():
                for i in vdw_pair:
                    if i: _ = i.convert2c612()
        else:
            for vdw_pair in self.get_intDB().vdw.values():
                for i in vdw_pair:
                    if i: _ = i.convert2se()


    def add_vdw(self, at_id1, at_id2, vdw_normal = None, vdw_pair = None, replace=False, **kwargs):
        """adds vdw parameters for a pair of atoms in self.vdw - key is a tuple from self.__at_pair"""
        for at_id in (at_id1, at_id2):
            if at_id not in self.get_intDB().a_type:
                raise Exception('atom type ' + str(at_id) + ' not defined')
        at_pair = self.__at_pair(self.get_intDB().a_type[at_id1], self.get_intDB().a_type[at_id2])
        if vdw_normal:
            vdw_normal.id = at_pair
            vdw_normal.add_atom(*at_pair)
        if vdw_pair:
            vdw_pair.id = at_pair
            vdw_pair.add_atom(*at_pair)
        vdw_atpair = self.vdWAtomPairType((vdw_normal, vdw_pair))
        self.get_intDB().add2container(vdw_atpair, replace=replace, item_id=at_pair, create=True, **kwargs)

    def get_vdw(self, at_id1, at_id2):
        """gets vdw parameters (both normal and pair) for a pair of atom types"""
        if not hasattr(self.get_intDB(),'vdw'):
            self.generate_vdw()
        at1 = self.get_intDB().a_type[at_id1]
        at2 = self.get_intDB().a_type[at_id2]
        at_pair = self.__at_pair(at1, at2)
        if at_pair in self.get_intDB().vdw:
            return self.get_intDB().vdw[at_pair]
        else:
            raise Exception("vdw parameters not defined - atoms: " + str(at_id1) + " " + str(at_id2))

    def add_m_type(self, mass_id, mass, atom_name, replace=False):
        self.get_intDB().add2container(self.MassType(mass_id, mass, atom_name), replace=replace)

# additional functions
    def atom_radius(self, at_type, at_type2=None):
        if at_type2 is None:
            at_type2 = at_type
        return self.get_vdw(at_type, at_type2)[0].get_radius()

# converter part - matching the atom types with the mass types
    @staticmethod
    def match_at_m_name(at_name, m_names):
        if "DUM" in at_name:
            return 0.0
        if at_name in m_names:
            return m_names[at_name]
        for i in range(1, len(at_name)):
            if at_name[:-i] in m_names:
                return m_names[at_name[:-i]]
        return False

    def find_matches_at_m_types(self):
        matches = OrderedDict()
        m_names = {}
        for i in self.get_intDB().m_type:
            m_names[self.get_intDB().m_type[i].at_id] = self.get_intDB().m_type[i].m
        for i in self.get_intDB().a_type:
            for j in self.get_intDB().a_type[i].name.split(","):
                temp_match = self.match_at_m_name(j, m_names)
                if temp_match:
                    matches[j] = [i, temp_match]
        self.matches = matches

    def find_imp_pairs(self):################################################### sort of fine
        self.get_intDB().imp_tetra_pair = {}
        self.get_intDB().imp_flat_pair = {}
#        for imp_p_id in combinations(self.get_intDB().impropers, 2):
        imp_cont = self.get_intDB().get_container('gr_impropers', allow_not_found = True)
        if imp_cont:
            for imp_p in combinations(imp_cont.values(), 2):
                if imp_p[0].p[0] == imp_p[1].p[0] and imp_p[0].p[1] == -imp_p[1].p[1] and 34 < abs(imp_p[0].p[1]) < 37:
                    self.get_intDB().imp_tetra_pair[imp_p[0]] = imp_p[1]
                    self.get_intDB().imp_tetra_pair[imp_p[1]] = imp_p[0]
                if imp_p[0].p[0] == imp_p[1].p[0] and abs(abs(imp_p[0].p[1] - imp_p[1].p[1]) - 180) < 0.1:
                    if abs(imp_p[0].p[1]) < 0.1 or abs(imp_p[1].p[1]) < 0.1:
                        self.get_intDB().imp_flat_pair[imp_p[0]] = imp_p[1]
                        self.get_intDB().imp_flat_pair[imp_p[1]] = imp_p[0]

    def __check_atom_int_types_match_pert(self, atoms, int_type_atoms, perturbation):
        flag = False
        if hasattr(self.get_intDB(), 'b_a_type'):
            flag = True
            for i in range(len(atoms)):
                if int_type_atoms[i]!=self._b_a_atom_wildcard and atoms[perturbation[i]].b_id != int_type_atoms[i]:
                    flag = False
        if flag:
            return True
        for i in range(len(atoms)):
            if int_type_atoms[i] != self._b_a_atom_wildcard and atoms[perturbation[i]] != int_type_atoms[i]:
                return False
        return True

    def __check_atom_int_types_match(self, atoms, int_type, allowed_perturbations, atoms_reversed, **kwargs):
        int_type_atoms = list(int_type.atoms)
        if kwargs.get('flag_dih_2_atoms', True):
            if isinstance(int_type, DihedralType) and len(int_type.atoms)==int_type.na - 2:
                int_type_atoms.insert(0, self._b_a_atom_wildcard)
                int_type_atoms.append(self._b_a_atom_wildcard)
        if len(atoms)!=len(int_type_atoms) and len(int_type_atoms)!=int_type.na:return False
        if allowed_perturbations is None:
            allowed_perturbations = [list(range(int_type.na))]
        if atoms_reversed:
            allowed_perturbations.append(list(reversed(range(int_type.na))))
        for at_pert in allowed_perturbations:
            if self.__check_atom_int_types_match_pert(atoms, int_type_atoms, at_pert):
                return True
        return False

    def __find_exact_atom_perturbation(self, atoms, int_type_mathes):
        new_int_type_matches = []
        for it in int_type_mathes:
            if atoms == it.atoms:
                new_int_type_matches.append(it)
        return new_int_type_matches

    def find_interaction_type(self, atoms, int_type_klass, fnc_type = None, allowed_perturbations = None,
                              atoms_reversed = True, use_vdw = True, flag_fnc_type = True, **kwargs):
        """
        :param atoms: list of atom types!
        :param int_type_klass: BondTypes for instance
        :param fnc_type: e.g. '1' - see int_type_klass.fnc
        :param allowed_perturbations:
        :param atoms_reversed:
        :param kwargs:
            int_container
            allow_not_found
            allow_multiple_matches
        :return:
        """
        atom_types = []
        for at in atoms:
            if isinstance(at, AtomType):
                atom_types.append(at)
            else:
                atom_types.append(at.a_type)
        atoms = atom_types
        if use_vdw and int_type_klass in (vdWType,PairType):
            return [self.get_vdw(*[at.id for at in atoms])[(vdWType,PairType).index(int_type_klass)]]
        kwargs['cont_pref'] = kwargs.get('cont_pref', self._gm_cont_pref)
        int_type_container = kwargs.get('int_type_container', self.get_intDB(int_type_klass, **kwargs))
        int_type_matches =[]
        if not int_type_container:return int_type_matches
        for temp_int_type in int_type_container:
            #if temp_int_type.fnc_type == fnc_type:
            if not flag_fnc_type or temp_int_type.fnc_type == fnc_type:
                if (self.__check_atom_int_types_match(atoms, temp_int_type, allowed_perturbations, atoms_reversed,
                                                      **kwargs)):
                    int_type_matches.append(temp_int_type)
        if not int_type_matches:
            if kwargs.get('allow_not_found'):return
            raise Exception("parameters missing; int_type:", int_type_klass, "atoms:", atoms, 'fnc_type:', fnc_type)
        else:
            if not kwargs.get('allow_multiple_matches') and len(int_type_matches)>1:
                raise Exception("multiple matches; int_type:", int_type_klass, "atoms:", atoms, 'fnc_type:', fnc_type)
        return int_type_matches

    def check_eq_int_type(self, int_type, **kwargs):
        kwargs['flag_item'] = 1
        kwargs['allow_not_found'] = 1
        temp_cont = self.get_intDB().get_container(int_type, **kwargs)
        if temp_cont:
            if isinstance(temp_cont, dict):
                temp_cont_list = temp_cont.values()
            else:
                temp_cont_list = temp_cont
            for temp_int_type in temp_cont_list:
                if int_type.fnc_type==temp_int_type.fnc_type and int_type.check_eq_params(temp_int_type.p):
                    return temp_int_type
        return False

_FF_defs = {}
_FF_defs['_int_db'] = ''
_FF_defs['_b_a_atom_wildcard'] = 'X'
_FF_defs['__gm_pair_fnc_default'] = '1'
FF._add_defaults(_FF_defs, flag_set=1)


class InteractionContainer(GeneralContainer):
    def get_interaction_containers(self, EP_cont_exclude = False, **kwargs):
        """
        finds all interaction containers (e.g. bonds, angels, etc.) - generator
        :param EP_cont_exclude:
        :param kwargs:
             get_container_name
             flag_list_only (search for lists only - e.g. skip dictionaries) - False by default
        :return:
        """
        for container_name in getattr(self, '_containers', []):
            temp_cont = self.get_container(container_name)
            if temp_cont:
                if isinstance(temp_cont, list) and isinstance(temp_cont[0], Interaction):
                    if EP_cont_exclude and temp_cont[0].int_type in (PairType, ExclusionType):
                        continue
                    if kwargs.get('get_container_name'):
                        yield temp_cont, container_name
                    else:
                        yield temp_cont
                if not kwargs.get('flag_list_only', True):
                    if not isinstance(temp_cont, list):
                        try:
                            for temp_int in temp_cont:
                                if isinstance(temp_int, Interaction):
                                    if kwargs.get('get_container_name'):
                                        yield temp_cont, container_name
                                    else:
                                        yield temp_cont
                                break
                        except:pass
                        try:
                            for temp_int in temp_cont.values():
                                if isinstance(temp_int, Interaction):
                                    if kwargs.get('get_container_name'):
                                        yield temp_cont.values(), container_name
                                    else:
                                        yield temp_cont.values()
                                break
                        except:pass

def _get_interaction_containers(self, EP_cont_exclude = False, **kwargs):
    for container_name in getattr(self, '_containers', []):
        temp_cont = self.get_container(container_name)
        if temp_cont:
            if isinstance(temp_cont, list) and isinstance(temp_cont[0], Interaction):
                if EP_cont_exclude and temp_cont[0].int_type in (PairType, ExclusionType):
                    continue
                if kwargs.get('get_container_name'):
                    yield temp_cont, container_name
                else:
                    yield temp_cont
            if not kwargs.get('flag_list_only', True):
                if not isinstance(temp_cont, list):
                    try:
                        for temp_int in temp_cont:
                            if isinstance(temp_int, Interaction):
                                if kwargs.get('get_container_name'):
                                    yield temp_cont, container_name
                                else:
                                    yield temp_cont
                            break
                    except:pass
                    try:
                        for temp_int in temp_cont.values():
                            if isinstance(temp_int, Interaction):
                                if kwargs.get('get_container_name'):
                                    yield temp_cont.values(), container_name
                                else:
                                    yield temp_cont.values()
                            break
                    except:pass


class GeneralAtom(GeneralContainer):
    container2write = 'atoms'
    #__slots__ = ('id', 'name', 'a_type', 'p_ch', 'm')

    def __init__(self, atom_id=None, **kwargs):
        self.id = atom_id
        self.name = kwargs.get('name', 'ANN')

    def gr_get_element(self):
        if self.a_type.name == 'P,SI':
            if round(self.m_type.m) == 31:
                return _atomtype_name__element_map['P']
            else:
                return _atomtype_name__element_map['SI']
        else:
            return _atomtype_name__element_map[self.a_type.name]

    @property
    def __prf(self):
        return str(self.id) + ' ' + str(self.name)

    def __str__(self):
        return self.__prf

    def __repr__(self):
        return self.__prf

    def add_pair(self, *atoms, **kwargs):
        for at in atoms:
            self.add2container(at, container_name='p_l', db_type=list, create=True, **kwargs)
            pass

    def add_excl(self, *atoms, **kwargs):
        for at in atoms:
            self.add2container(at, container_name='e_l', db_type=list, create=True, **kwargs)
            pass

    def add_atom_state(self, atom_type, p_ch, mass):
        self.a_type_states.append(atom_type)
        self.p_ch_states.append(p_ch)
        self.m_states.append(mass)

    def _generate_self_state(self):
        self.add_atom_state(self.a_type, self.p_ch, self.m)

    get_interaction_containers = _get_interaction_containers

class BBAtom(GeneralAtom, grBBAtomWriting):
    def __init__(self, atom_id=None, **kwargs):
        self.id = atom_id
        self.flag_bb = False
        self.e_l = []
        self.flag_excl = True
        self.name = kwargs.get('name', 'ANN')


class TopAtom(GeneralAtom, grTOPAtomWriting):
    def __init__(self, atom_id=None, **kwargs):
        self.id = atom_id
        self.name = kwargs.get('name', 'ANN')
        self.a_type = kwargs.get('a_type', None)
        self.p_ch = kwargs.get('p_ch', 0.0)
        self.m = kwargs.get('m', 0.0)
        self.a_type_states = []
        self.p_ch_states = []
        self.m_states = []
        self.mk_cg = 0
        self.p_l = []
        self.e_l = []
#        self.EP_l = set()
#        self.st = []

    @classmethod
    def __create_new_atom(cls, **kwargs):
        return cls(**kwargs)


class ChargeGroup():
    """net_ch - net charge
    atoms - list of atoms"""
    container2write = 'cg'

    def __init__(self, *args, **kwargs):
        gm_num = kwargs.get('gm_num', None)
        if gm_num is not None:
            self.n = gm_num
        self.net_ch = 0.0
        self.atoms = []

    def add_atom(self, atom, **kwargs):
        """adds an atom and updates the net_ch"""
        self.atoms.append(atom)
        self.net_ch += atom.p_ch
        atom.cg = self

    def update(self):
        """updates the net_ch - if p_ch of an atom changes"""
        self.net_ch = 0
        for at in self.atoms:
            self.net_ch += at.p_ch
            at.mk_cg = 0
        at.mk_cg = 1


class GeneralTopology(DataDumping, InteractionContainer, GraphDirected):
    """general class for any type of topology, including BB, molecule or whole system
    contains atoms, interactions and other attributes"""

    def __just_for_readability(self):
        self.atoms = OrderedDict()
        self.residues = OrderedDict()
        self.excl_pair = None

    Atom = TopAtom
    ChargeGroup = ChargeGroup
    Interaction = Interaction

    @staticmethod
    def _get_atom_kwargs(**kwargs):
        a_type_kwargs = dict(kwargs)
        a_type_kwargs.update(dict(ND_states = None, states = 'a_type_states'))
        m_kwargs = dict(kwargs)
        m_kwargs.update(dict(skip_None_states = True, ND_states = 'ND_m_states', states = 'm_states'))
        pch_kwargs = dict(kwargs)
        pch_kwargs.update(dict(skip_None_states = False, ND_states = 'ND_pch_states', states = 'p_ch_states'))
        return a_type_kwargs, m_kwargs, pch_kwargs

    def get_atom_ptp_states(self, at, a_type_ptp, m_ptp, p_ch_ptp, **kwargs):
        flag_ptp = a_type_ptp or m_ptp or p_ch_ptp
        if not a_type_ptp:
            a_type = self._get_state(at, 'a_type', 'a_type_states', **kwargs)
            a_type_ptp = [a_type] * 2
        if not m_ptp:
            m = self._get_state(at, 'm', 'm_states', **kwargs)
            m_ptp = [m]*2
        if not p_ch_ptp:
            p_ch = self._get_state(at, 'p_ch', 'p_ch_states', **kwargs)
            p_ch_ptp = [p_ch] * 2
        return a_type_ptp, m_ptp, p_ch_ptp, flag_ptp


    def _get_top_other_state(self, top_state = None, other_state = None, **kwargs):
        if top_state is None:
            top_state = getattr(self, 'top_state', 0)
        if other_state is None:
            other_state = getattr(self, 'other_state', 1)
        return top_state, other_state

    def get_ptp_states(self, at_int, top_state = None, other_state = None, **kwargs):
        skip_None_states = kwargs.get('skip_None_states', True)
        top_state, other_state = self._get_top_other_state(top_state, other_state)
        if 'ND_states' in kwargs:
            if kwargs.get('ND_states'):
                ND_states = getattr(at_int, kwargs.get('ND_states'), None)
            else:
                ND_states = None
        else:
            ND_states = getattr(at_int, kwargs.get('ND_states', 'ND_states'), None)
        if ND_states:
            flags = [False, False]
            states = [None, None]
            for ND_state in ND_states:
                for i in range(2):
                    if (top_state, other_state)[i] in ND_states[ND_state]:
                        if skip_None_states and ND_state is None:continue
                        states[i] = ND_state
                        flags[i] = True
            if sum(flags)==2 and states[0]!=states[1]:
                return states
            else:
                return
        states = getattr(at_int, kwargs.get('states', 'states'), None)
        try:
            if skip_None_states and (states[top_state] is None or states[other_state] is None):
                return
            if hasattr(states[top_state], 'p') and hasattr(states[other_state], 'p'):
                if states[top_state].check_eq_params(states[other_state].p, **kwargs):
                    return
                else:
                    return states[top_state], states[other_state]
            else:
                if check_if_eq(states[top_state], states[other_state], **kwargs):
                    return
                else:
                    return states[top_state], states[other_state]
        except:
            return

    def _get_state(self, at_int, state_attr = 'state', states_attr = 'states', **kwargs):
        if state_attr and hasattr(at_int, state_attr):
            return getattr(at_int, state_attr)
        return self.get_state(getattr(at_int, states_attr), **kwargs)

    @staticmethod
    def _check_if_physical_state(state):
        return state not in (None, Dummy)

    def get_state(self, states, top_state = None, other_state = None, find_other_state = False, **kwargs):
        """
        :param states: list of n-states (e.g. at.m_states or interaction.states)
        :param top_state: state to get (by default 0)
        :param other_state: the other state for ptp
        :param find_other_state: find the other state for ptp or EDS ('any', 1, -1)
                'any' - any non-None; 1 - first non-None from top_state; -1 first non-None from top_state with the step -1
        :param kwargs:
        :return:
        """
        if top_state is None:
            top_state = getattr(self, 'top_state', 0)
        res_state = states[top_state]
        if self._check_if_physical_state(res_state):
            return res_state
        elif self._check_if_physical_state(other_state):
            return states[other_state]
        if find_other_state:
            assert find_other_state in ('any', 1, -1)
            if find_other_state == 'any':
                for st in states:
                    if self._check_if_physical_state(st):
                        return st
            if find_other_state == 1:
                for st_i in range(top_state + 1, top_state + len(states)):
                    st = states[st_i % len(states)]
                    if self._check_if_physical_state(st):
                        return st
            if find_other_state == -1:
                for st_i in range(top_state - 1, -1, -1):
                    st = states[st_i]
                    if self._check_if_physical_state(st):
                        return st
                for st_i in range(top_state + 1, len(states)):
                    st = states[st_i]
                    if self._check_if_physical_state(st):
                        return st
        return res_state

    def set_top_state(self, top_state = 0, **kwargs):
        """
        :param top_state:
        :param kwargs:
            flag_EDS_mass will set the mass as the average of all matched atoms
            otherwise - see  get_state function (e.g. other_state)
        :return:
        """
        self.top_state = top_state
        if kwargs.get('other_state'):
            self.other_state = kwargs['other_state']
        for at in self.get_atoms():
            try:
                if kwargs.get('flag_EDS_mass'):
                    temp_masses = []
                    temp_w = []
                    for m in at.ND_m_states:
                        if m is not None:
                            temp_masses.append(m)
                            temp_w.append(len(at.ND_m_states[m]))
                    at.m = np.average(temp_masses, weights=temp_w)
                else:
                    at.m = self.get_state(at.m_states, **kwargs)
            except:pass
            try:at.p_ch = self.get_state(at.p_ch_states, **kwargs)
            except:pass
            try:at.a_type = self.get_state(at.a_type_states, **kwargs)
            except:pass
        for int_cont in self.get_interaction_containers(EP_cont_exclude=False, flag_list_only = False):
            for temp_int in int_cont:
                temp_int.state = self.get_state(temp_int.states, **kwargs)

    def add_atom(self, atom_id=None, **kwargs):
        if atom_id is None:
            self.find_next_code(self.Atom, flag_class=True)
        temp_at = self.Atom(atom_id=atom_id)
        self.add2container(temp_at, **kwargs)
        return temp_at

    def get_atoms(self):
        if dict in self.atoms.__class__.__mro__:
            return self.atoms.values()
        else:
            return self.atoms

    def get_residues(self):
        if dict in self.residues.__class__.__mro__:
            return self.residues.values()
        else:
            return self.residues

    def sort_atoms(self, attrib = 'id', **kwargs):
        self.sort_container(GeneralAtom.container2write, attrib=attrib, **kwargs)

    def create_adj_list(self, containers=('bonds',), adj_type=dict, **kwargs):
        """
        :param containers: ('bonds',)
        :param kwargs:
            adj_type                dict, defaultdict (not so good for interactive use), OrderedDict...
            adj_format_new_vertex   set; list...
            set_adj_format_type     set; tuple...
        :return:
        """
        self.adj = defaultdict(set)
        for container in containers:
            temp_cont = self.get_container(container, allow_not_found = True)
            if temp_cont:
                for temp_interaction in temp_cont:
                    self.add_edge(temp_interaction.atoms)
        for at in self.get_atoms():
            if at not in self.adj:
                self.adj[at]
        kwargs['adj_type'] = adj_type
        self._set_adj_type_format(**kwargs)
        self.set_adj_type_format(**kwargs)

    def gen_graph(self, containers = ('bonds',)):
        self.create_adj_list(containers)

    def add_exclusions_neigh(self, nexcl=2):
        self.create_adj_list()
        for at in self.get_atoms():
            for temp_at in self.BFS_d(at, nexcl):
                if at!=temp_at:
                    #at.add_excl(temp_at, replace = -1)
                    temp_interaction = Interaction(ExclusionPairType, atoms=(at, temp_at))
                    temp_interaction.add_state(fnc_type = None, params = (True,))
                    self.add2container(temp_interaction, create=True, item_id=frozenset((at, temp_at)), replace = -1)
                    self.add_atom_pair2EP_l(at, temp_at)

    @property
    def EP_l(self):
        if not hasattr(self, '_EP_l'):
            self._EP_l = defaultdict(set)
        return self._EP_l

    def get_excl_pair(self):
        return self.get_container('excl_pair', create=True)

    def add_atom_pair2EP_l(self, at, *atoms, **kwargs):
        for at2 in atoms:
            self.EP_l[at].add(at2)
            self.EP_l[at2].add(at)

    def get_HH(self):
        """get hydrogen and heavy atoms"""
        HH = ([], [])
        for at in self.get_atoms():
            if at.m > 1.5:
                HH[1].append(at)
            else:
                HH[0].append(at)
        return HH

    def improper_type(self, imp):
        bond_c = []
        for at in imp.atoms:
            c = 0
            for bond_at in self.adj[at]:
                if bond_at in imp.atoms:
                    c += 1
            bond_c.append(c)
        return bond_c

    #get_interaction_containers = _get_interaction_containers

    def reduce(self, atom_ids = None, **kwargs):
        """
        :param atom_ids:
        :param kwargs:
            flag_inplace
        :return:
        """
        if kwargs.get('flag_inplace'):
            new_obj = self
        else:
            new_obj = copy.deepcopy(self)
        if atom_ids is None:
            if hasattr(self, 'atoms'):
                atom_ids = self.atoms
            else:
                return new_obj
        at_id2del = []
        for at_id in new_obj.atoms:
            if at_id not in atom_ids:
                at_id2del.append(at_id)
        at_id2del = set(at_id2del)
        at2del = set(new_obj.atoms[at_id] for at_id in at_id2del)
        for int_cont in new_obj.get_interaction_containers():
            int2del = []
            for i, temp_int in enumerate(int_cont):
                for temp_at in temp_int.atoms:
                    if temp_at in at2del:
                        int2del.append(i)
                        break
            for i in reversed(int2del):
                del(int_cont[i])
        excl_pair = new_obj.get_excl_pair()
        for at_id in at_id2del:
            at = new_obj.atoms[at_id]
            temp_res = getattr(at, 'res')
            if temp_res:
                temp_res.atoms.remove(at)
            temp_cg = getattr(at, 'cg')
            if temp_cg:
                temp_cg.atoms.remove(at)
            del(new_obj.atoms[at_id])
            # EP_l
            if at in new_obj.EP_l:
                for at2 in new_obj.EP_l[at]:
                    temp_pair = frozenset({at, at2})
                    if temp_pair in excl_pair:
                        del(excl_pair[temp_pair])
                del(new_obj.EP_l[at])
        # EP_l of atoms to keep
        for at in new_obj.get_atoms():
            if at in new_obj.EP_l:
                EP_at2del = []
                for at2 in new_obj.EP_l[at] & at2del:
                    temp_pair = frozenset({at, at2})
                    if temp_pair in excl_pair:
                        del(excl_pair[temp_pair])
                new_obj.EP_l[at] = new_obj.EP_l[at] - at2del
        residues = getattr(new_obj, 'residues', {})
        res2del = []
        for res_id in residues:
            res = residues[res_id]
            if not res.atoms:
                res2del.append(res_id)
        for res_id in res2del:
            del(residues[res_id])
        charge_groups = getattr(new_obj, 'cg', [])
        cg2del = []
        for i, cg in enumerate(charge_groups):
            if not cg.atoms:
                cg2del.append(i)
        for i in reversed(cg2del):
            del(charge_groups[i])
        for cg in charge_groups:
            cg.update()
        return new_obj

    def add_interactions2atoms_dict(self, EP_cont_exclude = True, **kwargs):
        atom_interactions = self.get_container('atom_interactions', create=True)
        for int_container in self.get_interaction_containers(EP_cont_exclude = EP_cont_exclude):
            for temp_int in int_container:
                for at in temp_int.atoms:
                    if at not in atom_interactions:
                        atom_interactions[at] = InteractionContainer()
                    atom_interactions[at].add2container(temp_int, create=True, db_type=list, **kwargs)


    def get_bb(self, bb_atom_names=('C', 'N', 'CA')):
        bb_at = []
        for at_id in self.atoms:
            if self.atoms[at_id].name in bb_atom_names:
                bb_at.append(at_id)
        return bb_at

    def find_interactions(self, atoms, container2search, **kwargs):
        int_matches = []
        for temp_interaction in container2search:
            if temp_interaction.check_atoms_int_match(atoms, **kwargs):
                int_matches.append(temp_interaction)
        if not int_matches:
            if kwargs.get('allow_not_found'):
                return int_matches
            raise Exception("interaction not found, atoms:", atoms)
        else:
            if not kwargs.get('allow_multiple_matches') and len(int_matches)>1:
                raise Exception("multiple matches; atoms:", atoms)
        return int_matches

    def _generate_atom_excl_pair_list(self, **kwargs):
        e_l = {}
        p_l = {}
        for at in self.get_atoms():
            e_l[at] = []
            p_l[at] = []
        atom_order_map = kwargs.get('atom_order_map', self._get_atom_order_map())
        for atom_pair in self.get_excl_pair():
            temp_EP = self.excl_pair[atom_pair]
            at1, at2 = atom_pair
            if hasattr(temp_EP, 'state'):
                temp_state = temp_EP.state
            else:
                temp_state = self.get_state(temp_EP.states, **kwargs)
            if temp_state.p == (1,):
                if atom_order_map[at1] < atom_order_map[at2]:
                    e_l[at1].append(at2)
                else:
                    e_l[at2].append(at1)
            elif temp_state.p == (2,):
                if atom_order_map[at1] < atom_order_map[at2]:
                    p_l[at1].append(at2)
                else:
                    p_l[at2].append(at1)
        for at in self.get_atoms():
            e_l[at].sort(key=lambda item:atom_order_map[item])
            p_l[at].sort(key=lambda item:atom_order_map[item])
        return e_l, p_l

    def _gr_get_excl_pairs(self, at, **kwargs):
        atom_order_map = kwargs.get('atom_order_map', self._get_atom_order_map())
        e_l = []
        p_l = []
        for at2 in self.EP_l[at]:
            if atom_order_map[at] < atom_order_map[at2]:
                atom_pair = frozenset((at, at2))
                temp_EP = self.excl_pair[atom_pair]
                if hasattr(temp_EP, 'state'):
                    temp_state = temp_EP.state
                else:
                    temp_state = self.get_state(temp_EP.states, **kwargs)
                if temp_state.p == (1,):
                    e_l.append(at2)
                if temp_state.p == (2,):
                    p_l.append(at2)
        e_l.sort(key=lambda item:atom_order_map[item])
        p_l.sort(key=lambda item:atom_order_map[item])
        return e_l, p_l

    def _get_atom_order_map(self):
        return dict((at,i) for i, at in enumerate(self.get_atoms()))


class BuildingBlock(BBParsing, BBWriting, GeneralTopology):
    """Building block class"""
    Atom = BBAtom

    container2write = 'bb'

    def __init__(self, bb_id=None, parse_from=None, ff=None, **kwargs):
        self.id = bb_id
        if parse_from:
            self.parse_bb(parse_from, ff, **kwargs)

    def parse_bb(self, parse_from, parse_from_file=True, ff = None):
        self._parse_gr(parse_from, parse_from_file=parse_from_file, ff = ff)

    def write_bb(self, gs = None, get_str = False):
        if gs is None:
            gs = GromosString("")
        blocks2write = ('MTBUILDBLSOLUTE',)
        return self.write_gromos_format(gs, *blocks2write, get_str = get_str)


class TopBBdb(Defaults):
    """setting default for _int_db for TOP and BBdb (MTB/RTP) classes"""
    pass

_TopBBdb_defs = {}
_TopBBdb_defs['_int_db'] = 'ff'
TopBBdb._add_defaults(_TopBBdb_defs, flag_set=1)


class BBdb(TopBBdb, FF, MTBBlocksParser, GromosWriter):
    """Building Block database class"""
    BuildingBlock = BuildingBlock

    def __init__(self, mtb_file=None, ifp_file=None, int_db = None):
        if int_db is not None:
            self._int_db = int_db
        if ifp_file:
            self._parse_gr(ifp_file)
        if mtb_file:
            self._parse_gr(mtb_file)

    def add_bb(self, bb_id=None, parse_from = None, parse_from_file = True, **kwargs):
        """adds building blocks to the MTB; see add2container for the variables"""
        if parse_from:
            bb = self.BuildingBlock(bb_id)
            bb.parse_bb(parse_from, parse_from_file, ff = self)
            self.add2container(bb, create=True)
        else:
            bb = self.get_item(bb_id, self.BuildingBlock, create=True)
        return bb

    def write_mtb(self, f_path):
        temp_f = GromosFile(f_path, write=True)
        self.write_gromos_format(
            temp_f, 'TITLE', 'FORCEFIELD', 'MAKETOPVERSION', 'PHYSICALCONSTANTS', 'LINKEXCLUSIONS')
        for bb in self.bb:
            self.bb[bb].write_bb(temp_f)
#            temp_f.write_block('MTBUILDBLSOLUTE', self.bb[bb].write_bb)
        temp_f.f.close()


class Residue(GeneralContainer):
    #__slots__ = ('atoms', 'name')
    container2write = 'residues'

    def __init__(self, res_id=None, res_name='res', **kwargs):
        self.id = res_id
        self.name = res_name

    def add_atom(self, atom, **kwargs):
        self.add2container(atom, db_type=list, create = True, **kwargs)
        atom.res = self

    def find_atom(self, atom_name):
        for at in self.atoms:
            if at.name == atom_name:
                return at


class MolTop(GeneralTopology):
    Atom = TopAtom
    Residue = Residue

    def add_atom(self, atom_id=None, atom_name = None, **kwargs):
        at = self.get_item(atom_id, self.Atom, create=True, create_container=True, **kwargs)
        at.name = atom_name
        return at

    def add_residue(self, res_id=None, res_name=None, **kwargs):
        res = self.get_item(res_id, self.Residue, create=True, create_container=True, **kwargs)
        res.name = res_name
        return res

    def add_interactions(self, int_type, **kwargs):
        return self.get_item(item_id=None, klass = int_type, create=True, db_type = list, **kwargs)

    def generate_constraints(self, constraints='bonds', excl_bonds=None, flag_check_form=True, **kwargs):
        in_kwargs = dict(kwargs)
        if excl_bonds is None:
            excl_bonds = []
        for bond in self.get_container(constraints):
            flag = True
            for excl_bond in excl_bonds:
                if hasattr(excl_bond, 'atoms'):
                    excl_bond = excl_bond.atoms
                if excl_bond[0] in bond.atoms and excl_bond[1] in bond.atoms:
                    flag = False
                    break
            if flag:
                temp_const = Interaction(ConstraintType, **in_kwargs)
                temp_const.add_atom(*bond.atoms)
                for state in bond.states:
                    if state and state != Dummy:
                        bond_k = state.get_gr_gm_params()[0]
                        if flag_check_form:
                            temp_form = state.get_form()
                            if temp_form == 'gr':
                                fnc_type = 'gr_fnc'
                            elif temp_form == 'gm':
                                fnc_type = '1'
                        temp_const.add_state(params = bond_k, fnc_type=fnc_type, int_code=state.id)
                    else:
                        temp_const.states.append(state)
                self.add2container(temp_const, db_type=list, create=True)

    def get_HH(self):
        """get hydrogen and heavy atoms"""
        HH = ([], [])
        for at in self.get_atoms():
            if at.m > 1.5:
                HH[1].append(at)
            else:
                HH[0].append(at)
        return HH


class MoleculeType(MolTop, gmFragmentMoleculeIO, MolType_g2g):
    container2write = 'molecule_types'
    id_pref = 'mol_'

    def __str__(self):
        txt = str(self.id)
        if self.mol_name:
            txt += ' ' + str(self.mol_name)
        return txt

    def __repr__(self):
        return self.__str__()

    def __init__(self, molecule_id = None, nrexcl = 0, ff=None, **kwargs):
        self.atoms = OrderedDict()
        self.id = molecule_id
        self.mol_name = kwargs.get('molecule_name', molecule_id)
        self.nrexcl = nrexcl
        if ff is None:
            self.ff = FF()
        else:
            self.ff = ff

    def _gm_generate_excl_pairs(self, **kwargs):
        e_l, p_l = self._generate_atom_excl_pair_list(**kwargs)
        if kwargs.get('reset', True):
            db = self.get_container(PairType, flag_class=True, create=True, db_type=list)
            try:
                db.clear()
            except:
                del db[:]
            db = self.get_container(ExclusionType, flag_class=True, create=True, db_type=list)
            try:
                db.clear()
            except:
                del db[:]
        for at in self.get_atoms():
            for at2 in p_l[at]:
                at_pair = (at, at2)
                temp_state = self.ff.find_interaction_type(at_pair, PairType, flag_fnc_type=False, **kwargs)
                assert len(temp_state)==1, 'found more than one vdw pair type given at_pair'
                temp_interaction = Interaction(PairType, atoms=at_pair)
                temp_interaction.add_state(temp_state[0])
                self.add2container(temp_interaction, create=True, db_type=list)
                if kwargs.get('flag_excl_from_p_l', True):
                    temp_interaction = Interaction(ExclusionType, atoms=at_pair)
                    temp_interaction.add_state(fnc_type = None, params = (True,))
                    self.add2container(temp_interaction, create=True, db_type=list)
            for at2 in e_l[at]:
                at_pair = (at, at2)
                temp_interaction = Interaction(ExclusionType, atoms=at_pair)
                temp_interaction.add_state(fnc_type = None, params = (True,))
                self.add2container(temp_interaction, create=True, db_type=list)


class Molecule:
    container2write = 'molecules'

    def __str__(self):
        return str(self.mol_type) + ' ' + str(self.n)

    def __repr__(self):
        return self.__str__()


    def __init__(self, mol_type=None, num_mols=None, **kwargs):
        self.mol_type = mol_type
        self.n = num_mols


class Topology(Top_g2g, TopBBdb, MolTop, FF, topBlocksParser, topBlocksWriter, PTP_EDS_BlocksParser, gmTopologyIO):
    """Topology class - all atoms of the system + interactions + FF
    contains also Molecule class, which is very similar"""

    MoleculeType = MoleculeType
    Molecule = Molecule

    def __init__(self, parse_from = None, parse_from_file=True, format_type = 'gr', **kwargs):
        self.atoms = OrderedDict()
        if parse_from is not None:
            self.parse_top(parse_from, parse_from_file=parse_from_file, format_type = format_type, **kwargs)

    def parse_top(self, parse_from, parse_from_file = True, format_type = 'gr', **kwargs):
        """
        :param parse_from:
        :param parse_from_file:
        :param format_type:
        :param kwargs:
        :return:
        """
        fnc2call = {'gr':self.parse_top_gr, 'gm':self.parse_top_gm}
        assert format_type in fnc2call
        fnc2call[format_type](parse_from, parse_from_file, **kwargs)

    def write_top(self, f_path = None, format_type = 'gr', **kwargs):
        fnc2call = {'gr':self.write_top_gr, 'gm':self.write_top_gm}
        assert format_type in fnc2call
        fnc2call[format_type](f_path = f_path, **kwargs)

    '''
    # still missing - waiting to be implemented...
    def parse_ptp(self, parse_from, parse_from_file = True, **kwargs):
        """
        parse ptp state
        :param parse_from:
        :param parse_from_file:
        :param kwargs:
        :return:
        """
        self._parse_gr(parse_from, parse_from_file, **kwargs)
    '''

    def parse_eds(self, parse_from, parse_from_file = True, **kwargs):
        """
        parse eds states
        :param parse_from:
        :param parse_from_file:
        :param kwargs:
        :return:
        """
        self._parse_gr(parse_from, parse_from_file, **kwargs)

    def add_molecule(self, mol_type=None, num_mols=None, create = True, create_container = True, **kwargs):
        kwargs['create'] = kwargs.get('create', True)
        kwargs['create_container'] = kwargs.get('create_container', True)
        kwargs['db_type'] = kwargs.get('db_type', list)
        temp_mol = self.get_item(None, self.Molecule, mol_type=mol_type, num_mols=num_mols, **kwargs)
        return temp_mol

    def add_molecule_type(self, molecule_id = None, nrexcl = 0, ff=None, **kwargs):
        kwargs['create'] = kwargs.get('create', True)
        kwargs['create_container'] = kwargs.get('create_container', True)
        if ff is None:
            ff = self
        mol_type = self.get_item(molecule_id, self.MoleculeType, nrexcl = nrexcl, ff=ff, **kwargs)
        return mol_type

    def get_reduced_state(self, top_state):
        self.set_top_state(top_state)
        ff_dum = self.get_DUM_type
        atoms2keep = [at.id for at in self.get_atoms() if at.a_type != ff_dum]
        return self.reduce_top(atoms2keep)

    def reduce_top(self, atom_ids = None, **kwargs):
        """
        :param atom_ids:
        :param kwargs:
            flag_inplace
        :return:
        """
        txt = '\treduced topology: ' + str(atom_ids)
        top_red = self.reduce(atom_ids, **kwargs)
        sys_title = top_red.get_container(DescriptionPart, flag_class=1, allow_not_found = 1)
        if sys_title:
            sys_title[0].lines.insert(0,txt)
        else:
            dp = DescriptionPart()
            dp.add_line(txt)
            top_red.add2container(dp, list_index=0, create=True)
        return top_red

    def get_rev_top_ptp(self):########################################################
        new_top = self.reduce_top(list(self.atoms.keys()))
        new_top.add_perturbation()
        for at in self.ptp.atoms:
            at_id = at.id
            new_at = new_top.atoms[at_id]
#            rev_name = getattr(at, 'rev_name', False)
#            if rev_name:
#                new_at.name = rev_name
#                new_at.res.name = at.rev_res_name
            tp0, tp1, tp2 = new_top.ff.a_type[new_at.a_type.id], new_at.m, new_at.p_ch
            new_at.a_type = new_top.ff.a_type[self.ptp.atoms[at][0].id]
            new_at.m = float(self.ptp.atoms[at][1])
            new_at.p_ch = float(self.ptp.atoms[at][2])
            new_top.ptp.atoms[new_at] = (tp0, tp1, tp2, self.ptp.atoms[at][3])
        for pair_id in self.ptp.pair_vdw:
            new_pair_id = (new_top.atoms[pair_id[0].id], new_top.atoms[pair_id[1].id])
            new_top.ptp.pair_vdw[new_pair_id] = tuple(reversed(self.ptp.pair_vdw[pair_id]))
            new_top.ptp.pair_vdw[new_pair_id] = (self.ptp.pair_vdw[pair_id][1], self.ptp.pair_vdw[pair_id][0])
            if self.ptp.pair_vdw[pair_id][0] == 0:
                if new_pair_id[0] in new_pair_id[1].e_l:
                    new_pair_id[1].e_l.pop(new_pair_id[1].e_l.index(new_pair_id[0]))
                if new_pair_id[1] in new_pair_id[0].e_l:
                    new_pair_id[0].e_l.pop(new_pair_id[0].e_l.index(new_pair_id[1]))
            if self.ptp.pair_vdw[pair_id][0] == 2:
                if new_pair_id[0] in new_pair_id[1].p_l:
                    new_pair_id[1].p_l.pop(new_pair_id[1].p_l.index(new_pair_id[0]))
                if new_pair_id[1] in new_pair_id[0].p_l:
                    new_pair_id[0].p_l.pop(new_pair_id[0].p_l.index(new_pair_id[1]))
            if self.ptp.pair_vdw[pair_id][1] == 0:
                new_pair_id[0].e_l.append(new_pair_id[1])
            if self.ptp.pair_vdw[pair_id][1] == 2:
                new_pair_id[0].p_l.append(new_pair_id[1])
        for bond_type in ('bonds', 'angles', 'impropers', 'dihedrals'):
            if hasattr(self.ptp, bond_type):
                temp_ptp_ints = getattr(self.ptp, bond_type)
                temp_ints = getattr(self, bond_type)
                new_temp_ptp_ints = getattr(new_top.ptp, bond_type)
                new_temp_ints = getattr(new_top, bond_type)
                for temp_int in temp_ptp_ints:
                    temp_int_i = temp_ints.index(temp_int)
                    new_temp_int = new_temp_ints[temp_int_i]
                    new_temp_ptp_ints[new_temp_int] = tuple(reversed(temp_ptp_ints[temp_int]))
                    new_temp_int.p = getattr(new_top.ff, bond_type)[temp_ptp_ints[temp_int][1]]
        new_top._gr_exp_red_sort_excl_pair()
        return new_top

    def get_molecules(self, at=None, flag_sort=False, flag_mol_atom_sort = False, **kwargs):
        """
        :param at: search in a subset of atoms
        :param flag_sort: sort atoms
        :param flag_mol_atom_sort: sort molecules
        :param kwargs:
            flag_disres - use distance restraints in addition to bonds
        :return:
        """
        if flag_mol_atom_sort:flag_sort = True
        self.create_adj_list()
        if at:
            return tuple(self.DFS(at))
        mols = []
        atoms = list(self.get_atoms())
        if kwargs.get('flag_disres'):
            G = self.sub_graph(self.atoms.values())
            disres_cont = self.get_container(self.Distance_r, flag_class=True,  allow_not_found = True)
            if disres_cont:
                for dr in disres_cont:
                    G.add_edge(dr.atoms)
        else:
            G = self
        while atoms:
            at = atoms[0]
            mols.append([])
            for at in G.DFS(at):
                mols[-1].append(at)
                atoms.remove(at)
            if flag_sort:
                new_mol = []
                for at in self.get_atoms():
                    if at in mols[-1]:
                        new_mol.append(at)
                assert len(new_mol)==len(mols[-1])
                mols[-1] = new_mol
            mols[-1] = tuple(mols[-1])
        if flag_mol_atom_sort:
            kwargs['flag_mol_atom_sort'] = True
            _, _, mols = self.check_molecule_atoms_order(molecules = mols, **kwargs)
        return tuple(mols)

    def check_molecule_atoms_order(self, molecules = None, print_warn = True, **kwargs):
        if molecules is None:
            molecules = self.get_molecules(flag_sort = 1)
        top_atoms = list(self.atoms.values())
        last_at_index = [top_atoms.index(mol[-1]) for mol in molecules]
        flag_check = last_at_index==sorted(last_at_index)
        if print_warn and not flag_check:
            do_warn('\n\tmolecule atoms are not sorted!!!\nlast atom index of each molecule:\n'+str(last_at_index))
        if kwargs.get('flag_mol_atom_sort'):
            molecules = []
            temp_max = -1
            for i, lai in enumerate(last_at_index):
                if temp_max > lai:
                    continue
                molecules.append(tuple(top_atoms[temp_max+1:lai+1]))
                temp_max = lai
        return flag_check, last_at_index, tuple(molecules)

    def get_G_mol(self, m=None, **kwargs):
        mols = self.get_molecules(**kwargs)
        if m is None:
            mG = []
            for tm in mols:
                mG.append(self.sub_graph(tm))
            return mG
        else:
            return self.sub_graph(mols[m])

    """
    def improper_type(self, imp):
        bond_c = []
        for at in imp.atoms:
            c = 0
            for bond_at in self.adj[at]:
                if bond_at in imp.atoms:
                    c += 1
            bond_c.append(c)
        return bond_c
    """

    def _check_interaction_atoms(self, bi, atom_ids):
        flag = True
        for at in bi.atoms:
            if at.id not in atom_ids:
                flag = False
                break
        return flag

    def renumber(self, **kwargs):
        """
        :param kwargs:
            c_init = 1
            attrib = 'id'
            flag_id_map = False
        :return:
            id_map
        """
        id_map_atoms = self.renumber_container('atoms', **kwargs)
        id_map_residues = self.renumber_container('residues', **kwargs)
        return id_map_atoms, id_map_residues


class ConfAtom():
    container2write = 'atoms'
    _attributes = ('id', 'name', 'res_id', 'res_name', 'coord', 'vel', 'latshift')


class Box():
    def __init__(self, a = None, b = None, c = None, vec = None, **kwargs):
        self.euler_angles = np.zeros(3)
        self.origin =  np.zeros(3)
        if a is None:
            if vec is None:
                self.gr_box_type = 0
                self.abc = np.zeros(3)
                self.angles = np.zeros(3)
                self.euler_angles = np.zeros(3)
                self.origin =  np.zeros(3)
                self.vec = np.zeros((3,3))
            else:
                self.vec = np.zeros((3,3))
                if len(vec) == 3:
                    for i in range(3):
                        self.vec[i,i] = vec[i]
                elif len(vec) == 9:
                    self.vec = np.array(vec).reshape(3,3)
                else:
                    self.vec = vec
                self._gm_vec = vec
                    ######### there are more cases... make sure to 
        else:
            self.gr_box_type = int(kwargs.get('gr_box_type', 1))
            self.euler_angles = kwargs.get('euler_angles', np.zeros(3))
            self.origin = kwargs.get('origin', np.zeros(3))
            self.angles = np.array(kwargs.get('angles', np.ones(3) * 90))
            assert self.euler_angles.shape == (3,)
            assert self.origin.shape == (3,)
            if b is None:
                self.abc = np.array([a,a,a])
                self.abc2vec()
            else:
                try:
                    if len(a) == 3:
                        assert c is not None
                        self.vec = np.array([a,b,c])
                        assert self.vec.shape == (3,3)
                        self.vec2abc()
                except TypeError:
                    assert c is not None
                    self.abc = np.array([a,b,c])
                    assert self.abc.shape == (3,)
                    self.abc2vec()

    def vec2abc(self):
        self.abc = np.linalg.norm(self.vec, axis = 1)
        self.angles = []
        for p in ((1,2), (0,2), (0,1)):
            temp_cos = self.vec[p[0]].dot(self.vec[p[1]]) / np.prod(self.abc[p,])
            temp_ang = np.degrees(np.arccos(temp_cos))
            self.angles.append(temp_ang)
        self.angles = np.array(self.angles)

    def abc2vec(self):
        v1 = np.zeros(3)
        v1[0] = self.abc[0]
        v2 = rot(v1, self.angles[2]) * self.abc[1] / self.abc[0]
        x3 = np.cos(np.radians(self.angles[1])) * self.abc[2]
        y3 = (np.cos(np.radians(self.angles[0])) * np.prod(self.abc[:2]) - v2[0] * x3) / v2[1]
        z3 = np.sqrt(self.abc[2]**2 - x3**2 - y3**2)
        self.vec = np.array([v1, v2, [x3, y3, z3]])


class Configuration(cnfBlocksParser, cnfBlocksWriter, gmConfigurationIO):
    """
    Configuration class
        main container: atoms
        _N_dim (usually 3)
        _coord - np.array to store coordinate
        _vel - np.array to store velocities
    """
    Atom = ConfAtom
    Box = Box

    def __init__(self, f_path = None, N_dim=3, **kwargs):
        self._N_dim = N_dim
        if f_path and f_path.endswith('cnf'):
            self.parse_cnf(f_path, **kwargs)
            self._generate_coord()
        if f_path and f_path.endswith('gro'):
            self.parse_gro(f_path)

    def _generate_coord(self):
        """
        generate self._coord and link to each atom object
        """
        self._coord = np.empty((len(self.atoms), self._N_dim))
        for i, at in enumerate(self.atoms):
            self._coord[i] = at.coord
            at.coord = self._coord[i]

    def _d2(self, c1, c2):
        temp_d2 = c1 - c2
        return np.dot(temp_d2, temp_d2)

    def d2(self, at1, at2, periodic=True):
        if periodic:
            c1 = self.atoms[at1].coord % self.box.abc
            c2 = self.atoms[at2].coord % self.box.abc
            d2_min = list(self.box.abc)
            for coord_i in range(len(d2_min)):
                for i in range(-1, 2):
                    temp = abs(c1[coord_i] + self.box.abc[coord_i] - c2[coord_i])
                    if d2_min[coord_i] > temp:
                        d2_min[coord_i] = temp
            d2_min = np.dot(d2_min, d2_min)
            return d2_min
        else:
            return self._d2(self.atoms[at1].coord, self.atoms[at2].coord)

    def d2_within_c2(self, at1, at2, cutoff2, periodic=True):
        temp_d2 = self.d2(at1, at2, periodic)
        if temp_d2 < cutoff2:
            return True
        else:
            return False

    def d2_within(self, at1, at2, cutoff, periodic=True):
        self.d2_within_c2(at1, at2, cutoff**2, periodic)

    def find_within_c2(self, atom, atoms_others, cutoff2, periodic=True):
        at_within = []
        for at in atoms_others:
            if self.d2_within_c2(atom, at, cutoff2, periodic):
                at_within.append(at)
        return at_within

    def find_within(self, atom, atoms_others, cutoff,  periodic=True):
        return self.find_within_c2(atom, atoms_others, cutoff**2, periodic)

    def find_within_box(self, atoms, cutoff, periodic=True):
        if len(atoms) == 1:
            return self.find_within(atoms[0], list(range(len(self.atoms))), cutoff, periodic)
        cutoff2 = cutoff**2
        n_at = len(atoms)
        m_d2 = np.empty((n_at, n_at))
        for i in range(n_at - 1):
            for j in range(i, n_at):
                m_d2[i, j] = self.d2(atoms[i], atoms[j])
                m_d2[j, i] = m_d2[i, j]
        cent_at_i = 0
        mmd2 = max(m_d2[0])
        for i in range(1, n_at):
            temp_m = max(m_d2[i])
            if temp_m < mmd2:
                mmd2 = temp_m
                cent_at_i = i
        cent_at_i = atoms[i]
        cutoff2_2 = (cutoff + math.sqrt(mmd2) + 0.00001)**2
        list_cutoff22 = set(self.find_within_c2(cent_at_i, list(range(len(self.atoms))), cutoff2_2))
        at_within = []
        for at in atoms:
            temp_at_within = self.find_within_c2(at, list_cutoff22, cutoff2)
            at_within.extend(temp_at_within)
            list_cutoff22 -= set(temp_at_within)
        return at_within

    def cog(self, atoms):
        cog = np.array(self.atoms[atoms[0]].coord)
        n_at = len(atoms)
        for i in range(1, n_at):
            cog += self.atoms[atoms[i]].coord
        return cog / n_at

    def parse_cnf(self, f_path):
        self._parse_gr(f_path)

    def write_cnf(self, f_path):
        temp_f = GromosFile(f_path, write=True)
        self.write_gromos_format(temp_f, 'TITLE')
        int_blocks2write = ('TIMESTEP', 'POSITION', 'VELOCITY', 'GENBOX')
        self.write_gromos_format(temp_f, *int_blocks2write)

    def write_rpr(self, f_path):
        temp_f = GromosFile(f_path, write=True)
        temp_title = list(self.title)
        temp_title.insert(0, 'reference positions for position restraining based on')
        temp_f.write_block('TITLE', temp_title)
        self.write_gromos_format(temp_f, 'REFPOSITION')

    def write_por(self, f_path, *args, **kwargs):
        temp_f = GromosFile(f_path, write=True)
        temp_title = list(self.title)
        temp_title.insert(0, str(kwargs))
        temp_title.insert(0, 'position restraints specifications')
        temp_f.write_block('TITLE', temp_title)
        self.write_gromos_format(temp_f, 'POSRESSPEC', *args, **kwargs)

    def renumber(self):
        for i in range(len(self.atoms)):
            self.atoms[i].id = str(i + 1)

    def reduce_cnf(self, atoms):
        new_cnf = Configuration()
        new_cnf.title = ['reduced_cnf']
        new_cnf.box = self.box
        new_cnf.atoms = []
        for at_i in range(len(self.atoms)):
            if at_i in atoms:
                old_at = self.atoms[at_i]
                new_at = new_cnf.Atom()
                for attr in new_at._attributes:
                    if hasattr(old_at, attr):
                        setattr(new_at, attr, getattr(old_at, attr))
                new_cnf.add2container(new_at, db_type=list)
                new_at.id = str(len(new_cnf.atoms))
        return new_cnf

class Trajectory(TrjCnfBlocksParser, TrjCnfBlocksWriter):
    def __init__(self, f_path = None, N_atoms = None, int_num_dtype = np.int32, real_num_dtype = np.float32, **kwargs):
        self.int_num_dtype = int_num_dtype
        self.real_num_dtype = real_num_dtype
        self.time_real_num_dtype = kwargs.get('time_dtype', np.float32)
        self.__adjust_real_num_format()
        self.N_atoms = N_atoms
        self._trj = []
        self.trj = None
        if N_atoms:
            self.get_frame_dtype()
        if f_path:
            self.parse_trc(f_path, **kwargs)

    def __adjust_num_format(self, dtype, atr):
        it_size = dtype().itemsize
        if it_size<8:
            num_dig = it_size * 2 - 1
            num_spaces = str(15 - 9 + num_dig)
            setattr(self, atr, '{:' + num_spaces + '.' + str(num_dig) + 'f}' + ' ' * (9 - num_dig))

    def __adjust_real_num_format(self):
        self.__adjust_num_format(self.real_num_dtype, '_real_num_format')
        self.__adjust_num_format(self.time_real_num_dtype, '_time_real_num_format')

    def _get_frame_dtype_original(self):
        frame_dtype = [('time_step', self.int_num_dtype, (1,)), ('time', np.float32, (1,))]
        frame_dtype.append(('coord', self.real_num_dtype, (self.N_atoms, 3)))
        frame_dtype.extend([('box_type', np.int8, (1,)), ('box', self.real_num_dtype, (4,3))])
        self.frame_dtype = np.dtype(frame_dtype)
        self.trj = np.empty(0, dtype=self.frame_dtype)

    def _get_frame_dtype_simplified(self):
        frame_dtype = [('time_step', self.int_num_dtype), ('time', np.float32)]
        frame_dtype.append(('coord', self.real_num_dtype, (self.N_atoms, 3)))
        frame_dtype.extend([('box_type', np.int8), ('box', self.real_num_dtype, (4,3))])
        self.frame_dtype = np.dtype(frame_dtype)
        self.trj = np.empty(0, dtype=self.frame_dtype)
    
    def get_frame_dtype(self, fnc2call=None):
        if fnc2call:
            fnc2call = getattr(self,fnc2call)
            fnc2call()
        else:
            self._get_frame_dtype_original()

    def set_N_frames(self, N_frames):
        self.N_frames = N_frames
        self.trj.resize(N_frames)

    def reduce_N_atoms(self, N_atoms=None, atom_sel=None):
        """
        reduce number of atoms
        :param f_path: N_atoms
        :param atom_sel: selection of atoms e.g. (1,2,5,8) - will be used to index atoms in the trj array
        :return: new Trajectory instance
        """
        if atom_sel:
            N_atoms=len(atom_sel)
            atom_sel = tuple(atom_sel)
        assert N_atoms is not None
        new_trj = Trajectory(N_atoms=N_atoms, int_num_dtype=self.int_num_dtype, real_num_dtype=self.real_num_dtype,
                            time_dtype=self.time_real_num_dtype)
        new_trj.set_N_frames(self.trj.shape[0])
        for var_name in self.frame_dtype.names:
            if var_name=='coord':
                if atom_sel is None:
                    new_trj.trj[var_name][:] = self.trj[var_name][:,:N_atoms,:]
                else:
                    new_trj.trj[var_name][:] = self.trj[var_name][:,atom_sel,:]
            else:
                new_trj.trj[var_name][:] = self.trj[var_name][:]
        return new_trj.trj

    def _gr_add_frame(self):
        temp_frame = np.array((self.time_step, self.time, self.coord, self.box_type, self.box), dtype=self.frame_dtype)
        self._trj.append(temp_frame)

    def load_trc_npz(self, f_path):
        temp_title, self.trj = np.load(f_path).values()
        self.add_gr_title([str(temp_title)])
        self.N_atoms = self.trj['coord'].shape[1]
        self.frame_dtype = self.trj.dtype
        self.int_num_dtype = self.frame_dtype[0].subdtype[0].type
        self.time_real_num_dtype = self.frame_dtype[1].subdtype[0].type
        self.real_num_dtype = self.frame_dtype[2].subdtype[0].type
        self.__adjust_real_num_format()

    def parse_trc(self, f_path, **kwargs):
        if f_path.endswith('.gz'):
            self._parse_gr(f_path, fnc2open = gzip.open)
        else:
            self._parse_gr(f_path, **kwargs)
        if not self._trj:
            self.N_atoms = 0
            self.get_frame_dtype()
        self._trj = np.array(self._trj, dtype=self.frame_dtype)
        self.trj = self._trj


class MolSystem(Topology, Configuration):
    def __init__(self):
        self.atoms = OrderedDict()

from SMArt.md.ana.incl import get_lset

class MD_Parameters(GromosParser, GromosWriter):
    def __init__(self, template_md_in):
        self.md_in = template_md_in
        if self.md_in.endswith('imd'):
            self._parse_gr(self.md_in)
    
    def change_md_in(self, **kwargs):
        if self.md_in.endswith('mdp'):
            return self.change_mdp(**kwargs)

    def change_mdp(self, md_kw, **kwargs):
        s = ''
        mdp_in = self.md_in
        f = open(mdp_in)
        for l in f:
            temp = l.strip().lower().split()
            if temp:
                temp = temp[0]
            else:continue
            if temp in md_kw:
                s += '{:<25}= {:}\n'.format(temp, md_kw.pop(temp))
            else:
                s += l
        for temp in md_kw:
            s += '{:<25}= {:}\n'.format(temp, md_kw[temp])
        mdp_out = kwargs.get('mdp_out')
        if mdp_out:
            f=open(mdp_out, 'w')
            f.write(s)
            f.close()
        else:
            return s

    def __get_bl_params(self, bl):
        params = []
        params_per_line = []
        for l in self.undefined_bl[bl]:
            temp = l.split()
            params_per_line.append(len(temp))
            params.extend(temp)
        return params, params_per_line

    def __generate_bl_txt(self, bl, params, params_per_line):
        bl_lines = []
        bl_txt = ''
        c_param = 0
        for N_params in params_per_line:
            line = ''
            for i in range(N_params):
                line += '{:>9} '.format(params[i + c_param])
            bl_lines.append(line[:-1] + '\n')
            c_param += N_params
        self.undefined_bl[bl] = bl_lines

    def change_imd(self, md_kw, **kwargs):
        for bl in md_kw:
            bl_params, bl_params_per_line = self.__get_bl_params(bl)
            for param, value in md_kw[bl].items():
                bl_params[param] = value
            self.__generate_bl_txt(bl, bl_params, bl_params_per_line)

    def get_LPs_pred(self):
        if self.md_in.endswith('mdp'):
            mdp_in = self.md_in
            f = open(mdp_in)
            for l in f:
                if l.startswith('fep-lambdas'):
                    temp = l.split('=')[1]
                    LPs_pred = get_lset(temp.split())
            return LPs_pred
    
    def write_imd(self, f_path = None, **kwargs):
        gs = self._get_grs_from_fpaht(f_path)
        self.write_gromos_format(gs, 'TITLE')
        blocks2write = list(self.undefined_bl)
        self.write_gromos_format(gs, *blocks2write, **kwargs)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()



from SMArt.incl import copy, Defaults, do_warn, np, OrderedDict, permutations

_atomtype_name__element_map = {'O': '8', 'OM': '8', 'OA': '8', 'OE': '8', 'OW': '8', 'N': '7', 'NT': '7', 'NL': '7',
    'NR': '7', 'NZ': '7', 'NE': '7', 'C': '6', 'CH0': '6', 'CH1': '6', 'CH2': '6', 'CH3': '6', 'CH4': '6', 'CH2r': '6',
    'CR1': '6', 'HC': '1', 'H': '1', 'DUM': '0', 'S': '16', 'CU1+': '29', 'CU2+': '29', 'FE': '26', 'ZN2+': '30',
    'MG2+': '12', 'CA2+': '20', 'P': '15', 'SI': '14', 'AR': '18', 'F': '9', 'CL': '17', 'BR': '35', 'CMet': '6',
    'OMet': '8', 'NA+': '11', 'CL-': '17', 'CChl': '6', 'CLChl': '17', 'HChl': '1', 'SDmso': '16', 'CDmso': '6',
    'ODmso': '8', 'CCl4': '6', 'CLCl4': '17', 'FTFE': '9', 'CTFE': '6', 'CHTFE': '6', 'OTFE': '8', 'CUrea': '6',
    'OUrea': '8', 'NUrea': '7', 'CH3p': '6'}

_element__AtomicNum_map = dict(H=1, C=6, N=7, O=8, F=9, NA=11, MG=12, AL=13, SI=14, P=15, S=16, CL=17, AR=18, \
                                K=19, CA=20, FE=26, CU=29, ZN=30, BR=35)
_AtomicNum__element_map = dict(item[::-1] for item in _element__AtomicNum_map.items())


def gr_get_atom_element(a_type):
    if a_type.name == 'P,SI':
        if round(a.m_type.m) == 31:
            return _atomtype_name__element_map['P']
        else:
            return _atomtype_name__element_map['SI']
    else:
        return _atomtype_name__element_map[a_type.name]

class ContPref_g2g(Defaults):
    ### GROMOS thingies
    def _set_gr_interaction_type_container(self, int_type):
        int_type.cont_pref = self._gr_cont_pref

    ### GROMACS thingies
    def _set_gm_interaction_type_container(self, int_type):
        int_type.cont_pref = self._gm_cont_pref

_ContPref_g2g_defs = {}
_ContPref_g2g_defs['_gr_cont_pref'] = 'gr_'
_ContPref_g2g_defs['_gm_cont_pref'] = 'gm_'
ContPref_g2g._add_defaults(_ContPref_g2g_defs, flag_set=1)


class IntType_g2g(ContPref_g2g):
    @classmethod
    def get_int_type(cls, **kwargs):
        return cls

    def get_gr_gm_params(self):
        gr_gm_params = []
        for i in self._gr_gm_pos_fnc_map[self.fnc_type]:
            if i is None:
                gr_gm_params.append(None)
            else:
                gr_gm_params.append(self.p[i])
        return gr_gm_params

    def generate_params_gr_gm(self, gr_gm_params):
        params = [0] * len(self.flag) # PTP flag
        for i1, i2 in enumerate(self._gr_gm_pos_fnc_map[self.fnc_type]):
            if i2 is not None and gr_gm_params[i1] is not None:
                params[i2] = gr_gm_params[i1]
        self.add_params(params)

# gr2gm
    @classmethod
    def _get_gr2gm_int_type_class(cls, dv = None, **kwargs):
        if dv is None:
            dv = cls._get_current_class_dvalues(**kwargs)
        int_type = cls
        int_type_def_key = '_' + cls.__name__ + '__gr2gm_int_type'
        if int_type_def_key in dv:
            int_type = dv[int_type_def_key]
        return int_type

    @classmethod
    def get_gr2gm_int_type(cls, **kwargs):
        dv = cls._get_current_class_dvalues(**kwargs)
        int_type = cls._get_gr2gm_int_type_class(dv = dv)
        kwargs['fnc_type'] = kwargs.get('fnc_type', dv['_' + cls.__name__ + '__gr2gm_fnc'])
        return int_type(**kwargs)

    @classmethod
    def get_gr2gm_def_pref_class(cls, **kwargs):
        """
        gets gm_def_pref for a given class, e.g. gb_ for BondType
        :param kwargs:
            flag_not_found
        :return:
        """
        dv = kwargs.get('dv', cls._get_current_class_dvalues(**kwargs))
        return dv.get('_' + cls.__name__ + '__gr2gm_def_pref', kwargs.get('flag_not_found' ,cls.__name__ + '_'))

    def get_gr2gm_def_pref(self, **kwargs):
        """
        gets gm_def_pref for a given class, e.g. gb_ for BondType
        :param kwargs:
            flag_not_found
        :return:
        """
        cls = self.__class__
        dv = cls._get_current_class_dvalues(**kwargs)
        temp_k = '_' + cls.__name__ + '__gr2gm_def_pref_fnc_type'
        if temp_k in dv and self.fnc_type in dv[temp_k]:
            return dv[temp_k][self.fnc_type]
        kwargs['dv'] = dv
        return cls.get_gr2gm_def_pref_class(**kwargs)

    def gr2gm(self, **kwargs):
        kwargs['int_code'] = kwargs.get('int_code', self.id)
        kwargs['atoms'] = kwargs.get('atoms', self.atoms)
        if kwargs.get('flag_params') or not hasattr(self, '_gr_gm_pos_fnc_map'):
            gm_int_type = self.get_gr2gm_int_type(params = self.p, **kwargs)
            return gm_int_type
        gm_int_type = self.get_gr2gm_int_type(**kwargs)
        gr_gm_params = self.get_gr_gm_params()
        if hasattr(self, 'gm2gr_param_trans_fnc'):
            for i, temp_fnc in enumerate(self.gr2gm_param_trans_fnc):
                if temp_fnc:
                    gr_gm_params[i] = temp_fnc(gr_gm_params[i])
        gm_int_type.generate_params_gr_gm(gr_gm_params)
        return gm_int_type

# gm2gr
    def _get_gm2gr_int_type_class(self, **kwargs):
        cls = self.__class__
        int_type = cls
        dv = cls._get_current_class_dvalues(**kwargs)
        temp_k = '_' + cls.__name__ + '__gm2gr_int_type_fnc_type'
        if temp_k in dv and self.fnc_type in dv[temp_k]:
            int_type = dv[temp_k][self.fnc_type]
        return int_type

    def get_gm2gr_int_type(self, **kwargs):
        int_type = self._get_gm2gr_int_type_class(**kwargs)
        kwargs['fnc_type'] = 'gr_fnc'
        return int_type(**kwargs)

    def gm2gr(self, **kwargs):
        kwargs['int_code'] = kwargs.get('int_code', self.id)
        kwargs['atoms'] = kwargs.get('atoms', self.atoms)
        gr_int_type = self.get_gm2gr_int_type(**kwargs)
        gr_gm_params = self.get_gr_gm_params()
        if hasattr(self, 'gm2gr_param_trans_fnc'):
            for i, temp_fnc in enumerate(self.gm2gr_param_trans_fnc):
                if temp_fnc:
                    if gr_gm_params[i] is None:
                        gr_gm_params[i] = None
                    else:
                        gr_gm_params[i] = temp_fnc(gr_gm_params[i])
        if kwargs.get('verbose'):
            print(gr_gm_params)
        gr_int_type.generate_params_gr_gm(gr_gm_params)
        return gr_int_type


class FFg2g(ContPref_g2g):

    def __just_for_readability(self):
        self.get_intDB = None  # from FF
        self._gr_block_names = None # from gromos FF
        self.gm_FF_Defaults = None # from gromacs FF
        self.add_define = None # from gromacs FF
        self.self._check_vdw_pair_gm_params = None # from gmFFWriter

    def __check_ff_defs(self):
        for d in getattr(self, 'defines',{}):
            if d.startswith('_FF_'):
                return d
        return

# gr2gm
    def gr2gm(self, **kwargs):
        self.ff_gr2gm(**kwargs)

    def ff_gr2gm(self, **kwargs):
        """
        :param kwargs:
            none
        :return:
        """
        if not self.__check_ff_defs():
            self.add_define(['_FF_GROMOS96'])
            ff_name = getattr(self, self._gr_block_names.FORCEFIELD)
            if ff_name:
                self.add_define(['_FF_GROMOS' + ff_name.upper()])
        if not hasattr(self, 'gm_ff_defaults'):
            self.gm_ff_defaults = self.gm_FF_Defaults(['1', '1', 'no', 1.0, 1.0])

        # add SI from P,SI

        for at_id in self.get_intDB().a_type:
            at = self.get_intDB().a_type[at_id]
            at.m = at.p_ch = 0.
            at.p_type = 'A'
            at.vdw = list(self.get_vdw(at.gr_id, at.gr_id)[0].p)
            anames = at.name.split(',')
            at.gm_id = anames[0]
            if not hasattr(at, 'element'):
                at.element = '0'
                flag = 1
                for at_name1 in _atomtype_name__element_map:
                    for at_name2 in anames:
                        if at_name1 == at_name2:
                            at.element = _atomtype_name__element_map[at_name1]
                            flag = 0
                            break
                    else:
                        continue
                    break
                if flag:
                    for at_name1 in _atomtype_name__element_map:
                        for at_name2 in anames:
                            if at_name1.upper() == at_name2.upper():
                                at.element = _atomtype_name__element_map[at_name1]
                                break
                        else:
                            continue
                        break
        ff = self.get_intDB()
        for vdw_pair in ff.vdw.values():
            for i in range(2):
                vdw_pair[i].gm_int_type = vdw_pair[i].gr2gm(**kwargs)
        self._check_vdw_pair_gm_params()
        for container_name in ff._containers:#############################################
            if container_name.startswith(self._gr_cont_pref):
                temp_cont = ff.get_container(container_name)
                for int_type in temp_cont.values():
                    gm_int_type = int_type.gr2gm(**kwargs)
                    temp_def = gm_int_type.gm_params2define()
                    self.add2container(temp_def, replace = kwargs.get('flag_def_replace'))
                    int_type.gm_int_type = gm_int_type

# gm2gr
    def __get_gm_def_pref(self):
        for int_type_name in self._available_int_types:
            int_type = getattr(self, int_type_name)
            if hasattr(int_type, 'get_gr2gm_def_pref_class'):
                gm_def_pref = int_type.get_gr2gm_def_pref_class(flag_not_found=False, flag_do_warn = 0)
                if gm_def_pref:
                    yield gm_def_pref, int_type

    def ff_gm2gr(self, **kwargs): ######################################################## could be better with reading gr ff
        if not hasattr(self.get_intDB(), 'vdw'):
            self.generate_vdw()
        for at_id in self.get_intDB().a_type:
            at = self.get_intDB().a_type[at_id]
            temp_vdw = self.get_intDB().vdw[(at,)*2][0].gm2gr()
            at.c612 = list(np.sqrt(temp_vdw.p))
            for i in range(2):
                at.c612.append(at.c612[-1]) ############################################## could be better with reading gr ff
            temp_pair = self.get_intDB().vdw[(at,)*2][1]
            if temp_pair:
                temp_pair = temp_pair.gm2gr()
                at.c612.extend(np.sqrt(temp_pair.p))
            else:
                at.c612.extend(at.c612[:2])
        self.get_intDB().renumber_container('a_type', attrib = 'gr_id')
        flag_ff_defs = self.__check_ff_defs()
        if flag_ff_defs and 'GROMOS' in flag_ff_defs.upper():
            gm_def_pref_int_type = dict(self.__get_gm_def_pref())
            for d in self.defines:
                for gm_def_pref in gm_def_pref_int_type:
                    if d.startswith(gm_def_pref):
                        temp_id = d[len(gm_def_pref):]
                        temp_gr_int_type = gm_def_pref_int_type[gm_def_pref](fnc_type = 'gr_fnc')
                        temp_gm_int_type = temp_gr_int_type.get_gr2gm_int_type(params = self.defines[d].v)
                        new_gr_int_type = temp_gm_int_type.gm2gr(int_code=temp_id)
                        new_gr_int_type.cont_pref = self._gr_cont_pref
                        self.get_intDB().add2container(new_gr_int_type, create = True)
            for gr_int_type in gm_def_pref_int_type.values():
                temp_kw = dict(flag_class = 1, allow_not_found = True)
                temp_gm_conts = [self.get_intDB().get_container(gr_int_type, cont_pref = self._gm_cont_pref, **temp_kw)]
                #temp_gm_conts.append(self.get_intDB().get_container(gr_int_type, **temp_kw))
                for temp_gm_cont in temp_gm_conts:
                    if temp_gm_cont:
                        for temp_gm_int_type in temp_gm_cont:
                            temp_gr_int_type = temp_gm_int_type.gm2gr()
                            temp_gr_int_type.cont_pref = self._gr_cont_pref
                            if not self.check_eq_int_type(temp_gr_int_type):
                                self.get_intDB().add2container(temp_gr_int_type, create = True)

    def gm2gr(self, **kwargs):
        self.ff_gm2gr(**kwargs)


class Top_g2g:
    """
    g2m (gr2gm & gm2gr) functionality for topology
    """
    def __just_for_readability(self):
        self.ff_gr2gm = None
        self.molecules = None
        self.molecule_types = None

# gr2gm
    def top_gr2gm(self, **kwargs):
        """

        :param kwargs:
            ff_gr2gm(**kwargs)
            generate_mol_types(*mol_ids, **kwargs)
            mol_ids - names for each molecule type
            flag_add_molecule - adding molecule 1 for each molecule type
            num_mols - adding more than 1 molecule for each molecule type
        :return:
        """
        self.ff_gr2gm(**kwargs)
        """
        if kwargs.get('reset_excl_pairs', True):
            self._gr_exp_red_sort_excl_pair()
            self.get_exclusions_from_e_l(reset = 1)
            self.get_pairs_from_p_l(reset = 1)
        """
        if kwargs.get('flag_clear_mols_types'):
            self.molecules.clear()
            self.molecule_types.clear()
        mol_ids = kwargs.get('mol_ids', [])
        self.generate_mol_types(*mol_ids, **kwargs)

    def gr2gm(self, **kwargs):
        self.top_gr2gm(**kwargs)

    def generate_mol_types(self, *mol_ids, **kwargs):
        """

        :param mol_ids:  names for each molecule type
        :param kwargs:
            flag_add_molecule - adding molecule 1 for each molecule type
            num_mols - adding more than 1 molecule for each molecule type

        :return:
        """
        if kwargs.get('flag_add_molecule', True):
            kwargs['num_mols'] = kwargs.get('num_mols', 1)
        mols = self.get_molecules(flag_mol_atom_sort = 1, **kwargs)
        for i,m in enumerate(mols):
            atom_ids = [at.id for at in m]
            if i<len(mol_ids):
                mol_id = mol_ids[i]
            else:
                mol_id = None
            mol_type = self.generate_mol_type(atom_ids, molecule_id=mol_id, **kwargs)
            if kwargs.get('flag_add_molecule', True):
                self.add_molecule(mol_type = mol_type, **kwargs)
        flag_sol = kwargs.get('flag_sol', '#include "gromos54a7.ff/spc.itp"')
        if flag_sol:
            self.parse_gm(flag_sol, parse_from_file = 0)
        flag_ion = kwargs.get('flag_ion', '#include "gromos54a7.ff/ions.itp"')
        if flag_ion:
            self.parse_gm(flag_ion, parse_from_file = 0)

    def generate_mol_type(self, atom_ids, molecule_id = None, **kwargs):
        mol_type = self.add_molecule_type(molecule_id = molecule_id, **kwargs)
        if kwargs.get('flag_decouple'):
            do_warn('this also decouples the force-field parameters')
            temp_source = copy.deepcopy(self)
        else:
            temp_source = self
        for at in temp_source.get_atoms():
            if at.id in atom_ids:
                mol_type.add2container(at, create=True)
                mol_type.add2container(at.res, replace=-1, create=True)
                mol_type.add2container(at.cg, replace=-1, create=True, db_type=list)
        self.renumber_container('cg', attrib = 'n', id_type = int)
        for int_cont in temp_source.get_interaction_containers():
            flag_ds = False # distance restraints
            if int_cont and int_cont[0].int_type == self.Distance_r:
                do_warn('distance restraint parameters decoupled - only state 0')
                flag_ds = True
                c_index = 0
            for temp_int in int_cont:
                flag_add = 1
                for at in temp_int.atoms:
                    if at.id not in atom_ids:
                        flag_add = 0
                        break
                if flag_add:
                    if flag_ds:
                        new_int_state0 = temp_int.states[0].gr2gm(index = c_index)
                        new_int = self.Interaction(self.Distance_r, atoms = temp_int.atoms, states = [new_int_state0])
                        for temp_state in temp_int.states[1:]:
                            new_int.add_state(temp_state)
                        mol_type.add2container(new_int, create=True, db_type=list, container_name=new_int_state0, flag_class=1)
                        c_index+=1
                    else:
                        gm_int_type = temp_int.states[0]._get_gr2gm_int_type_class()
                        mol_type.add2container(temp_int, create=True, db_type=list, container_name=gm_int_type, flag_class=1)
        if hasattr(temp_source, 'excl_pair'):
            for EP_pair, temp_int in temp_source.excl_pair.items():
                flag_add = 1
                for at in temp_int.atoms:
                    if at.id not in atom_ids:
                        flag_add = 0
                        break
                if flag_add:
                    mol_type.add2container(temp_int, create=True, item_id=EP_pair, replace=-1)
                    for perm_EP_pair in permutations(EP_pair):
                        mol_type.EP_l[perm_EP_pair[0]].add(perm_EP_pair[1])
                    
        mol_type._gm_generate_excl_pairs(reset=True)
        if kwargs.get('flag_renumber_atoms', True):
            mol_type.renumber_container('atoms', attrib = 'gm_id', id_type = int)
        if kwargs.get('flag_renumber_residues', False):
            mol_type.renumber_container('residues', attrib = 'gm_id', id_type = int)
        else:
            for res in mol_type.residues.values():
                res.gm_id = res.id
        return mol_type

# gm2gr
    def generate_sys_atoms_from_mols_types(self, **kwargs):
        sys_atoms = self.get_container('atoms', create = True)
        sys_residues = self.get_container('residues', create = True)
        curr_res = None
        for m in self.molecules:
            del(m.mol_type.ff)
            for i in range(m.n):
                temp_mol_type = copy.deepcopy(m.mol_type)
                temp_mol_type.ff = self
                #temp_mol_type.generate_e_p_l()
                for at in temp_mol_type.get_atoms():
                    at.id = self.find_next_code('atoms')
                    at.gr_id = at.id
                    if at.res != curr_res:
                        at.res.id = self.find_next_code('residues')
                        at.res.gr_id = at.res.id
                        curr_res = at.res
                    self.add2container(at)
                    self.add2container(at.res, replace=-1, create=True)
                    self.add2container(at.cg, replace=-1, create=True, db_type=list)
                self.EP_l.update(temp_mol_type.EP_l)
                for atom_pair in getattr(temp_mol_type, 'excl_pair', []):
                    self.add2container(temp_mol_type.excl_pair[atom_pair], item_id = atom_pair, create = True, replace = -1)
                for int_cont in temp_mol_type.get_interaction_containers():
                    if int_cont and int_cont[0].states[0].__class__ not in (self.PairType, self.ExclusionType) and \
                            hasattr(int_cont[0].states[0], 'fnc') and 'gr_fnc' in int_cont[0].states[0].fnc:
                        for temp_int in int_cont:
                            gr_int_type_class = temp_int.states[0]._get_gm2gr_int_type_class()
                            for i, temp_int_type in enumerate(temp_int.states):
                                temp_gr_int_type = temp_int_type.gm2gr()
                                temp_gr_int_type.cont_pref = self._gr_cont_pref
                                temp_ff_gr_int_type = self.check_eq_int_type(temp_gr_int_type)
                                if not self.check_eq_int_type(temp_ff_gr_int_type):
                                    self.get_intDB().add2container(temp_gr_int_type, create = True)
                                    temp_ff_gr_int_type = temp_gr_int_type
                                temp_int.states[i] = temp_ff_gr_int_type
                            self.add2container(temp_int, create=True, db_type=list, container_name=gr_int_type_class,
                                               flag_class=1)
            m.mol_type.ff = self
        for cg in self.cg:
            cg.update()

    def __get_SOL_gm2gr(self):
        if 'SOLVENTATOM' not in self.undefined_bl:
            self.undefined_bl['SOLVENTATOM'] = []
            sol_m_type = self.molecule_types.get('SOL')
            if sol_m_type:
                self.undefined_bl['SOLVENTATOM'].append('{:}\n'.format(len(sol_m_type.atoms)))
                for i,at in enumerate(sol_m_type.get_atoms()):
                    temp_l = '{:>5} {:>5} {:>4} {:9.5f} {:9.5f}\n'.format(str(i+1), at.name, at.a_type.gr_id, at.m, at.p_ch)
                    self.undefined_bl['SOLVENTATOM'].append(temp_l)
        if 'SOLVENTCONSTR' not in self.undefined_bl:
            self.undefined_bl['SOLVENTCONSTR'] = []
            sol_m_type = self.molecule_types.get('SOL')
            if sol_m_type:
                if hasattr(sol_m_type, 'bonds'):
                    self.undefined_bl['SOLVENTCONSTR'].append('{:}\n'.format(len(sol_m_type.bonds)))
                    for b in sol_m_type.bonds:
                        temp_l = '{:>5}{:>5}{:9.5f}\n'.format(b.atoms[0].id + 1, b.atoms[1].id + 1, b.states[0].p[0])
                        self.undefined_bl['SOLVENTCONSTR'].append(temp_l)
                elif hasattr(sol_m_type, 'settles'):
                    self.undefined_bl['SOLVENTCONSTR'].append('3\n')
                    for i in range(2):
                        temp_l = '{:>5}{:>5} {:9.5f}\n'.format(1, i + 2, sol_m_type.settles[0].states[0].p[0])
                        self.undefined_bl['SOLVENTCONSTR'].append(temp_l)
                    temp_l = '{:>5}{:>5} {:9.5f}\n'.format(2, 3, sol_m_type.settles[0].states[0].p[1])
                    self.undefined_bl['SOLVENTCONSTR'].append(temp_l)

    def top_gm2gr(self, **kwargs):
        if not hasattr(self, 'undefined_bl'):
            self.undefined_bl = OrderedDict()
        self.ff_gm2gr(**kwargs)
        self.__get_SOL_gm2gr()
        self.generate_sys_atoms_from_mols_types(**kwargs)

    def gm2gr(self, **kwargs):
        self.top_gm2gr(**kwargs)

class MolType_g2g:
# gm2gr
    """
    def generate_e_p_l(self):
        db = self.get_container(self.ff.PairType.find_int_container2write(), allow_not_found = True)
        if db:
            for temp_pair in db:
                at_pair = temp_pair.atoms
                at_pair[0].add_pair(at_pair[1], replace = -1)
        self.add_exclusions_neigh(self.nrexcl)
        db = self.get_container(self.ff.ExclusionType.find_int_container2write(), allow_not_found = True)
        if db:
            for temp_pair in db:
                at_pair = temp_pair.atoms
                at_pair[0].add_excl(at_pair[1], replace = -1)
        self._gr_expand_excl_pair()
        self._gr_reduce_excl_pair()
        self._gr_sort_excl_pair()
        self.remove_excl_from_p_l(flag_reduced=True)

    def remove_excl_from_p_l(self, **kwargs):
        if kwargs.get('flag_reduced'):
            N4for = 1
        else:N4for = 2
        for at in self.get_atoms():
            for at2 in at.p_l:
                at_pair = (at, at2)
                for i in range(N4for):
                    try:
                        at_pair[i].e_l.remove(at_pair[(i+1)%2])
                    except ValueError:
                        pass
    """
    pass

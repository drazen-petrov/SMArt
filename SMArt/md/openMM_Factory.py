import numpy as np
from itertools import combinations
from dataclasses import dataclass, fields

from SMArt.md.data_st import MD_Parameters

import openmm as mm
from openmm import app

### non-bonded
# for LJ
# normal vdw for normal pairs
# 1-4 vdw for 1-4 pairs
# 0 for exclusions

# for reaction field
# E = coul + RF + RFC
# for normal pairs and 1-4 pair we need all 3
# if contribution of excluded atoms is included:
#    for exclusions we need RF + RFC
#    for self-pairs we need only RFC

## custom GROMOS RF adds could + RF + RFC for normal pairs only (exclusions and 1-4 pairs will be excluded);
# we need to add coul + RF + RFC for 1-4 pairs, RF + RFC for exclusions and RFC for self-pairs

## openMM native adds coul + RF + RFC for normal;
# if exceptions used (and qq_product added), this keeps the direct coul term only!!!
#    if pairs are implemented as exceptions, we need to add RF + RFC them
# we need to add RF + RFC for exclusions and RFC for self-pairs
# not sure what happens if PME is used...


@dataclass
class WaterModel:
    C6: float
    C12: float
    sigma: float
    epsilon: float
    OH: float
    HH: float
    O_m:float
    H_m:float
    O_pch: float
    H_pch: float
    O_name: str='OW'
    H_name: str='HW'
    
    def get_SE_atom(self, pair_SE):
        sigma = pair_SE[0] * 2 - self.sigma
        epsilon = pair_SE[1] ** 2 / self.epsilon
        return sigma, epsilon
    
    def get_LJ_from_vdw_pair(self, OW_OW_vdw_pair):
        self.C6, self.c12 = OW_OW_vdw_pair.p
        self.sigma, self.epsilon = OW_OW_vdw_pair.convert2se()
        
        
SPC = WaterModel(C6=0.002617346, C12=2.634129e-06, sigma=0.3165648203088436, epsilon=0.6501674826589737, 
                 OH=0.1, HH=0.1632990, O_m=15.99940, H_m=1.00800, O_pch=-0.82, H_pch=0.41)

@dataclass
class ELE_params:
    RF_eps: float
    R_cutoff: float
    eps1: float=1.
    kappa: float=0
    flag_RF_excl: bool=True
    FPEPSI: float=138.9354

    def set_FPEPSI_from_top(self, top):
        self.FPEPSI = float(top._gr_physconst[0])

    def set_Crf(self):
        fac = (2*self.eps1 - 2*self.RF_eps)*(1 +self.kappa*self.R_cutoff)
        fac -= self.RF_eps*((self.kappa*self.R_cutoff)**2)
        fac2 = (self.eps1 + 2*self.RF_eps)*(1 +self.kappa*self.R_cutoff) 
        fac2 += self.RF_eps*((self.kappa*self.R_cutoff)**2)
        self.Crf = fac/fac2 
    
    def set_eps_const_recip(self):
        self.eps_const_recip = self.FPEPSI * self.eps1 

    def set_eps_const_recip_R_3__C_rf_half(self):
        R_cutoff = self.R_cutoff
        self.eps_const_recip_R_3__C_rf_half = self.eps_const_recip * (-0.5 * self.Crf)  / (R_cutoff*R_cutoff*R_cutoff)
    
    def set_eps_const_recip_R_1__C_rf_half__1(self):
        self.eps_const_recip_R_1__C_rf_half__1 = self.eps_const_recip * (0.5 * self.Crf - 1) / self.R_cutoff
    
    def set_all_constants(self, top=None):
        if top:
            self.set_FPEPSI_from_top(top)
        self.set_Crf()
        self.set_eps_const_recip()
        self.set_eps_const_recip_R_3__C_rf_half()
        self.set_eps_const_recip_R_1__C_rf_half__1()


class openMM_Factory:
    class BaseNBFlags:
        @classmethod
        def prep_dict_4_init(cls, init_dict):
            new_dict = dict()
            for field in fields(cls):
                k = field.name
                if k in init_dict:
                    new_dict[k] = init_dict[k]
            return new_dict
    
        def _check_linked_attr(self, parent_attr, child_attr, flag_and):
            assert hasattr(self, parent_attr) and hasattr(self, child_attr)
            parent_attr_value = getattr(self, parent_attr)
            child_attr_value = getattr(self, child_attr)
            if getattr(self, child_attr) is None:
                return getattr(self, parent_attr)
            else:
                if flag_and:
                    return getattr(self, child_attr) and getattr(self, parent_attr)
                else:
                    return getattr(self, child_attr)

        @property
        def _use_Q_14(self):
            return self._check_linked_attr('use_charges', 'use_Q_14', flag_and=True)


    @dataclass
    class Native_NB_Flags(BaseNBFlags):
        use_charges: bool=False
        use_lj_water: bool=False
        use_lj_solute_solute: bool=False
        use_Q_14: bool=None
        use_lj_14: bool=None
        @property
        def _use_lj_14(self):
            return self._check_linked_attr('use_lj_solute_solute', 'use_lj_14', flag_and=False)


    @dataclass
    class GROMOS_NB_Flags(BaseNBFlags):
        use_charges: bool=False # coul + RF + RFC
        use_lj: bool=False
        coul: bool=None
        RF: bool=None
        RFC: bool=None
        use_Q_excl: bool=None # for exclusions, use_charges (e.g. for including RF for exclusions)
        use_Q_14: bool=None # for 14 pair
        @property
        def _coul(self):
            return self._check_linked_attr('use_charges', 'coul', flag_and=True)
        @property
        def _RF(self):
            return self._check_linked_attr('use_charges', 'RF', flag_and=True)
        @property
        def _RFC(self):
            return self._check_linked_attr('use_charges', 'RFC', flag_and=True)
        
        @property
        def _use_Q_excl(self):
            return self._check_linked_attr('use_charges', 'use_Q_excl', flag_and=True)

        def assert_one_term(self):
            assert self.use_lj or self._coul or self._RF or self._RFC

        def assert_self_term(self):
            assert (self.use_lj, self._coul, self._RF, self._RFC) == (False,False,False,True)

    RAD2DEG_2 = (180 / np.pi)**2

    FG_map = dict(FG_BONDS=0, FG_ANGLES=1, FG_IMPROPERS=2, FG_DIHEDRALS=3, FG_NB_NATIVE_GROUP=4)
    FG_NB_CUSTOM_GROUP_START = 5
    FG_map2 = dict(FG_NB_LJ_CUSTOM_GROUP=20, FG_NB_RF_CUSTOM_GROUP=21, FG_B_LJ14_CUSTOM_GROUP=22)
    
    exclusions_Q = [0.]
    exclusions_LJ = [0.5, 0.]
    
    
    def _get_ELE_params(self, MD_params):
        if MD_params is None:
            MD_params={}
        if isinstance(MD_params, dict):
            R_cutoff = MD_params.get('R_cutoff')
            if R_cutoff is None:
                print('setting R_cutoff = 1.4')
                R_cutoff = 1.4
            RF_eps = MD_params.get('RF_eps')
            if RF_eps is None:
                print('setting RF_eps = 62.')
                RF_eps = 62.
            eps1 = MD_params.get('eps1')
            if eps1 is None:
                print('setting eps1 = 1.')
                eps1 = 1.
            kappa = MD_params.get('kappa')
            if kappa is None:
                print('setting kappa = 0')
                kappa = 0
            flag_RF_excl = MD_params.get('flag_RF_excl')
            if flag_RF_excl is None:
                print('setting flag_RF_excl = True')
                flag_RF_excl = True
            RF_params = ELE_params(RF_eps, R_cutoff, eps1, kappa, flag_RF_excl)
            return RF_params
        elif isinstance(MD_params, ELE_params):
            return MD_params
        elif isinstance(MD_params, MD_Parameters):
            assert MD_params.NLRELE == 1 # this is RF
            self.RF_flag = True
            eps1 = 1. # no idea where is this defined... so it's HARDCODED here...
            print('setting eps1 = 1.')
            return ELE_params(MD_params.EPSRF, MD_params.RCRF, eps1, MD_params.KAPPA, MD_params.NSLFEXCL)
        else:
            raise Exception('MD_params of unknown type' )
        return
    
    def generate_water_top(self, N_water, water_res_name='SOLV', **kwargs):
        for at in self.top.get_atoms():
            N_states = len(at.m_states)
            break
        top_w = self.top.reduce()
        self._top_water = top_w
        self.top = top_w
        OW_atype = self._get_OW_atype(**kwargs)
        HW_atype = self._get_HW_atype(**kwargs)
        OH_constraint_type = top_w.ConstraintType(params=[self.water_model.OH])
        HH_constraint_type = top_w.ConstraintType(params=[self.water_model.HH])
        top_w._add_excl_pair_types()
        excl_type = top_w.get_intDB().excl_pair['excl']
        self.N_water=N_water
        for i in range(N_water):
            new_res = top_w.add_residue(res_name=water_res_name)
            OW = top_w.add_atom(atom_name=self.water_model.O_name)
            HW1 = top_w.add_atom(atom_name=self.water_model.H_name + '1')
            HW2 = top_w.add_atom(atom_name=self.water_model.H_name + '2')
            for temp_at in (OW, HW1, HW2):
                temp_at.gr_id = temp_at.id
                new_res.add_atom(temp_at)
                temp_at._SOLV_flag = True
            # add atom params
            OW.m_states = [self.water_model.O_m] * N_states
            OW.p_ch_states = [self.water_model.O_pch] * N_states
            OW.a_type_states = [OW_atype] * N_states
            #OW.p_ch = self.water_model.O_pch
            for HW in (HW1, HW2):
                HW.m_states = [self.water_model.H_m] * N_states
                HW.p_ch_states = [self.water_model.H_pch] * N_states
                #HW.p_ch = self.water_model.H_pch
                HW.a_type_states = [HW_atype] * N_states
            # skip making a charge group
            # add constraints
            OH1_c = top_w.Interaction(top_w.ConstraintType, atoms=[OW, HW1], states=[OH_constraint_type]*N_states)
            OH2_c = top_w.Interaction(top_w.ConstraintType, atoms=[OW, HW2], states=[OH_constraint_type]*N_states)
            HH_c = top_w.Interaction(top_w.ConstraintType, atoms=[HW1, HW2], states=[HH_constraint_type]*N_states)
            for temp_constraint in (OH1_c, OH2_c, HH_c):
                top_w.add2container(temp_constraint, create=True, db_type=list)
            # add exclusions
            OH1_excl = top_w.Interaction(top_w.ExclusionPairType, atoms=[OW, HW1], states=[excl_type]*N_states)
            OH2_excl = top_w.Interaction(top_w.ExclusionPairType, atoms=[OW, HW2], states=[excl_type]*N_states)
            HH_excl = top_w.Interaction(top_w.ExclusionPairType, atoms=[HW1, HW2], states=[excl_type]*N_states)
            for temp_excl in (OH1_excl, OH2_excl, HH_excl):
                top_w.add2container(temp_excl, create=True, item_id=frozenset(temp_excl.atoms))
        return

    @staticmethod
    def _check_if_water_atom(at):
        return hasattr(at, '_SOLV_flag') and at._SOLV_flag
    
    def _get_solute_atoms(self):
        for at in self.top.get_atoms():
            if not self._check_if_water_atom(at):
                yield at
    
    def __init__(self, top=None, MD_params=None, top_state=None, water_model=SPC, 
                 RF_flag=True, NB_cutoff=None, **kwargs):
        self.FG_gr_nb_count = self.FG_NB_CUSTOM_GROUP_START
        if top_state:
            top.set_top_state(top_state)
        self._top = top
        self.top = top.reduce()
        self.MD_params = MD_params
        self.RF_flag = RF_flag
        self.ELE_params = self._get_ELE_params(MD_params)
        self.ELE_params.set_all_constants(top=top)
        if NB_cutoff:
            self.NB_cutoff = NB_cutoff
        else:
            self.NB_cutoff = self.ELE_params.R_cutoff
        self.exclusions_params = list(self._get_openMM_NB_params(*self.exclusions_Q, *self.exclusions_LJ))
        self.exclusions_params_list = self._get_openMM_NB_params(*self.exclusions_Q, *self.exclusions_LJ, replace=True)
        self.exclusions_params_list = list(self.exclusions_params_list)
        self.water_model = water_model
    
    ########### Non-bonded energy ###########
    def _get_bonded_container(self, container_name):
        bi_container = self.top.get_container(container_name, allow_not_found=True)
        if bi_container:
            return bi_container
        else:
            return list()
    
    def _generate_openmm_gromos_bonds(self):
        energy = 'k_bond_quarter * (r*r - r0_2)^2'
        force = mm.CustomBondForce(energy)
        force.addPerBondParameter('k_bond_quarter') # force constant (already multiplied with 0.25)
        force.addPerBondParameter('r0_2') # eq values (already squared)
        force.setName('GromosBond')
        force.setForceGroup(self.FG_map['FG_BONDS'])
        return force

    def add_constraints(self, bond_constraints=True):
        if bond_constraints:
            for bi in self._get_bonded_container('bonds'): # bi - bonded interaction
                atoms_idx = [at._openmm_id for at in bi.atoms]
                state = self.top.get_state(bi.states)
                self.system.addConstraint(*(atoms_idx + [state.p[2]]))
        for bi in self._get_bonded_container('constraints'):
            atoms_idx = [at._openmm_id for at in bi.atoms]
            state = self.top.get_state(bi.states)
            self.system.addConstraint(*(atoms_idx + list(state.p)))
        return

    def add_bonds(self):
        force = self._generate_openmm_gromos_bonds()
        for bi in self._get_bonded_container('bonds'): # bi - bonded interaction
            atoms_idx = [at._openmm_id for at in bi.atoms]
            state = self.top.get_state(bi.states)
            k_bond_quarter = state.p[0] * 0.25
            r0_2 = state.p[2] * state.p[2]
            b_params = [(k_bond_quarter, r0_2)]
            force.addBond(*(atoms_idx + b_params))
        self.system.addForce(force)

    def _generate_openmm_gromos_angles(self):
        energy = 'k_ang_half*(cos(theta)-cos_theta0)^2'
        force = mm.CustomAngleForce(energy)
        force.addPerAngleParameter('k_ang_half') # force constant (already multiplied with 0.5)
        force.addPerAngleParameter('cos_theta0') # eq values (cos of it)
        force.setName('GromosAngle')
        force.setForceGroup(self.FG_map['FG_ANGLES'])
        return force

    def add_angles(self):
        force = self._generate_openmm_gromos_angles()
        for bi in self._get_bonded_container('angles'): # bi - bonded interaction
            atoms_idx = [at._openmm_id for at in bi.atoms]
            state = self.top.get_state(bi.states)
            k_ang_half = state.p[0] * 0.5
            cos_theta0 = np.cos(np.radians(state.p[2]))
            b_params = [(k_ang_half, cos_theta0)]
            force.addAngle(*(atoms_idx + b_params))
        self.system.addForce(force)

    def _generate_openmm_gromos_impropers(self):
        energy = 'k_imp_rad_half * (theta - theta0)^2' # d_thera is in rad!
        force = mm.CustomTorsionForce(energy)
        force.addPerTorsionParameter('k_imp_rad_half') # force constant (already multiplied with 0.5)
        force.addPerTorsionParameter('theta0') # eq value (in rad)
        force.setName('GromosImproper')
        force.setForceGroup(self.FG_map['FG_IMPROPERS'])
        return force

    def add_impropers(self):
        force = self._generate_openmm_gromos_impropers()
        for bi in self._get_bonded_container('impropers'): # bi - bonded interaction
            atoms_idx = [at._openmm_id for at in bi.atoms]
            state = self.top.get_state(bi.states)
            k_imp_rad_half = state.p[0] * 0.5 * self.RAD2DEG_2 # convert from deg^-2 to rad^-2
            theta0 = np.radians(state.p[1])
            b_params = [(k_imp_rad_half, theta0)]
            force.addTorsion(*(atoms_idx + b_params))
        self.system.addForce(force)

    def _generate_openmm_gromos_dihedrals(self):
        energy = 'k_dih*(1 + cos_delta0 * cos(m*theta))'
        force = mm.CustomTorsionForce(energy)
        force.addPerTorsionParameter('k_dih') # force constant
        force.addPerTorsionParameter('cos_delta0') # shift (already cos of it)
        force.addPerTorsionParameter('m') # multiplicity
        force.setName('GromosDihedral')
        force.setForceGroup(self.FG_map['FG_DIHEDRALS'])
        return force

    def add_dihedrals(self):
        force = self._generate_openmm_gromos_dihedrals()
        for bi in self._get_bonded_container('dihedrals'): # bi - bonded interaction
            atoms_idx = [at._openmm_id for at in bi.atoms]
            state = self.top.get_state(bi.states)
            k_dih = state.p[0]
            cos_delta0 = np.cos(np.radians(state.p[1]))
            m = state.p[2]
            b_params = [(k_dih, cos_delta0, m)]
            force.addTorsion(*(atoms_idx + b_params))
        self.system.addForce(force)

    ########### Non-bonded energy ###########
    ## helper functions
    @staticmethod
    def _get_openMM_NB_params(qq_prod, sigma, epsilon, replace=None):
        params = (qq_prod, sigma, epsilon)
        if replace is not None:
            params += (replace, )
        return params

    def get_SE_param__LJ(self, a_type1, a_type2, flag_pair=False):
        vdw_pair = self.top.get_vdw_pair(a_type1, a_type2)
        if flag_pair:
            return vdw_pair[1].convert2se()
        else:
            return vdw_pair[0].convert2se()

    def _openMM_SE_convention(self, sigma, epsilon):
        if epsilon==0.:
            return tuple(self.exclusions_LJ)
        else:
            return sigma, epsilon

    def get_SE_param__LJ_openMM(self, a_type1, a_type2, flag_pair=False, **kwargs):
        SE = self.get_SE_param__LJ(a_type1, a_type2, flag_pair=flag_pair, **kwargs)
        return self._openMM_SE_convention(*SE)

    def _get_OW_atype(self, **kwargs):
        OW_atype_id = kwargs.get('OW_atype_id')
        if OW_atype_id:
            return self.top.get_intDB().a_type[OW_atype_id]
        else:
            for atype in self.top.get_intDB().a_type.values():
                if atype.name=='OW':
                    return atype

    def _get_HW_atype(self, **kwargs):
        HW_atype_id = kwargs.get('HW_atype_id')
        if HW_atype_id:
            return self.top.get_intDB().a_type[HW_atype_id]
        else:
            for atype in self.top.get_intDB().a_type.values():
                if atype.name=='H':
                    return atype

    ## set pbc parameters (and cutoff) to the forces
    def __set_pbc_nb_force(self, temp_nb_force, pbc, flag_native, openMM_cutoff=None):
        if pbc:
            if flag_native:
                temp_nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            else:
                temp_nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        else:
            if flag_native:
                temp_nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            else:
                temp_nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        if openMM_cutoff is None:    
            temp_nb_force.setCutoffDistance(self.NB_cutoff)
        else:
            temp_nb_force.setCutoffDistance(openMM_cutoff)
        return

    def __set_pbc_b_force(self, temp_b_force, pbc):
        if pbc:
            temp_b_force.setUsesPeriodicBoundaryConditions(True)
        else:
            temp_b_force.setUsesPeriodicBoundaryConditions(False)
        return
    
    def set_pbc(self, pbc=True, nb_force=None, custom_nb_forces=None, custom_b_forces=None, openMM_cutoff=None):
        """
        sets parameters like PBC and cutoff
        """
        if nb_force:
            if self.RF_flag:
                print('setting eps =', self.ELE_params.RF_eps)
                nb_force.setReactionFieldDielectric(self.ELE_params.RF_eps)
            self.__set_pbc_nb_force(nb_force, pbc, True, openMM_cutoff=openMM_cutoff)# True -> to treat it as a native openMM NB force
        if custom_nb_forces:
            for temp_nb_force in custom_nb_forces:
                self.__set_pbc_nb_force(temp_nb_force, pbc, False, openMM_cutoff=openMM_cutoff)# False -> to treat it as a custom openMM NB force
        if custom_b_forces:
            for temp_b_force in custom_b_forces:
                self.__set_pbc_b_force(temp_b_force, pbc)

    ## openMM native non-bonded interactions
    def _generate_openmm_native_nonbonded(self, flag_add_params=False, **kwargs):
        """
        generats native openMM NB force object
        :param flag_add_params: if True, adds per atom parameters
        :kwargs:
            flag_DispersionCorrection (False): used for setUseDispersionCorrection
        :return: nb_force object
        """
        nb_force = mm.NonbondedForce()
        flag_DispersionCorrection = kwargs.get('flag_DispersionCorrection', False)
        nb_force.setUseDispersionCorrection(flag_DispersionCorrection)
        nb_force.setForceGroup(self.FG_map['FG_NB_NATIVE_GROUP'])
        if flag_add_params:
            self.add_particles2nb_force(nb_force, **kwargs)
        return nb_force
    
    def add_particles2nb_force(self, force, **kwargs):
        """
        addes per particle parameters to the openMM native non-bonded force object
        number of water molecules takes from self.N_water 
        and NB_flags (what to add to native openMM NB; e.g. use_charges or use_lj_water) from self.native_NB_flags
        if use_lj_solute_solute == True; all vs all pairs of solute atoms will be used for exceptions!
        and in this case all non-bonded interactions have to be defined via the openMM native NB force object
        
        :param force: - openMM NB force object
        :param kwargs:
            passed on to self._get_OW_atype
            passed on to self._get_HW_atype
        :return: None
        """
        if self._native_NB_flags.use_charges and self._native_NB_flags.use_lj_solute_solute:
            print('WARNING: use_lj_solute_solute is implemented via exceptions for all vs all solute atoms, which removes RF contribution!')
            assert kwargs.get('use_charges__lj_solute_solute')
        if self._native_NB_flags.use_lj_water:
            assert self._native_NB_flags.use_lj_solute_solute
            OW_atype = self._get_OW_atype(**kwargs)
            HW_atype = self._get_HW_atype(**kwargs)
        FF = self.top.get_intDB()
        for at in self.top.get_atoms():
            p_ch, sigma, epsilon = self.exclusions_params# sigma=0.5, epsilon=0 -> ignore LJ potential
            if self._native_NB_flags.use_charges:
                p_ch = self.top.get_state(at.p_ch_states)
            if self._native_NB_flags.use_lj_water: # calculate "fake" sigma and epsilon that it gets correct LJ potential with OW atom type
                atype = self.top.get_state(at.a_type_states)
                SE_pair = self.get_SE_param__LJ_openMM(atype, OW_atype, flag_pair=False)
                SE_atom = self.water_model.get_SE_atom(SE_pair)
                sigma, epsilon = self._openMM_SE_convention(*SE_atom) # make sure to have it openMM compatible (for 0)
            force.addParticle(p_ch, sigma, epsilon)
        if self._native_NB_flags.use_lj_solute_solute:
            for at_pair in combinations(self._get_solute_atoms(), 2):
                atom_idx = [at._openmm_id for at in at_pair]
                a_type1, a_type2 = (self.top.get_state(at.a_type_states) for at in at_pair)
                SE_pair = self.get_SE_param__LJ_openMM(a_type1, a_type2, flag_pair=False)
                qq_product = 0.
                if self._native_NB_flags.use_charges:
                    qq_product = self.top.get_state(at_pair[0].p_ch_states) * self.top.get_state(at_pair[1].p_ch_states)
                pair_param_list = list(self._get_openMM_NB_params(qq_product, *SE_pair, replace=True))
                force.addException(*(atom_idx + pair_param_list))
        return


    ## GROMOS custom energy terms
    ## GROMOS custom LJ potential
    LJ_energy_expression = 'c12*r_12 - c6*r_6; '
    LJ_C6_C12_expression = 'c6 = get_c6(atype1, atype2); c12 = get_c12(atype1, atype2); '
    ## GROMOS custom RF potential
    coul_direct_expression = 'eps_const_recip * qq_prod * r_1; ' # this should be done per particle
    RF_dist_dep_expression = 'eps_const_recip_R_3__C_rf_half * qq_prod * r2; ' # this should be done per pair (as bond)
    RF_dist_indep_expression = 'eps_const_recip_R_1__C_rf_half__1 * qq_prod; ' # this should be done per pair (as bond) - include self term as well!
    qq_prod_expression = 'qq_prod = q1*q2; '

    #qq_prod = q1*q2
    #eps_const_recip = 1 / (4 pi eps0 eps1)
    #eps_const_recip_R_3 = eps_const_recip * R_cutoff^-3
    #C_rf = ((2eps1 - 2eps2)*(1+kappa*R_cutoff) - eps2*(kappa*R_cutoff)^2) / ((esp1 + 2eps2)*(1+kappa*R_cutoff) + eps2*(kappa*R_cutoff)^2)
    #C_rf_half = -1/2 * C_rf 
    #eps_const_recip_R_3__C_rf_half = eps_const_recip_R_3 * C_rf_half
    #eps_const_recip_R_1 = eps_const_recip * R_cutoff^-1
    #C_rf_half__1 = C_rf_half - 1
    #eps_const_recip_R_1__C_rf_half__1 = eps_const_recip_R_1 * C_rf_half__1
    # usually: eps1=1; kappa=0;

    @staticmethod
    def _get_r_expression(energy_expression, **kwargs):
        #kwargs r_12=False, r_6=False, r_2=False, r_1=False, r2=False
        r_map = dict(kwargs)
        for r_term in ['r_12', 'r_6', 'r_2', 'r2', 'r_1']:
            r_map[r_term] = kwargs.get(r_term, False)
            if r_term in energy_expression:
                r_map[r_term] = True
        r_expression = ''
        if r_map['r_12']:
            r_map['r_6'] = True
            r_expression += 'r_12 = r_6*r_6; '
        if r_map['r_6']:
            r_map['r_2'] = True
            r_expression += 'r_6 = r_2*r_2*r_2; '
        if r_map['r_2']:
            if not r_map['r_1'] and r_map['r2']:
                r_expression += 'r_2 = 1 / r2; '
            else:
                r_map['r_1'] = True
                r_expression += 'r_2 = r_1*r_1; '
        if r_map['r_1']:
            r_expression += 'r_1 = 1/r; '
        if r_map['r2']:
            r_expression += 'r2 = r*r; '
        return r_expression

    def _get_GROMOS_NB_energy_expression(self, LJ=False, coul=False, RF=False, RFC=False):
        assert sum((LJ, coul, RF, RFC))>0
        if sum((LJ, coul, RF, RFC))>1:
            energy = ''
            ene_terms = ''
            if LJ:
                energy += 'lj + '
                ene_terms += 'lj = {LJ_energy_expression}'
            if coul:
                energy += 'coul + '
                ene_terms += 'coul = {coul_direct_expression}'
            if RF:
                energy += 'RF + '
                ene_terms += 'RF = {RF_dist_dep_expression}'
            if RFC:
                energy += 'RFC + '
                ene_terms += 'RFC = {RF_dist_indep_expression}'
            energy = energy[:-3] + '; ' + ene_terms
        else:
            if LJ:energy = '{LJ_energy_expression}'
            if coul:energy = '{coul_direct_expression}'
            if RF:energy = '{RF_dist_dep_expression}'
            if RFC:energy = '{RF_dist_indep_expression}'
        energy = energy.format(LJ_energy_expression=self.LJ_energy_expression,
                               coul_direct_expression=self.coul_direct_expression,
                               RF_dist_dep_expression=self.RF_dist_dep_expression,
                               RF_dist_indep_expression=self.RF_dist_indep_expression)
        return energy

    def __generate_FG_gr_name(self, gromos_NB_flags, FG_template='FG_NB{:}_CUSTOM', **kwargs):
        k = str(FG_template)
        addition = ''
        if gromos_NB_flags.use_lj:
            addition+='_LJ'
        if gromos_NB_flags._coul:
            addition+='_coul'
        if gromos_NB_flags._RF:
            addition+='_RF'
        if gromos_NB_flags._RFC:
            addition+='_RFC'
        FG_name_base = k.format(addition)
        suf = kwargs.get('FG_name_suf')
        if not suf:
            if kwargs.get('flag_excl_14'):
                suf = '_excl_14'
        if suf:
            FG_name_base += suf
        FG_name = str(FG_name_base)
        c = 1
        while FG_name in self.FG_map:
            c+=1
            FG_name = FG_name_base + '_' + str(c)
        return FG_name
    
    def add_force_group(self, force, gromos_NB_flags=None, **kwargs):
        """
        adds force groups number to the force object and updates the counter if needed
        :param force: openMM force object
        :param gromos_NB_flags: flags for energy terms (if name and num is not defined)
        :kwargs:
            FG_number: defines which number is used for force group
            FG_name: defines which name is used for force group (for FG_map)
            suf: adds a user specified sufix to the FG_name
            flag_excl_14: adds "_excl_14" sufix to the FG_name (for exclusions / 14 pairs)
        :return None
        """
        FG_number = kwargs.get('FG_number')
        if not FG_number:
            FG_number = self.FG_gr_nb_count
            self.FG_gr_nb_count += 1
        assert FG_number not in set(self.FG_map.values())
        FG_name = kwargs.get('FG_name')
        if not FG_name:
            FG_name = self.__generate_FG_gr_name(gromos_NB_flags, **kwargs)
        self.FG_map[FG_name] = FG_number
        force.setForceGroup(FG_number)
    
    def _generate_openmm_gromos_nonbonded(self, gromos_NB_flags, flag_add_params=False, flag_add_FG=True, **kwargs):
        """
        defines energy expression and creates mm.CustomNonbondedForce object
        :param gromos_NB_flags: defines LJ, coul, RF and RFC flags (which terms to include)
        :param flag_add_params: if True, adds parameters per particle parameters
        :param flag_add_FG: if True, adds force group to the object
        :kwargs:
            kwargs for force group:
                FG_number: defines which number is used for force group
                FG_name: defines which name is used for force group (for FG_map)
            kwargs for add_force_group function
        :return mm.CustomNonbondedForce object
        """
        # define energy expression
        #LJ = gromos_NB_flags.use_lj
        # coul = gromos_NB_flags._coul
        # RF = gromos_NB_flags._RF
        # RFC = gromos_NB_flags._RFC
        gromos_NB_flags.assert_one_term() # at least one of the terms has to be True
        energy = self._get_GROMOS_NB_energy_expression(gromos_NB_flags.use_lj, gromos_NB_flags._coul, 
                                                       gromos_NB_flags._RF, gromos_NB_flags._RFC)
        if gromos_NB_flags.use_lj:
            energy += self.LJ_C6_C12_expression
        if gromos_NB_flags.use_charges:
            energy +=  self.qq_prod_expression

        # add definitions for r terms
        r_expression = self._get_r_expression(energy)
        energy += r_expression
        # generate the force object
        gromos_nb_force = mm.CustomNonbondedForce(energy)

        # add parameters
        # LJ parameters
        if gromos_NB_flags.use_lj:
            FF = self.top.get_intDB()
            ff_a_type = self.top.get_intDB().a_type
            N_a_types = len(ff_a_type)
            c6_matrix = np.empty((N_a_types, N_a_types))
            c12_matrix = np.empty((N_a_types, N_a_types))
            for i, at_type_i in enumerate(ff_a_type.values()):
                for j, at_type_j in enumerate(ff_a_type.values()):
                    vdw_pair = self.top.get_vdw_pair(at_type_i, at_type_j)
                    c6_matrix[i,j] = vdw_pair[0].p[0]
                    c12_matrix[i,j] = vdw_pair[0].p[1]
            gromos_nb_force.addTabulatedFunction('get_c6', mm.Discrete2DFunction(N_a_types, N_a_types, c6_matrix.flatten()))
            gromos_nb_force.addTabulatedFunction('get_c12', mm.Discrete2DFunction(N_a_types, N_a_types, c12_matrix.flatten()))
            gromos_nb_force.addPerParticleParameter('atype')
        # RF parameters
        if gromos_NB_flags.use_charges:
            gromos_nb_force.addPerParticleParameter('q')
            if gromos_NB_flags._coul:
                gromos_nb_force.addGlobalParameter('eps_const_recip', self.ELE_params.eps_const_recip)
            if gromos_NB_flags._RF:
                eps_const_recip_R_3__C_rf_half = self.ELE_params.eps_const_recip_R_3__C_rf_half
                gromos_nb_force.addGlobalParameter('eps_const_recip_R_3__C_rf_half', eps_const_recip_R_3__C_rf_half)
            if gromos_NB_flags._RFC:
                eps_const_recip_R_1__C_rf_half__1 = self.ELE_params.eps_const_recip_R_1__C_rf_half__1
                gromos_nb_force.addGlobalParameter('eps_const_recip_R_1__C_rf_half__1', eps_const_recip_R_1__C_rf_half__1)

        if flag_add_params:
            self.add_GROMOS_NB_atom_parameters(gromos_nb_force, gromos_NB_flags)
        
        if flag_add_FG:
            self.add_force_group(gromos_nb_force, gromos_NB_flags, **kwargs)
        return gromos_nb_force

    def add_GROMOS_NB_atom_parameters(self, gromos_nb_force, gromos_NB_flags):
        for at in self.top.get_atoms():
            params = ()
            if gromos_NB_flags.use_lj:
                at_type = self.top.get_state(at.a_type_states)
                params += (int(at_type.gr_id) - 1, ) # atom type number (should start from 0)
            if gromos_NB_flags.use_charges:
                p_ch = self.top.get_state(at.p_ch_states)
                params += (p_ch, )
            gromos_nb_force.addParticle(params)
        return

    

##################################################################################################################
    @staticmethod
    def _fix_half_self_RFC_energy(energy):
        return energy.replace('eps_const_recip_R_1__C_rf_half__1', 'half_eps_const_recip_R_1__C_rf_half__1')

    def _generate_openmm_gromos_excl_pairs(self, gromos_NB_flags, flag_add_params=False, flag_self_term=False,
                                           flag_add_FG=True, **kwargs):
        """
        defines energy expression for excl and or 14 pairs, and creates mm.CustomBondForce object
        :param gromos_NB_flags: defines LJ, coul, RF and RFC flags (which terms to include)
        :param flag_add_params: if True, adds parameters for individual bonds (pairs/exclusions)
        :param flag_self_term: if True, self-pairs added (particles with themselves)
        :param flag_add_FG: if True, adds force group to the object
        :kwargs:
            half_self_RFC: if True, self term is devided by 2
            kwargs for force group:
                FG_number: defines which number is used for force group
                FG_name: defines which name is used for force group (for FG_map)
            kwargs for add_force_group function
            kwargs for add_gromos_excl_pairs_parameters function
        :return mm.CustomBondForce object
        """
        gromos_NB_flags.assert_one_term() # at least one of the terms has to be True
        #assert sum((LJ, coul, RF, RFC))>0
        temp_kwargs = dict(kwargs)
        
        energy = self._get_GROMOS_NB_energy_expression(gromos_NB_flags.use_lj, gromos_NB_flags._coul, 
                                                       gromos_NB_flags._RF, gromos_NB_flags._RFC)
        if gromos_NB_flags._RFC and flag_self_term and kwargs.get('half_self_RFC', True):
            energy = self._fix_half_self_RFC_energy(energy)

        r_expression = self._get_r_expression(energy)
        energy += r_expression

        force = mm.CustomBondForce(energy)
        if gromos_NB_flags.use_lj:
            force.addPerBondParameter('c6')
            force.addPerBondParameter('c12')
        if gromos_NB_flags.use_charges:
            force.addPerBondParameter('qq_prod')
        if gromos_NB_flags._coul:
            force.addGlobalParameter('eps_const_recip', self.ELE_params.eps_const_recip)
        if gromos_NB_flags._RF:
            force.addGlobalParameter('eps_const_recip_R_3__C_rf_half', self.ELE_params.eps_const_recip_R_3__C_rf_half)
        if gromos_NB_flags._RFC:
            eps_const_recip_R_1__C_rf_half__1 = self.ELE_params.eps_const_recip_R_1__C_rf_half__1
            if flag_self_term and kwargs.get('half_self_RFC', True):
                eps_const_recip_R_1__C_rf_half__1 = 0.5 * eps_const_recip_R_1__C_rf_half__1
                force.addGlobalParameter('half_eps_const_recip_R_1__C_rf_half__1', eps_const_recip_R_1__C_rf_half__1)
            else:
                force.addGlobalParameter('eps_const_recip_R_1__C_rf_half__1', eps_const_recip_R_1__C_rf_half__1)
        
        if flag_add_params:
            if flag_self_term:
                gromos_NB_flags.assert_self_term() # (LJ, coul, RF, RFC) == (False, False, False, True)
            self.add_gromos_excl_pairs_parameters(force, gromos_NB_flags, flag_self_term=flag_self_term, **kwargs)
        
        if flag_add_FG:
            temp_kwargs = dict(kwargs)
            if 'flag_excl_14' not in temp_kwargs:
                temp_kwargs['flag_excl_14'] = True
            self.add_force_group(force, gromos_NB_flags, **temp_kwargs)

        return force

    def add_gromos_excl_pairs_parameters(self, force, gromos_NB_flags, flag_self_term=False, **kwargs):
        """
        adds parameters for each pair of atoms (exclusions and 14 pairs)
        :param force: openMM force object
        :param gromos_NB_flags: defines LJ, coul, RF and RFC flags (which terms to include) - it's called here NB, but it's implemented as bonds!
        :param flag_self_term: if True, self-pairs added (particles with themselves)
        :return: None
        """
        FF = self.top.get_intDB()
        for excl_pair in self.top.get_excl_pair().values():
            atom_idx = sorted([excl_pair.atoms[0]._openmm_id, excl_pair.atoms[1]._openmm_id])
            state = self.top.get_state(excl_pair.states)
            if gromos_NB_flags.use_charges:
                q1 = self.top.get_state(excl_pair.atoms[0].p_ch_states)
                q2 = self.top.get_state(excl_pair.atoms[1].p_ch_states)
                qq_prod = q1*q2
            params = ()
            if state.p[0] == 1:# this is exclusion
                if gromos_NB_flags.use_lj:
                    params += tuple(self.exclusions_LJ)
                elif gromos_NB_flags._use_Q_excl:
                    params = (qq_prod, )
            elif state.p[0] == 2:# this is pair
                if gromos_NB_flags.use_lj:
                    a_type1 = self.top.get_state(excl_pair.atoms[0].a_type_states)
                    a_type2 = self.top.get_state(excl_pair.atoms[1].a_type_states)
                    vdw_pair = self.top.get_vdw_pair(a_type1, a_type2)
                    params += tuple(vdw_pair[1].p)
                if gromos_NB_flags._use_Q_14:
                    params +=  (qq_prod, )
            if params and params!=tuple(self.exclusions_LJ):
                force.addBond(*atom_idx, params)
        if flag_self_term:
            #this assumes (LJ, coul, RF, RFC) == (False, False, False, True)
            gromos_NB_flags.assert_self_term()
            for at_idx, at in enumerate(self.top.get_atoms()):
                q = self.top.get_state(at.p_ch_states)
                force.addBond(at_idx, at_idx, (q*q,))
        return

    ##### for native NB: use_charges=False, use_lj_solute_solute=False, use_lj_water=False
    # if use charges - this is for all (solute-solute, solute-water, water-water)
    # if use_lj_solute_solute - this works only with exceptions
    # if use_lj_water - this covers solute-water and water-water
    # 1-4 pairs could be done with exceptions
    
    ##### for GROMOS NB
    # any part of LJ or coul for solute-solute; solute-water, water-water can be split
    # 1-4 pairs (both LJ and coul) can only be done with bonds!
    
    ## this has to be in, regardless of how other things are
    # RF contribution of excluded atoms included with bonds!
    # RF distance independent contribution (also including with bonds!)

    def create_openMM_system(self):
        """generates openMM system object and stores it in self.system"""
        self.system = mm.System()
    
    def add_particles(self):
        "add atoms/particles to the system (masses)"
        for at_idx, at in enumerate(self.top.get_atoms()):
            m = self.top.get_state(at.m_states)
            self.system.addParticle(m)
            at._openmm_id = at_idx
        return

    # add exclusions / pairs
    def add_exclusions_pairs(self, nb_force, custom_nb_forces, **kwargs):
        """
        adds exclusions / exceptions to all non-bonded forces in the system
        :param nb_force: native openMM NB force object
        :param custom_nb_forces: list of custom NB force objects
        :kwargs:
        :return: None
        """
        for excl_pair in self.top.get_excl_pair().values():
            atom_idx = sorted([excl_pair.atoms[0]._openmm_id, excl_pair.atoms[1]._openmm_id])
            state = self.top.get_state(excl_pair.states)
            if state.p[0] == 1:# this is exclusion
                nb_force.addException(*(atom_idx + self.exclusions_params_list))
                for temp_nb_force in custom_nb_forces:
                    temp_nb_force.addExclusion(*atom_idx)
            elif state.p[0] == 2:# this is pair
                # get sigma and eps for the pair
                a_type1 = self.top.get_state(excl_pair.atoms[0].a_type_states)
                a_type2 = self.top.get_state(excl_pair.atoms[1].a_type_states)
                SE_pair = self.get_SE_param__LJ_openMM(a_type1, a_type2, flag_pair=True)
                q1 = self.top.get_state(excl_pair.atoms[0].p_ch_states)
                q2 = self.top.get_state(excl_pair.atoms[1].p_ch_states)
                # add exclusion with normal q*q, 1-4 pair sigma and eps, and True for replace
                coul_params, LJ_params = self.exclusions_Q, self.exclusions_LJ # default params for exclusions
                if self._native_NB_flags._use_Q_14:
                    coul_params = [q1*q2]
                if self._native_NB_flags._use_lj_14:
                    LJ_params = list(SE_pair)
                pair_param_list = self._get_openMM_NB_params(*coul_params, *LJ_params, replace=True)
                nb_force.addException(*(atom_idx + list(pair_param_list)))
                # also exclude from lj_force
                for temp_nb_force in custom_nb_forces:
                    temp_nb_force.addExclusion(*atom_idx)
        return
    
    def make_openmm_system(self, N_water=0, bond_constraints=True, pbc=True, 
                           native_NB_flags=None, gromos_NB_flags=None, gromos_B_flags=None, **kwargs):
        """
        generates openMM system object
        :param N_water: number of water molecules
        :param bond_constraints: if True, constrain bonds
        :param native_NB_flags: defines which energy terms are used for native openMM NB force
        :param gromos_NB_flags: list of dicts that define which energy terms are used for custom NB force
        :param gromos_B_flags: list of dicts define which energy terms are used for custom bonded force (used for excl and 14 pairs)
        :kwargs:
            used for self.add_exclusions_pairs function
        :return: None
        """
        if native_NB_flags is None:
            native_NB_flags = {}
        self._native_NB_flags = self.Native_NB_Flags(**native_NB_flags)
        if N_water:
            self.N_water = N_water
            self.generate_water_top(N_water)

        self.create_openMM_system()
        ### add atoms (masses)
        self.add_particles()

        ### non-bonded
        # open MM native non-bonded interaction
        nb_force = self._generate_openmm_native_nonbonded(flag_add_params=True)

        # general non-bonded interaction (GROMOS form)
        self.FG_gr_nb_count = self.FG_NB_CUSTOM_GROUP_START
        custom_nb_forces = []
        if gromos_NB_flags is None:
            gromos_NB_flags = []
        for temp_gromos_NB_dict in gromos_NB_flags:
            if not temp_gromos_NB_dict.get('flag_add_params'):
                temp_gromos_NB_dict = dict(temp_gromos_NB_dict)
                temp_gromos_NB_dict['flag_add_params'] = True
            temp_gromos_NB_flags_dict = self.GROMOS_NB_Flags.prep_dict_4_init(temp_gromos_NB_dict)
            temp_gromos_NB_flags = self.GROMOS_NB_Flags(**temp_gromos_NB_flags_dict)
            temp_gr_nb_force = self._generate_openmm_gromos_nonbonded(temp_gromos_NB_flags, **temp_gromos_NB_dict)
            custom_nb_forces.append(temp_gr_nb_force)
        
        self.add_exclusions_pairs(nb_force, custom_nb_forces, **kwargs)
        
        # exclusion / pairs non-bonded interaction (GROMOS form) via custom bonded force
        custom_b_forces = []
        if gromos_B_flags is None:
            gromos_B_flags = []
        for temp_gromos_NB_dict in gromos_B_flags:# we are adding NB energy as bonds here (that's why NB naming is still used!)
            if not temp_gromos_NB_dict.get('flag_add_params'):
                temp_gromos_NB_dict = dict(temp_gromos_NB_dict)
                temp_gromos_NB_dict['flag_add_params'] = True
            temp_gromos_NB_flags_dict = self.GROMOS_NB_Flags.prep_dict_4_init(temp_gromos_NB_dict)
            temp_gromos_NB_flags = self.GROMOS_NB_Flags(**temp_gromos_NB_flags_dict)
            temp_gr_b_force = self._generate_openmm_gromos_excl_pairs(temp_gromos_NB_flags, **temp_gromos_NB_dict)
            custom_b_forces.append(temp_gr_b_force)
        
        # exclusion / pairs non-bonded interaction (GROMOS form) via custom bonded force

        # add NB forces to the system
        self.system.addForce(nb_force)
        for temp_nb_force in custom_nb_forces:
            self.system.addForce(temp_nb_force)
        for temp_b_force in custom_b_forces:
            self.system.addForce(temp_b_force)

        # set the pbc and cutoff
        self.set_pbc(pbc=pbc, nb_force=nb_force, custom_nb_forces=custom_nb_forces, custom_b_forces=custom_b_forces)

        ### bonded
        # add bonds / constraints
        if bond_constraints:
            self.add_constraints()
        else:
            self.add_bonds()    
        ### add angles; imp; dih
        self.add_angles()
        self.add_impropers()
        self.add_dihedrals()
        return
    
    def make_predefined_openmm_system(self, native_charges=False, N_water=0, bond_constraints=True, pbc=True, **kwargs):
        if native_charges:
            native_NB_flags = dict(use_charges=True)
            GROMOS_NB_flags = [dict(use_charges=False, use_lj=True)]
            # LJ 14
            GROMOS_B_flags = [dict(use_lj=True)]
            # 14 RF, RFC
            GROMOS_B_flags.append(dict(use_lj=False, use_charges=True, coul=False, use_Q_14=True, use_Q_excl=False, FG_name_suf='_14'))
        else:
            native_NB_flags = dict(use_charges=False, use_lj_water=False, use_lj_solute_solute=False)
            GROMOS_NB_flags = [dict(use_charges=True, use_lj=True)]
            # LJ 14
            GROMOS_B_flags = [dict(use_lj=True)]
            # 14 coul, RF, RFC
            GROMOS_B_flags.append(dict(use_lj=False, use_charges=True, use_Q_14=True, use_Q_excl=False, FG_name_suf='_14'))
        self.__add_excl_RF(GROMOS_B_flags)
        
        self.make_openmm_system(N_water=N_water, bond_constraints=bond_constraints, pbc=pbc, native_NB_flags=native_NB_flags,
                                gromos_NB_flags=GROMOS_NB_flags, gromos_B_flags=GROMOS_B_flags, **kwargs)
            
    def __add_excl_RF(self, GROMOS_B_flags):
        if self.ELE_params.flag_RF_excl:
            # excl RF, RFC
            GROMOS_B_flags.append(dict(use_lj=False, use_charges=True, coul=False, use_Q_14=False, use_Q_excl=True, FG_name_suf='_excl'))
            # self RFC
            temp_dict = dict(use_lj=False, use_charges=True, coul=False, RF=False, use_Q_14=False, use_Q_excl=False, FG_name_suf='_self', flag_self_term=True)
            GROMOS_B_flags.append(temp_dict)
        
    @staticmethod
    def _get_elem(atom):
        try:
            return app.element.Element.getByAtomicNumber(int(at.gr_get_element()))
        except:
            return

    def make_openmm_top(self, chains=None):
        top_atoms, top_atom_idx = list(self.top.get_atoms()), 0
        mm_top = app.Topology()
        mm_chain, chain_idx = mm_top.addChain(), 0
        for res in self.top.residues.values():
            if chains:
                if res not in chains[chain_idx]:
                    mm_chain, chain_idx = mm_top.addChain(), chain_idx + 1
                    assert res in chains[chain_idx]
            mm_res = mm_top.addResidue(res.name, mm_chain)
            for at in res.atoms:
                assert at is top_atoms[top_atom_idx] # ensures that order of atoms is the same as the order of residues x atoms_in_res
                top_atom_idx+=1
                mm_top.addAtom(at.name, self._get_elem(at), mm_res)
        for bi in self.top.bonds: # bi - bonded interaction
            atoms_idx = [at._openmm_id for at in bi.atoms]
            mm_top.addBond(*atoms_idx)
        return mm_top
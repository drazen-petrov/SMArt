#from .common_str_finder import find_common
from SMArt.alchemy.MCS import MCS
from SMArt.alchemy import top_matching_fnc

def _prep_tops(*tops):
    for t in tops:
        t.gen_graph()
        t.add_interactions2atoms_dict(replace = -1)
        try:
            t.find_imp_pairs()
        except:
            t.ff.imp_tetra_pair = {}
            t.ff.imp_flat_pair = {}

def _get_ptp_2_tops(top1, top2, **kwargs):
    _prep_tops(top1, top2)
    mcs = None
    mcs_kwargs = dict(kwargs)
    #temp_kwargs = dict(flag_top_update=True, flag_top_prune=kwargs.get('flag_top_prune', 'bond'))
    temp_kwargs = dict(flag_top_update=True)
    mcs_kwargs.update(temp_kwargs)
    enum_kwargs = dict(kwargs)
    #enum_kwargs['flag_score_fnc'] = kwargs.get('flag_score_fnc', 'bond')
    if kwargs.get('flag_get_res_common_atoms'):
        common_atoms, available_atoms_groups = top_matching_fnc.get_res_common_atoms(top1, top2)
        if available_atoms_groups:
            available_atoms_groups = available_atoms_groups[0] # this is only for 1 residue (first one that is not the same...)
        if common_atoms[0]:
            mcs = MCS(top1, top2, common_atoms = common_atoms, available_atoms = available_atoms_groups, **mcs_kwargs)
            s = mcs.initial_sol
            #mcs.make_estimates(s)
            #mcs.calc_score(s)
            #top_matching_fnc.generate_toptp(sol=s)
            mcs.enumerate_stepwise_sorted_call(s, **enum_kwargs)
    if mcs is None:
        mcs = MCS(top1, top2, **mcs_kwargs)
        if mcs.initial_sol._sol.shape[0]:
            mcs.enumerate_stepwise_sorted_call(mcs.initial_sol, **enum_kwargs)
        else:
            mcs.enumerate_stepwise_sorted(**enum_kwargs)
    return mcs

def point_mutation(top1, top2, **kwargs):
    """
    :param top1:
    :param top2:
    :param kwargs:
         flag_get_res_common_atoms - common atoms for residues with same names and len(atoms) - True as default
         flag_top_update - use topology information (True as default)
         flag_top_prune - which prune function to use (None as default; alternative 'bond' - makes sure no bonds are perturbed)
    :return:
        MCS instance with solutions (run generate_2state_top(mcs) to get the topology)
    """
    kwargs = dict(kwargs)
    kwargs['flag_get_res_common_atoms'] = kwargs.get('flag_get_res_common_atoms', True)
    return _get_ptp_2_tops(top1, top2, **kwargs)

def get_ptp_ligands(top1, top2, **kwargs):
    """
    :param top1:
    :param top2:
    :param kwargs:
         flag_get_res_common_atoms - common atoms for residues with same names and len(atoms) - False as default
         flag_top_update - use topology information (True as default)
         flag_top_prune - which prune function to use (None as default; alternative 'bond' - makes sure no bonds are perturbed)
    :return:
        MCS instance with solutions (run generate_2state_top(mcs) to get the topology)
    """
    kwargs = dict(kwargs)
    kwargs['flag_get_res_common_atoms'] = kwargs.get('flag_get_res_common_atoms', False)
    return _get_ptp_2_tops(top1, top2, **kwargs)

def generate_2state_top(mcs, solution=0, top_state=0, other_state=1, **kwargs):
    """
    generates a pairwise topology based on a given MCS solution
    :param mcs: MCS instance used to enumerate possible solutions
    :param solution: index of the solution (from mcs.solutions list)
    :param top_state: topology state (0 by default)
    :param other_state: perturbation state (1 by default)
    :param kwargs: kwargs to pass onto SMArt.alchemy.top_matching_fnc.generate_multi_state_top function
    :return:
    """
    sol = mcs.solutions[solution]
    gen_ms_top_kwargs = dict(kwargs)
    gen_ms_top_kwargs['top_state'] = top_state
    gen_ms_top_kwargs['other_state'] = other_state
    top_matching_fnc.generate_multi_state_top(sol, **gen_ms_top_kwargs)
    return sol


def _get_EDS_stepwise(mcs_kwargs, enum_kwargs, *tops, **kwargs):
    core_common_atoms = kwargs.get('core_common_atoms')
    if core_common_atoms:
        mcs = MCS(*tops[:2], common_atoms=core_common_atoms[:2], **mcs_kwargs)
    else:
        mcs = MCS(*tops[:2], **mcs_kwargs)
    mcs_steps = kwargs.get('mcs_steps', None)
    if mcs_steps is not None:mcs_steps.append(mcs)
    mcs.enumerate_stepwise_sorted(**enum_kwargs)
    for i in range(3, len(tops) + 1):
        if core_common_atoms:
            top_i = i-1
            next_top = tops[top_i]
            sol = mcs.solutions[0]
            df = sol._sol.copy()
            df[top_i] = None
            df.iloc[:len(core_common_atoms[top_i]),-1] = core_common_atoms[top_i] # update initial sol
            temp_mcs = MCS(*tops[:top_i+1], **mcs_kwargs, common_atoms=list(df.T.values))
            temp_mcs.enumerate_stepwise_sorted(**enum_kwargs)
        else:
            temp_mcs = MCS(*tops[:i], **mcs_kwargs)
            s = mcs.solutions[0].copy_add_state()
            s.tops = temp_mcs.tops
            s.available_atoms.append(temp_mcs.initial_sol.available_atoms[-1])
            temp_mcs.initial_sol = s
            temp_mcs.make_estimates(s)
            temp_mcs.calc_score(s)
            temp_mcs.enumerate_stepwise_sorted_call(s, **enum_kwargs)
        if mcs_steps is not None:mcs_steps.append(temp_mcs)
        mcs = temp_mcs
    return mcs

def _run_mcs_top0(tops):
    mcs_top0=[None]
    for i in range(1,len(tops)):
        mcs = get_EDS(tops[0], tops[i])
        mcs_top0.append(mcs)
    return mcs_top0

def _get_common_atms_mcs_top0(tops, mcs_top0):
    common_atoms_top0 = set(tops[0].get_atoms())
    for top_i in range(1,len(tops)):
        temp_common_set = set()
        temp_sol = mcs_top0[top_i].solutions[0]
        for sol_at in temp_sol.toptp.atoms.values():
            if None not in sol_at.ND_a_type_states:
                temp_common_set.add(temp_sol._sol.iloc[sol_at.id,0])
        common_atoms_top0 &= temp_common_set

    common_atoms = [list(common_atoms_top0)]
    for top_i in range(1,len(tops)):
        temp_common_atoms = []
        temp_sol = mcs_top0[top_i].solutions[0]
        for at_top0 in common_atoms[0]:
            row_idx = temp_sol.sol_atoms_map[(0, at_top0)]
            at_top_i = temp_sol._sol.iloc[row_idx, 1]
            temp_common_atoms.append(at_top_i)
        common_atoms.append(temp_common_atoms)
    return common_atoms

def get_common_atms_mcs_top0(tops):
    """
    runs MCS from tops[0] to all other tops[1:] andd finds all matches
    :param: tops - list of all tops
    :return:
        common_atoms - list of lists (already in the format to be used by MCS)
    """
    mcs_top0 = _run_mcs_top0(tops)
    return _get_common_atms_mcs_top0(tops, mcs_top0)

def get_EDS(*tops, **kwargs):
    """
    :param tops:
    :param kwargs:
        flag_stepwise - add topologies stepwise (True)
        flag_prune_EDS_match_mass (False)
        flag_get_res_common_atoms - uses residue and atom names to get common_atoms
        flag_get_core_common_atoms_top0 - uses get_common_atms_mcs_top0 to get core_common_atoms
        dihedral_match_v
        find_other_state
    :return:
        MCS instance with solutions (run generate_EDS_top(mcs) to get the EDS topology)
    """
    _prep_tops(*tops)
    mcs_kwargs = dict(kwargs)
    temp_mcs_kwargs = dict(flag_top_update = True, flag_top_prune = 'EDS', flag_score_fnc = 'EDS')
    mcs_kwargs.update(temp_mcs_kwargs)
    enum_kwargs = dict(kwargs)
    if 'flag_prune_EDS_match_mass' not in enum_kwargs:
       enum_kwargs['flag_prune_EDS_match_mass'] = False
    mcs = None
    if kwargs.get('flag_get_res_common_atoms'):
        common_atoms, available_atoms_groups = top_matching_fnc.get_res_common_atoms(*tops)
        if common_atoms[0]:
            mcs = MCS(*tops, common_atoms = common_atoms, available_atoms = available_atoms_groups[0], **mcs_kwargs)
            s = mcs.initial_sol
            mcs.make_estimates(s)
            mcs.calc_score(s)
            top_matching_fnc.generate_toptp(sol=s)
            mcs.enumerate_stepwise_sorted_call(s, **enum_kwargs)
    if mcs is None:
        core_common_atoms = kwargs.get('core_common_atoms')
        if core_common_atoms is None and kwargs.get('flag_get_core_common_atoms_top0'):
            core_common_atoms=get_common_atms_mcs_top0(tops)
            kwargs = dict(kwargs)
            kwargs['core_common_atoms'] = core_common_atoms
        if kwargs.get('flag_stepwise', True):
            mcs  = _get_EDS_stepwise(mcs_kwargs, enum_kwargs, *tops, **kwargs)
        else:
            if core_common_atoms:
                mcs = MCS(*tops, **mcs_kwargs, common_atoms=core_common_atoms)
            else:
                mcs = MCS(*tops, **mcs_kwargs)
            mcs.enumerate_stepwise_sorted(**enum_kwargs)
    return mcs

def generate_EDS_top(mcs, solution=0, flag_EDS_mass=True, find_other_state='any', **kwargs):
    gen_ms_top_kwargs = dict(kwargs)
    gen_ms_top_kwargs['flag_EDS_mass'] = flag_EDS_mass
    gen_ms_top_kwargs['find_other_state'] = find_other_state
    sol = mcs.solutions[solution]
    top_matching_fnc.generate_multi_state_top(sol, **gen_ms_top_kwargs)
    state_names = top_matching_fnc.get_state_names(sol)
    return sol, state_names

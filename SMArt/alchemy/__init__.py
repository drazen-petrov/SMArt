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
    temp_kwargs = dict(flag_top_update=True, flag_top_prune=kwargs.get('flag_top_prune', 'bond'))
    mcs_kwargs.update(temp_kwargs)
    enum_kwargs = dict(kwargs)
    enum_kwargs['flag_score_fnc'] = kwargs.get('flag_score_fnc', 'bond')
    if kwargs.get('flag_get_res_common_atoms'):
        common_atoms, available_atoms_groups = top_matching_fnc.get_res_common_atoms(top1, top2)
        if common_atoms[0]:
            mcs = MCS(top1, top2, common_atoms = common_atoms, available_atoms = available_atoms_groups[0], **mcs_kwargs)
            s = mcs.initial_sol
            mcs.make_estimates(s)
            mcs.calc_score(s)
            top_matching_fnc.generate_toptp(sol=s)
            mcs.enumerate_stepwise_sorted_call(s, **enum_kwargs)
    if mcs is None:
        mcs = MCS(top1, top2, **mcs_kwargs)
        mcs.enumerate_stepwise_sorted(**enum_kwargs)
    sol = mcs.solutions[0]
    gen_ms_top_kwargs = dict(kwargs)
    gen_ms_top_kwargs['other_state'] = 1
    top_matching_fnc.generate_multi_state_top(sol, **gen_ms_top_kwargs)
    return sol

def point_mutation(top1, top2, **kwargs):
    """
    :param top1:
    :param top2:
    :param kwargs:
         flag_get_res_common_atoms - common atoms for residues with same names and len(atoms) - True as default
         flag_top_update - use topology information (True as default)
         flag_top_prune - which prune function to use (bond as default - makes sure no bonds are perturbed)
    :return:
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
         flag_top_prune - which prune function to use (bond as default - makes sure no bonds are perturbed)
    :return:
    """
    kwargs = dict(kwargs)
    kwargs['flag_get_res_common_atoms'] = kwargs.get('flag_get_res_common_atoms', False)
    return _get_ptp_2_tops(top1, top2, **kwargs)

def _get_EDS_stepwise(mcs_kwargs, enum_kwargs, *tops):
    mcs = MCS(*tops[:2], **mcs_kwargs)
    mcs.enumerate_stepwise_sorted(**enum_kwargs)
    for i in range(3, len(tops) + 1):
        temp_mcs = MCS(*tops[:i], **mcs_kwargs)
        s = mcs.solutions[0].copy_add_state()
        s.tops = temp_mcs.tops
        s.available_atoms.append(temp_mcs.initial_sol.available_atoms[-1])
        temp_mcs.initial_sol = s
        temp_mcs.make_estimates(s)
        temp_mcs.calc_score(s)
        temp_mcs.enumerate_stepwise_sorted_call(s, **enum_kwargs)
        mcs = temp_mcs
    return mcs

def get_EDS(*tops, **kwargs):
    """
    :param tops:
    :param kwargs:
        flag_stepwise - add topologies stepwise (True)
        flag_prune_EDS_match_mass
        dihedral_match_v
        find_other_state
    :return:
    """
    _prep_tops(*tops)
    mcs_kwargs = dict(kwargs)
    temp_mcs_kwargs = dict(flag_top_update = True, flag_top_prune = 'EDS', flag_score_fnc = 'EDS')
    mcs_kwargs.update(temp_mcs_kwargs)
    enum_kwargs = dict(kwargs)
    if 'flag_prune_EDS_match_mass' not in enum_kwargs:
       enum_kwargs['flag_prune_EDS_match_mass'] = False
    if kwargs.get('flag_stepwise', True):
        mcs  = _get_EDS_stepwise(mcs_kwargs, enum_kwargs, *tops)
    else:
        mcs = MCS(tops, **mcs_kwargs)
        mcs.enumerate_stepwise_sorted(**enum_kwargs)
    gen_ms_top_kwargs = dict(kwargs)
    gen_ms_top_kwargs['flag_EDS_mass'] = True
    if 'find_other_state' not in gen_ms_top_kwargs:
        gen_ms_top_kwargs['find_other_state'] = 'any'
    sol = mcs.solutions[0]
    top_matching_fnc.generate_multi_state_top(sol, **gen_ms_top_kwargs)
    state_names = top_matching_fnc.get_state_names(sol)
    return sol, state_names


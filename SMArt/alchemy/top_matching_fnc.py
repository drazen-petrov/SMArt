from SMArt.incl import np, combinations, permutations, do_warn, deque, Counter
from SMArt.alchemy.incl import Dummy
from SMArt.md.incl import check_if_eq, DescriptionPart, Interaction, DihedralType, ImproperType, ExclusionPairType
from SMArt.md.data_st import Topology, MoleculeType


def get_res_common_atoms(*tops):
    N_res = len(tops[0].residues)
    flag_N_res = True
    for t in tops[1:]:
        if len(t.residues)!=N_res:
            flag_N_res = False
            break
    if flag_N_res:
        top_res = tuple(tuple(t.get_residues()) for t in tops)
        common_atoms = [[] for i in range(len(tops))]
        available_atoms_groups = []
        for res_i in range(N_res):
            temp_ND_res_states = {}
            for top_i in range(len(tops)):
                res_key = len(top_res[top_i][res_i].atoms), top_res[top_i][res_i].name # temp_N_at, res_name
                if res_key not in temp_ND_res_states:
                    temp_ND_res_states[res_key] = set()
                temp_ND_res_states[res_key].add(top_i)
            if len(temp_ND_res_states) < len(tops):
                for temp_ND_res_state in temp_ND_res_states:
                    for at_i in range(temp_ND_res_state[0]):
                        for top_i in range(len(tops)):
                            if top_i in temp_ND_res_states[temp_ND_res_state]:
                                common_atoms[top_i].append(top_res[top_i][res_i].atoms[at_i])
                            else:
                                common_atoms[top_i].append(None)
            else:
                #temp_available_atoms_group = [frozenset(top_res[top_i][res_i].atoms) for top_i in range(len(tops))]
                temp_available_atoms_group = [top_res[top_i][res_i].atoms for top_i in range(len(tops))]
                available_atoms_groups.append(temp_available_atoms_group)
        if not available_atoms_groups:
            available_atoms_groups.append([[] for top_i in range(len(tops))])
        return common_atoms, available_atoms_groups
    else:
        return None, None

def generate_toptp(N_tops = None, tops = None, sol = None, **kwargs):
    if N_tops is None:
        if tops is None:
            assert sol is not None
            tops = sol.tops
        N_tops = len(tops)
        klass = tops[0].__class__
    else:
        klass = kwargs.get('top_class', Topology)
        #klass = kwargs.get('top_class', MoleculeType)
    toptp = klass()
    toptp.get_container('excl_pair', create=True, db_type=dict)
    toptp.ptp_excl_pair = 0
    toptp.ptp_atom = np.zeros(3) # mass, p_charge, atom_type
    toptp.ptp_int = {} # for each interaction type (e.g. bonds), a counter of ptp matches
    toptp.ptp_make_break_int = {} #  for each interaction type (e.g. bonds), a counter of ptp make breaks
    toptp.interaction_states_map = tuple(dict() for _ in range(N_tops)) # for each top one dict;
        # in dict for each interaction map to its sol interaction
    if sol:
        sol.toptp = toptp
        # atoms
        atoms_in_sol = [list() for i in sol._sol.columns]
        for row_idx, row in sol._sol.iterrows():
            sol_at_i = None
            for top_i in sol._sol.columns:
                at_i = row[top_i]
                if at_i is not None:
                    if sol_at_i is None:
                        sol_at_i = add_new_toptp_atom(sol, row_idx)
                        ptp_atom = add_atom_to_row(sol_at_i, (top_i, at_i))
                    else:
                        ptp_atom = add_atom_to_row(sol_at_i, (top_i, at_i))
                        sol.toptp.ptp_atom += ptp_atom
                    if at_i != Dummy:
                        atoms_in_sol[top_i].append(row_idx)
        # exclusion vdw_pair pairs
        # 0 normal wdv; 1 excluded; 2 pair
        for row_idx, row in sol._sol.iterrows():
            for top_i in sol._sol.columns:
                #if row[top_i] not in (None, Dummy):
                     update_ptp_excl_pair_row_atom_pair(sol, [top_i, row[top_i]], sol.toptp.atoms[row_idx], **kwargs)
        """
        for top_i in sol._sol.columns:
            for row_idx in atoms_in_sol[top_i]:
                generate_ptp_excl_pair_single_atom(sol, row_idx, top_i, **kwargs)
        atoms_in_sol = [set(at_in_sol_top_i) for at_in_sol_top_i in atoms_in_sol]
        for top_pair in combinations(sol._sol.columns, 2):
            temp_atoms_in_sol = sorted(atoms_in_sol[top_pair[0]] & atoms_in_sol[top_pair[1]])
            for row_idx_pair in combinations(temp_atoms_in_sol, 2):
                sol_atom_pair = frozenset(sol.toptp.atoms[at_idx] for at_idx in row_idx_pair)
                if sol_atom_pair in sol.toptp.excl_pair:
                    sol_interaction = sol.toptp.excl_pair[sol_atom_pair]
                    for top_i in top_pair:
                        at_i_1, at_i_2 = tuple(sol._sol.loc[at_idx, top_i] for at_idx in row_idx_pair)
                        temp_state = get_exclusion_pair_state(at_i_1, at_i_2)
                        add_state2EP_ND_states(sol, sol_interaction, temp_state, top_i)
                #################################
                # only DEBUGGING ###############
                #################################
                else:
                    for top_i in top_pair:
                        at_i_1, at_i_2 = tuple(sol._sol.loc[at_idx, top_i] for at_idx in row_idx_pair)
                        temp_state = get_exclusion_pair_state(at_i_1, at_i_2)
                        assert temp_state.p == (0,)
                #################################
                # only DEBUGGING ###############
                #################################
        """
        # interactions
        for top_pos in range(1, len(sol.tops)):
            done_interactions = set()
            top2add_i = sol._sol.columns[top_pos] # this is actually the same
            for row_idx, row in sol._sol.iterrows():
                top_atom_list_01 = ([], [])
                if row[top2add_i] not in (None, Dummy):
                    top_atom_list_01[1].append((top2add_i, row[top2add_i]))
                    for top_i in sol._sol.columns[:top_pos]:
                        if row[top_i] not in (None, Dummy):
                            top_atom_list_01[0].append((top_i, row[top_i]))
                    update_ptp_interactions(sol, top_atom_list_01, done_interactions = done_interactions, **kwargs)
    return toptp

def get_generator_ptp_excl_pair_single_atom(sol, row_idx, top_i, done_sol_atom_pairs = None, **kwargs):
    # exclusion vdw_pair pairs
    # 0 normal wdv; 1 excluded; 2 pair
    if done_sol_atom_pairs is None:
        done_sol_atom_pairs = sol.toptp.excl_pair
    at_i = sol._sol.loc[row_idx, top_i]
    # from EP_l
    for other_at_i in sol.tops[top_i].EP_l[at_i]:
        if (top_i, other_at_i) in sol.sol_atoms_map:
            other_row_idx = sol.sol_atoms_map[(top_i, other_at_i)]
            sol_atom_pair = frozenset((sol.toptp.atoms[row_idx], sol.toptp.atoms[other_row_idx]))
            if sol_atom_pair not in done_sol_atom_pairs:
                ep_state = get_exclusion_pair_state_from_EP(at_i, other_at_i, sol.tops[top_i])
                yield sol_atom_pair, ep_state

def generate_ptp_excl_pair_single_atom(sol, row_idx, top_i, done_sol_atom_pairs = None, **kwargs):
    for sol_atom_pair, temp_state in get_generator_ptp_excl_pair_single_atom(sol, row_idx, top_i, **kwargs):
        add_new_excl_pair(sol, sol_atom_pair, temp_state, top_i)

def add_new_excl_pair(sol, sol_atom_pair, temp_state, top_i):
    sol_interaction = sol.toptp.Interaction(ExclusionPairType)
    sol_interaction.ND_states = {temp_state: {top_i}}
    sol.toptp.add2container(sol_interaction, item_id=sol_atom_pair, create=True)
    sol_at1, sol_at2 = sol_atom_pair
    sol.toptp.add_atom_pair2EP_l(sol_at1, sol_at2)
    return sol_interaction

def update_ptp(sol, top_atom_pair, old_sol, **kwargs):
    # get list of atoms that are being matched:
    # (1,1) for atom atom; (n,1) for row atom; (n,m) for row row
    top_atom_list_01 = ([], [])
    if top_atom_pair[1][1] != Dummy:
        if top_atom_pair[0][0] is None:
            row_idx0 = top_atom_pair[0][1]
            for top_idx_at_pair_0 in enumerate(old_sol._sol.loc[row_idx0]):
                if top_idx_at_pair_0[1] not in (None, Dummy):
                    top_atom_list_01[0].append(top_idx_at_pair_0)
            if top_atom_pair[1][0] is None:# both rows
                row_idx1 = top_atom_pair[1][1]
                for top_idx_at_pair_1 in enumerate(old_sol._sol.loc[row_idx1]):
                    if top_idx_at_pair_1[1] not in (None, Dummy):
                        top_atom_list_01[1].append(top_idx_at_pair_1)
            else:# row atom
                top_atom_list_01[1].append(top_atom_pair[1])
        else:#atom atom new row
            top_atom_list_01[0].append(top_atom_pair[0])
            top_atom_list_01[1].append(top_atom_pair[1])
    ########################
    # atom ptp and excl pair
    ########################
    if top_atom_pair[0][0] is None:
        row_idx = top_atom_pair[0][1]
        temp_at = sol.toptp.atoms[row_idx]
        if top_atom_pair[1][0] is None:# both rows
            row_idx2 = top_atom_pair[1][1]
            temp_at2 = sol.toptp.atoms.pop(row_idx2)
            for m in temp_at2.ND_m_states:
                add_atom_m_to_row(temp_at, m, temp_at2.ND_m_states[m], sol.toptp.ptp_atom)
            for p_ch in temp_at2.ND_pch_states:
                add_atom_pch_to_row(temp_at, p_ch, temp_at2.ND_pch_states[p_ch], sol.toptp.ptp_atom)
            for a_type in temp_at2.ND_a_type_states:
                add_atom_a_type_to_row(temp_at, a_type, temp_at2.ND_a_type_states[a_type], sol.toptp.ptp_atom)
            ########################################################
            # pair excl
            ########################################################
            update_ptp_excl_pair_row_row_pair(sol, (temp_at, temp_at2), top_atom_list_01, **kwargs)
            #del(sol.toptp.atoms[row_idx2])
        else: # atom 2 row
            sol.toptp.ptp_atom += add_atom_to_row(temp_at, top_atom_pair[1])
            update_ptp_excl_pair_row_atom_pair(sol, top_atom_pair[1], temp_at, **kwargs)
    else: # atom atom pair - new row
        row_idx = sol.sol_atoms_map[top_atom_pair[0]]
        temp_at = add_new_toptp_atom(sol, row_idx)
        for i in range(2):
            ptp_atom = add_atom_to_row(temp_at, top_atom_pair[i])
        sol.toptp.ptp_atom += ptp_atom
        update_ptp_excl_pair_atom_atom_pair(sol, row_idx, top_atom_pair, **kwargs)
    ################################
    # interactions #################
    ################################
    if top_atom_pair[1][1] == Dummy:
        return
    update_ptp_interactions(sol, top_atom_list_01, **kwargs)

"""
def get_exclusion_pair_state(at1, at2):
    if at2 in at1.e_l:
        assert at1 in at2.e_l
        return ExclusionType(params = [1])
    elif at2 in at1.p_l:
        assert at1 in at2.p_l
        return ExclusionType(params = [2])
    else:
        assert at1 not in at2.e_l and at1 not in at2.p_l
        return ExclusionType(params = [0])
"""

def get_exclusion_pair_state_from_EP(at1, at2, top, top_state = 0):
    if at2 in top.EP_l[at1]:
        assert at1 in top.EP_l[at2]
        return top.excl_pair[frozenset((at1, at2))].states[top_state]
    else:
        assert at1 not in top.EP_l[at2]
        return ExclusionPairType(params = (0,))

def add_state2EP_ND_states(sol, sol_interaction, temp_state, top_idx):
    match_ND_state = None
    for ND_state in sol_interaction.ND_states:
        if ND_state.p == temp_state.p:
            match_ND_state = ND_state
            break
    if match_ND_state is None:
        sol.toptp.ptp_excl_pair += 1
        sol_interaction.ND_states[temp_state] = {top_idx}
    else:
        sol_interaction.ND_states[match_ND_state].add(top_idx)

def merge_EP_ND_states(sol, sol_interaction, sol_interaction2):
    done_ND_states2 = set()
    for ND_state in sol_interaction.ND_states:
        flag_not_matched = True
        for temp_state in sol_interaction2.ND_states:
            if ND_state.p == temp_state.p:
                flag_not_matched = False
                done_ND_states2.add(temp_state)
                sol_interaction.ND_states[ND_state] |= sol_interaction2.ND_states[temp_state]
                break
        if flag_not_matched:
            sol.toptp.ptp_excl_pair += 1
    for ND_state in sol_interaction2.ND_states:
        if ND_state not in done_ND_states2:
            sol_interaction.ND_states[ND_state] = sol_interaction2.ND_states[ND_state]
            sol.toptp.ptp_excl_pair += 1

def update_ptp_excl_pair_atom_atom_pair(sol, row_idx, top_atom_pair, **kwargs):
    for i in range(2):
        #if top_atom_pair[i][1] not in (None, Dummy):
            update_ptp_excl_pair_row_atom_pair(sol, top_atom_pair[i], sol.toptp.atoms[row_idx], **kwargs)

def update_ptp_excl_pair_row_atom_pair(sol, top_atom, sol_atom, **kwargs):
    if top_atom[1] in (None, Dummy):
        return
    done_sol_atom_pair = set()
    for other_sol_at in sol.toptp.EP_l[sol_atom]:
    #for other_sol_at in sol_atom.EP_l:
        sol_atom_pair = frozenset((sol_atom, other_sol_at))
        sol_interaction = sol.toptp.excl_pair[sol_atom_pair]
        other_atom = sol._sol.loc[other_sol_at.id, top_atom[0]]
        if other_atom not in (None, Dummy):
            temp_state = get_exclusion_pair_state_from_EP(top_atom[1], other_atom, sol.tops[top_atom[0]])
            add_state2EP_ND_states(sol, sol_interaction, temp_state, top_atom[0])
        done_sol_atom_pair.add(sol_atom_pair)
    temp_gen = get_generator_ptp_excl_pair_single_atom(sol, sol_atom.id, top_atom[0], done_sol_atom_pair, **kwargs)
    for sol_atom_pair, sol_state in temp_gen:
        sol_interaction = add_new_excl_pair(sol, sol_atom_pair, sol_state, top_atom[0])
        for top_i in sol._sol.columns:
            if top_i!=top_atom[0]: # skip the top from which the added atom comes from
                row_idx_pair = [sol_at.id for sol_at in sol_atom_pair]
                atom_pair_i = sol._sol.loc[row_idx_pair, top_i].values
                if None not in atom_pair_i and Dummy not in atom_pair_i:
                    temp_state = get_exclusion_pair_state_from_EP(atom_pair_i[0], atom_pair_i[1], sol.tops[top_i])
                    add_state2EP_ND_states(sol, sol_interaction, temp_state, top_i)

def update_ptp_excl_pair_row_row_pair(sol, sol_atoms_pair, top_atom_list_01, **kwargs):
    ########################
    # at this point, in sol, we just poped out sol_atoms_pair[1] from toptp.atoms
    # so we need to be sure that searching with the id is not used...
    ########################
    done_sol_atom_pair = set()
    sol_atom_0, sol_atom_1 = sol_atoms_pair
    for other_sol_at in sol.toptp.EP_l[sol_atom_0]:
    #for other_sol_at in sol_atom_0.EP_l:
        sol_atom_pair_0 = frozenset((sol_atom_0, other_sol_at))
        sol_interaction_0 = sol.toptp.excl_pair[sol_atom_pair_0]
        sol_atom_pair_1 = frozenset((sol_atom_1, other_sol_at))
        if sol_atom_pair_1 in sol.toptp.excl_pair:
            sol_interaction_1 = sol.toptp.excl_pair.pop(sol_atom_pair_1)
            merge_EP_ND_states(sol, sol_interaction_0, sol_interaction_1)
        else:
            for top_atom_1 in top_atom_list_01[1]:
                other_atom_1 = sol._sol.loc[other_sol_at.id, top_atom_1[0]]
                if other_atom_1 not in (None, Dummy):
                    temp_state = get_exclusion_pair_state_from_EP(top_atom_1[1], other_atom_1, sol.tops[top_atom_1[0]])
                    add_state2EP_ND_states(sol, sol_interaction_0, temp_state, top_atom_1[0])
                    assert temp_state.p == (0,)
        done_sol_atom_pair.add(sol_atom_pair_1)
    EP_l_sol_atom_1 = sol.toptp.EP_l.pop(sol_atom_1)
    for other_sol_at in EP_l_sol_atom_1:
        sol_atom_pair_1 = frozenset((sol_atom_1, other_sol_at))
        if sol_atom_pair_1 in done_sol_atom_pair:continue
        sol_interaction_1 = sol.toptp.excl_pair[sol_atom_pair_1]
        sol_atom_pair_0 = frozenset((sol_atom_0, other_sol_at))
        assert sol_atom_pair_0 not in sol.toptp.excl_pair
        sol.toptp.excl_pair[sol_atom_pair_0] = sol_interaction_1
        for top_atom_0 in top_atom_list_01[0]:
            other_atom_0 = sol._sol.loc[other_sol_at.id, top_atom_0[0]]
            if other_atom_0 not in (None, Dummy):
                temp_state = get_exclusion_pair_state_from_EP(top_atom_0[1], other_atom_0, sol.tops[top_atom_0[0]])
                add_state2EP_ND_states(sol, sol_interaction_1, temp_state, top_atom_0[0])
                assert temp_state.p == (0,)
    for sol_at in sol.toptp.get_atoms():
        assert sol_at != sol_atom_1
#        if sol_at != sol_atom_1:
        if sol_atom_1 in sol.toptp.EP_l[sol_at]:
            sol.toptp.EP_l[sol_at].remove(sol_atom_1)
            assert sol_at!=sol_atom_0
            sol.toptp.add_atom_pair2EP_l(sol_at, sol_atom_0)

def add_atom_to_row(sol_at, top_atom_id):
    st_at, at = top_atom_id # st_at is the same top_idx
    st_at = {st_at}
    ptp_atom = np.zeros(3)
    add_atom_m_to_row(sol_at, at.m, st_at, ptp_atom)
    add_atom_pch_to_row(sol_at, at.p_ch, st_at, ptp_atom)
    add_atom_a_type_to_row(sol_at, at.a_type, st_at, ptp_atom)
    return ptp_atom

def update_ptp_interactions(sol, top_atom_list_01, done_interactions = None, **kwargs):
    if done_interactions is None:
        done_interactions = set()
    for top_at0 in top_atom_list_01[0]:
        for top_at1 in top_atom_list_01[1]:
            top_pair_ind = (top_at0[0], top_at1[0])
            top0 = sol.tops[top_pair_ind[0]]
            top1 = sol.tops[top_pair_ind[1]]
            top_01 = (top0, top1)
            at0 = top_at0[1]
            at1 = top_at1[1]
            for int_matches, int_ptp, container_name in find_interactions_ptp(at0, at1, top_01, top_pair_ind, sol,
                                                                              done_interactions, **kwargs):
                #############################
                # same (matched) interactions
                #############################
                for int_match in int_matches:
                    other_int_states = []
                    other_int_state = None
                    if int_match[0] in sol.toptp.interaction_states_map[top_pair_ind[0]]:
                        sol_interaction = sol.toptp.interaction_states_map[top_pair_ind[0]][int_match[0]]
                        sol_ND_state_tops_set = find_ND_state_tops_set(sol_interaction.ND_states, int_match[0].states[0], top_pair_ind[0])
                        if int_match[1] in sol.toptp.interaction_states_map[top_pair_ind[1]]:
                            other_sol_interaction = sol.toptp.interaction_states_map[top_pair_ind[1]][int_match[1]]
                            other_int_states = other_sol_interaction.int_states
                        else:
                            other_int_state = top_pair_ind[1], int_match[1]
                    else:
                        if int_match[1] in sol.toptp.interaction_states_map[top_pair_ind[1]]:
                            sol_interaction = sol.toptp.interaction_states_map[top_pair_ind[1]][int_match[1]]
                            sol_ND_state_tops_set = find_ND_state_tops_set(sol_interaction.ND_states,
                                                                           int_match[1].states[0], top_pair_ind[1])
                            other_int_state = top_pair_ind[0], int_match[0]
                        else:
                            new_interaction = top0.Interaction(int_match[0].int_type)
                            new_interaction.ND_states = {int_match[0].states[0]: set(top_pair_ind)}
                            new_interaction.int_states = [None] * len(sol.tops)
                            for i in range(2):
                                new_interaction.int_states[top_pair_ind[i]] = int_match[i]
                                sol.toptp.interaction_states_map[top_pair_ind[i]][int_match[i]] = new_interaction
                            #sol.toptp.add2container(new_interaction, create=True, db_type=list)
                            done_interactions.add(new_interaction)
                    if other_int_states:
                        for top_i_other, int_match_other in enumerate(other_int_states):
                            sol_interaction.int_states[top_i_other] = int_match_other
                            sol_ND_state_tops_set.add(top_i_other)
                            sol.toptp.interaction_states_map[top_i_other][int_match_other] = sol_interaction
                        done_interactions.add(sol_interaction)
                    if other_int_state:
                        top_i_other, int_match_other = other_int_state
                        sol_interaction.int_states[top_i_other] = int_match_other
                        sol_ND_state_tops_set.add(top_i_other)
                        sol.toptp.interaction_states_map[top_i_other][int_match_other] = sol_interaction
                        done_interactions.add(sol_interaction)
                ######################################
                # different (not matched) interactions
                ######################################
                for ii in range(len(int_ptp[0])):
                    sol_interaction = None
                    sol_interaction_2 = None
                    for i in range(2):
                        if int_ptp[i][ii] in sol.toptp.interaction_states_map[top_pair_ind[i]]:
                            if sol_interaction is None:
                                sol_interaction = sol.toptp.interaction_states_map[top_pair_ind[i]][int_ptp[i][ii]]
                                top_i_other = (i + 1) % 2
                            else:
                                sol_interaction_2 = sol.toptp.interaction_states_map[top_pair_ind[i]][int_ptp[i][ii]]
                    if sol_interaction_2 and sol_interaction is not sol_interaction_2:
                        sol_interaction.ND_states.update(sol_interaction_2.ND_states)
                        for top_set in sol_interaction_2.ND_states.values():
                            for top_i in top_set:
                                assert sol_interaction.int_states[top_i] is None
                                int_i = sol_interaction_2
                                sol_interaction.int_states[top_i] = int_i
                                sol.toptp.interaction_states_map[top_i][int_i] = sol_interaction
                        done_interactions.add(sol_interaction_2)
                    else:
                        if sol_interaction is None:
                            for i in range(2):
                                if int_ptp[i][ii] is not None and sol_interaction is None:
                                    sol_interaction = top_01[i].Interaction(int_ptp[i][ii].int_type)
                                    sol_interaction.ND_states = {int_ptp[i][ii].states[0]: {top_pair_ind[i]}}
                                    sol_interaction.int_states = [None] * len(sol.tops)
                                    sol_interaction.int_states[top_pair_ind[i]] = int_ptp[i][ii]
                                    sol.toptp.interaction_states_map[top_pair_ind[i]][int_ptp[i][ii]] = sol_interaction
                                    top_i_other = (i + 1) % 2
                        if int_ptp[top_i_other][ii] is not None:
                            i = top_i_other
                            sol_interaction.ND_states[int_ptp[i][ii].states[0]] = {top_pair_ind[i]}
                            sol_interaction.int_states[top_pair_ind[i]] = int_ptp[i][ii]
                            sol.toptp.interaction_states_map[top_pair_ind[i]][int_ptp[i][ii]] = sol_interaction
                        else:
                            if container_name not in sol.toptp.ptp_make_break_int:
                                sol.toptp.ptp_make_break_int[container_name] = 0
                            sol.toptp.ptp_make_break_int[container_name] += 1
                    done_interactions.add(sol_interaction)
                    if container_name not in sol.toptp.ptp_int:
                        sol.toptp.ptp_int[container_name] = 0
                    sol.toptp.ptp_int[container_name] += 1

def find_ND_state_tops_set(ND_states, state_top_i, top_i):
    if state_top_i in ND_states and top_i in ND_states[state_top_i]:
        return ND_states[state_top_i]
    for temp_state, tops_set in ND_states.items():
        if top_i in tops_set:
            return tops_set

def add_new_toptp_atom(sol, row_idx):
    sol_at = sol.toptp.add_atom(row_idx)
    #sol_at.sol_id = sol_at.id
    sol_at.ND_m_states = {}
    sol_at.ND_pch_states = {}
    sol_at.ND_a_type_states = {}
    #sol_at.EP_l = set()
    return sol_at

def add_atom_m_to_row(sol_at, at_m, st_at, ptp_atom):
    if at_m is None:
        if at_m not in sol_at.ND_m_states:
            sol_at.ND_m_states[at_m] = set()
            #ptp_atom[0] += 1
        sol_at.ND_m_states[at_m] |= st_at
    else:
        flag = True
        for m in sol_at.ND_m_states:
            if check_if_eq(m, at_m):
                sol_at.ND_m_states[m] |= st_at
                flag = False
                break
        if flag:
            sol_at.ND_m_states[at_m] = st_at
            ptp_atom[0] += flag

def add_atom_pch_to_row(sol_at, at_p_ch, st_at, ptp_atom):
    if at_p_ch is None:
        at_p_ch = 0.
        """
        if at_p_ch not in sol_at.ND_pch_states:
            sol_at.ND_pch_states[at_p_ch] = set()
            ptp_atom[1] += 1
        sol_at.ND_pch_states[at_p_ch] |= st_at
        """
    if at_p_ch is None:
        """
        if at_p_ch not in sol_at.ND_pch_states:
            sol_at.ND_pch_states[at_p_ch] = set()
            ptp_atom[1] += 1
        sol_at.ND_pch_states[at_p_ch] |= st_at
        """
        pass
    else:
        flag = True
        for pch in sol_at.ND_pch_states:
            if check_if_eq(pch, at_p_ch):
                sol_at.ND_pch_states[pch] |= st_at
                flag = False
                break
        if flag:
            sol_at.ND_pch_states[at_p_ch] = st_at
            ptp_atom[1] += flag

def add_atom_a_type_to_row(sol_at, at_a_type, st_at, ptp_atom):
    if at_a_type is None:
        if at_a_type not in sol_at.ND_a_type_states:
            sol_at.ND_a_type_states[at_a_type] = set()
            ptp_atom[2] += 1
        sol_at.ND_a_type_states[at_a_type] |= st_at
    else:
        flag = True
        for a_type in sol_at.ND_a_type_states:
            if a_type and a_type.id == at_a_type.id:
                sol_at.ND_a_type_states[a_type] |= st_at
                flag = False
                break
        if flag:
            sol_at.ND_a_type_states[at_a_type] = st_at
            ptp_atom[2] += flag

# interactions
def _sep_int_groups_ptp_flags(int_list):
    if len(int_list)==0:
        return []
    sep_int_lists = [[int_list.pop()]]
    while int_list:
        temp_int = int_list.pop()
        flag_no_match = True
        for temp_int_list in sep_int_lists:
            temp_int2 = temp_int_list[0]
            if temp_int.states[0].check_eq_ptp_flag_params(temp_int2.states[0].p):
                temp_int_list.append(temp_int)
                flag_no_match = False
                break
        if flag_no_match:
            sep_int_lists.append([temp_int])
    return sep_int_lists

def _check_imp_tetrahedral_order(at0_list_mapped_to_1, at1_list):  # for tetrahedral
    for i in range(3):
        if at0_list_mapped_to_1[0] == at1_list[i]:
            for j in range(1, 3):
                if at0_list_mapped_to_1[j] != at1_list[(i + j) % 3]:
                    return False
            return True

def _convert_improper_params_atom_order(top_01, improper_0, improper_1, atoms_01_map, **kwargs):
    orig_params =  improper_1.states[0]
    flag_same_order = True
    for i, at_0 in enumerate(improper_0.atoms):
        if improper_1.atoms[i] != atoms_01_map[0][at_0]:
            flag_same_order = False
            break
    if flag_same_order:
        return orig_params
    bond_c_0 = top_01[0].improper_type(improper_0)
    bond_c_1 = top_01[1].improper_type(improper_1)
    if bond_c_0 != bond_c_1:
        do_warn('different improper types with different atom order!!!')
        return None
    if max(bond_c_0)==3:
        c_pos_0 = bond_c_0.index(3)
        c_pos_1 = bond_c_1.index(3)
        if improper_0.states[0] in top_01[0].ff.imp_flat_pair and improper_1.states[0] in top_01[1].ff.imp_flat_pair:
            imp_flat_c_pos_map = {0:0, 1:1, 2:1, 3:1}
            if imp_flat_c_pos_map[c_pos_0] == imp_flat_c_pos_map[c_pos_1]:
                return orig_params
            else:
                return top_01[1].ff.imp_tetra_pair[orig_params]
        elif c_pos_0 in (0,3) and c_pos_1 in (0,3):
            if improper_0.states[0] in top_01[0].ff.imp_tetra_pair and improper_1.states[0] in top_01[1].ff.imp_tetra_pair:
                ################################### this needs one more check!!!!!!!!!!!!!!!!!
                ################################### calc energies / angles and compare!!!!!!!!!!!!
                if c_pos_0 == 0:
                    at0_rest_list = improper_0.atoms[1:]
                else:
                    at0_rest_list = improper_0.atoms[2::-1]
                if c_pos_1 == 0:
                    at1_rest_list = improper_1.atoms[1:]
                else:
                    at1_rest_list = improper_1.atoms[2::-1]
                at0_rest_mapped_to_1 = [atoms_01_map[0][at_0] for at_0 in at0_rest_list]
                if _check_imp_tetrahedral_order(at0_rest_mapped_to_1, at1_rest_list):
                    return orig_params
                else:
                    return top_01[1].ff.imp_tetra_pair[orig_params]
        else:
            do_warn('one of the imporpers is not in flat pair and has the central atom in the middle!!!')
            return None
    else:
        flag_rev_same_order = True
        for i, at_0 in enumerate(improper_0.atoms):
            if improper_1.atoms[-(i+1)] != atoms_01_map[0][at_0]:
                flag_rev_same_order = False
                break
        if flag_rev_same_order:
            return orig_params
        else:
            do_warn('1-2-3-4 improper types with different atom order!!!')
            return None

def _find_interactions_ptp_atom_set(top_01, int_list0, int_list1, atoms_01_map, ptp_score_kwargs = None, **kwargs):
    if ptp_score_kwargs is None:
        ptp_score_kwargs = {}
    int_matches = []
    int_ptp = ([], [])
    if kwargs.get('use_ptp_flags'):
        # separated into sub-groups of the interactions that are allowed to be matches (e.g. dihedrals with same mult)
        sep_int_lists0 = _sep_int_groups_ptp_flags(list(int_list0))
        sep_int_lists1 = _sep_int_groups_ptp_flags(list(int_list1))
        sep_int_lists_01 = []
        while sep_int_lists0:
            temp_sep_int_list0 = sep_int_lists0.pop()
            flag_no_match = True
            for i1, temp_sep_int_list1 in enumerate(sep_int_lists1):
                if temp_sep_int_list0[0].states[0].check_eq_ptp_flag_params(temp_sep_int_list1[0].states[0].p):
                    sep_int_lists_01.append((temp_sep_int_list0, temp_sep_int_list1))
                    flag_no_match = False
                    sep_int_lists1.pop(i1)
                    break
            if flag_no_match:
                sep_int_lists_01.append((temp_sep_int_list0, []))
        while sep_int_lists1:
            sep_int_lists_01.append(([], sep_int_lists1.pop()))
        ptp_score_kwargs['flag_ptp_order'] = True
        for sep_int_list_01 in sep_int_lists_01:
            temp_int_matches, temp_int_ptp = _find_interactions_ptp_atom_set(top_01, sep_int_list_01[0],
                                                    sep_int_list_01[1], atoms_01_map, ptp_score_kwargs=ptp_score_kwargs)
            int_matches.extend(temp_int_matches)
            int_ptp[0].extend(temp_int_ptp[0])
            int_ptp[1].extend(temp_int_ptp[1])
        return int_matches, int_ptp
# find eq pairs of int
    int_list_ptp0 = []
    flag_improper = True
    if int_list0 and int_list0[0].int_type != ImproperType:
        flag_improper = False
    elif int_list1 and int_list1[0].int_type != ImproperType:
        flag_improper = False
    while int_list0:
        int0 = int_list0.pop()
        flag_no_match = True
        for i1, int1 in enumerate(int_list1):
            if flag_improper:
                params2check = _convert_improper_params_atom_order(top_01, int0, int1, atoms_01_map, **kwargs)
            else:
                params2check = int1.states[0]
            if params2check and int0.states[0].check_eq_params(params2check.p):
                int_matches.append((int0, int1))
                del int_list1[i1]
                flag_no_match = False
                break
        if flag_no_match:
            int_list_ptp0.append(int0)
# find best ptp combination of matches
    int_list0 = int_list_ptp0
    int_ptp = find_best_ptp_matches(int_list0, int_list1, flag_improper, top_01, atoms_01_map, ptp_score_kwargs, **kwargs)
    return int_matches, int_ptp


def find_best_ptp_matches(int_list0, int_list1, flag_improper, top_01, atoms_01_map, ptp_score_kwargs, **kwargs):
    int_ptp = ([], [])
    if not kwargs.get('flag_best_ptp_matches'):
        len_int_list_01 = len(int_list0), len(int_list1)
        N_max = max(len_int_list_01)
        for i in range(2):
            for ii in range(N_max):
                if ii < len_int_list_01[i]:
                    int_ptp[i].append((int_list0, int_list1)[i][ii])
                else:
                    int_ptp[i].append(None)
        return int_ptp
    if len(int_list0) > len(int_list1):
        fac = 1
    else:
        fac = -1
    int_list_01 = (int_list0, int_list1)[::fac]
    if int_list_01[1]:
        best_score = None
        for perm_int_list_longer in permutations(int_list_01[0], len(int_list_01[1])):
            score = None
            temp_int_matches = list(zip(perm_int_list_longer, int_list_01[1]))
            for temp_int_match in temp_int_matches:
                int0, int1 = temp_int_match[::fac]
                if flag_improper:
                    state2check = _convert_improper_params_atom_order(top_01, int0, int1, atoms_01_map, **kwargs)
                    temp_score = int0.states[0].calc_ptp_score(state2check, **ptp_score_kwargs)
                else:
                    temp_score = int0.states[0].calc_ptp_score(int1.states[0], **ptp_score_kwargs)
                if score is None:
                    score = temp_score
                else:
                    for i, param_sc in enumerate(temp_score):
                        score[i] += param_sc
            if best_score:
                if best_score > score:
                    best_perm_match = temp_int_matches
                    best_score = score
            else:
                best_perm_match = temp_int_matches
                best_score = score
        for temp_int_match in best_perm_match:
            temp_int_match = temp_int_match[::fac]
            int_ptp[0].append(temp_int_match[0])
            int_ptp[1].append(temp_int_match[1])
            int_list_01[0].remove(temp_int_match[0])
    for temp_int in int_list_01[0]:
        temp_int_match = (temp_int, None)[::fac]
        int_ptp[0].append(temp_int_match[0])
        int_ptp[1].append(temp_int_match[1])
    return int_ptp

def _check_extra_dih_in_sol(temp_int_i, top_i, top_i_idx, top_ii_idx, sol):
    for i in range(2):
        temp_count = 0
        dih_at = temp_int_i.atoms[1:3][i]
        dih_at_other = temp_int_i.atoms[1:3][(i +1)%2]
        for at_neigh_i in top_i.adj[dih_at]:
            if at_neigh_i!=dih_at_other:
                if (top_i_idx, at_neigh_i) in sol.sol_atoms_map:
                    at_neigh_ii = sol.find_common_atom_pair((top_i_idx, at_neigh_i), top_ii_idx)
                    if at_neigh_ii not in (Dummy, None):
                        temp_count += 1
                        #if temp_count > 1:break
                        break
        if temp_count!=1:
            return False
    return True

def _find_int_list_01_to_compare(sol, top_01, top_pair_ind, int_container_01, done_interactions, flag_dih_params, **kwargs):
    flag_int_allow_not_found = kwargs.get('flag_int_allow_not_found', True)
    if int_container_01[0] or int_container_01[1]:
        flag_dih, atom_range, dihedral_match_v = flag_dih_params
        for i in range(2):
            fac = (1, -1)[i]
            ii = (i+1) % 2
            for temp_int_i in int_container_01[i]:
                if temp_int_i in sol.toptp.interaction_states_map[top_pair_ind[i]]:
                    temp_sol_interaction = sol.toptp.interaction_states_map[top_pair_ind[i]][temp_int_i]
                    if temp_sol_interaction in done_interactions:
                        continue
                """                
                if temp_int_i in done_int_01[i]:
                    continue
                """
                atoms_i = temp_int_i.atoms[atom_range[0]:atom_range[1]]
                atoms_ii = []
                for int_at_i in atoms_i:
                    if (top_pair_ind[i], int_at_i) in sol.sol_atoms_map:
                        int_at_ii = sol.find_common_atom_pair((top_pair_ind[i], int_at_i), top_pair_ind[ii])
                        if int_at_ii not in (Dummy, None):
                            atoms_ii.append(int_at_ii)
                        else:
                            break
                    else:
                        break
                if len(atoms_i) == len(atoms_ii):
                    if flag_dih:# flag_dih is only true when it's dihedral and dihedral_match_version == 1
                        if not _check_extra_dih_in_sol(temp_int_i, top_01[i], top_pair_ind[i], top_pair_ind[ii], sol):
                            continue
                    int_list_i = top_01[i].find_interactions(atoms_i, int_container_01[i], v=dihedral_match_v,
                                                allow_not_found=flag_int_allow_not_found, allow_multiple_matches = True)
                    int_list_ii = top_01[ii].find_interactions(atoms_ii, int_container_01[ii], v=dihedral_match_v,
                                                allow_not_found=flag_int_allow_not_found, allow_multiple_matches = True)
                    yield (int_list_i, int_list_ii)[::fac], (atoms_i, atoms_ii)[::fac]

def gen_int_container_01_list_new(at0, at1, top_01, top_pair_ind, sol, **kwargs):
    """generator for interaction containers to be checked for matches"""
# all but dihedrals
    done_containers = []
    flag_dih_params = False, (None, None), 2
    #for int_container0, container_name in at0.get_interaction_containers(get_container_name=True):
    for int_container0, container_name in top_01[0].atom_interactions[at0].get_interaction_containers(get_container_name=True):
        done_int_01 = [], []
        if kwargs.get('dihedral_match_v', 1) == 1:
            #if int_container0 and int_container0[0].int_type==top_01[0].DihedralType:continue
            if int_container0 and int_container0[0].int_type == DihedralType:continue
        int_container1 = top_01[1].atom_interactions[at1].get_container(container_name, allow_not_found = True)
        if int_container1 is None:
            int_container1 = []
        done_containers.append(container_name)
        yield (int_container0, int_container1), done_int_01, flag_dih_params, container_name
    #for int_container1, container_name in at1.get_interaction_containers(get_container_name=True):
    for int_container1, container_name in top_01[1].atom_interactions[at1].get_interaction_containers(get_container_name=True):
        if container_name in done_containers:
            continue
        done_int_01 = [], []
        if kwargs.get('dihedral_match_v', 1) == 1:
            #if int_container1 and int_container1[0].int_type==top_01[1].DihedralType:continue
            if int_container1 and int_container1[0].int_type == DihedralType:continue
        int_container0 = []
        yield (int_container0, int_container1), done_int_01, flag_dih_params, container_name
# dihedrals extra if dihedral_match_version == 1 (matching 2 middle atoms)
    container_name = 'dihedrals'
    if kwargs.get('dihedral_match_v', 1) == 1:
        flag_dih_params = True, (1, 3), 1
        int_container_01 = ([], [])
        done_int_01 = [], []
        for neigh_at0 in top_01[0].adj[at0]:
            if (top_pair_ind[0], neigh_at0) in sol.sol_atoms_map:
                neigh_at1 = sol.find_common_atom_pair((top_pair_ind[0], neigh_at0), top_pair_ind[1])
                if neigh_at1 not in (None, Dummy):
                    neigh_at_01 = neigh_at0, neigh_at1
                    for i in range(2):
                        #temp_dih_container_01 = neigh_at_01[i].get_container(top_01[i].DihedralType ,flag_class = True, allow_not_found = True)
                        #temp_dih_container_01 = neigh_at_01[i].get_container(container_name, allow_not_found = True)
                        temp_dih_container_01 = top_01[i].atom_interactions[neigh_at_01[i]].get_container\
                            (container_name, allow_not_found = True)
                        if temp_dih_container_01 is not None:
                            for temp_dih in temp_dih_container_01:
                                if neigh_at_01[i] in temp_dih.atoms[1:3]:
                                    int_container_01[i].append(temp_dih)
        yield int_container_01, done_int_01, flag_dih_params, container_name

def gen_int_container_01_list(at0, at1, top_01, top_pair_ind, sol, **kwargs):
    """generator for interaction containers to be checked for matches"""
# all but dihedrals
    done_containers = []
    flag_dih_params = False, (None, None), 2
    for int_container0, container_name in at0.get_interaction_containers(get_container_name=True):
    #for int_container0, container_name in top_01[0].atom_interactions[at0].get_interaction_containers(get_container_name=True):
        done_int_01 = [], []
        if kwargs.get('dihedral_match_v', 1) == 1:
            #if int_container0 and int_container0[0].int_type==top_01[0].DihedralType:continue
            if int_container0 and int_container0[0].int_type == DihedralType:continue
        int_container1 = at1.get_container(container_name, allow_not_found = True)
        if int_container1 is None:
            int_container1 = []
        done_containers.append(container_name)
        yield (int_container0, int_container1), done_int_01, flag_dih_params, container_name
    for int_container1, container_name in at1.get_interaction_containers(get_container_name=True):
    #for int_container1, container_name in top_01[1].atom_interactions[at1].get_interaction_containers(get_container_name=True):
        if container_name in done_containers:
            continue
        done_int_01 = [], []
        if kwargs.get('dihedral_match_v', 1) == 1:
            #if int_container1 and int_container1[0].int_type==top_01[1].DihedralType:continue
            if int_container1 and int_container1[0].int_type == DihedralType:continue
        int_container0 = []
        yield (int_container0, int_container1), done_int_01, flag_dih_params, container_name
# dihedrals extra if dihedral_match_version == 1 (matching 2 middle atoms)
    container_name = 'dihedrals'
    if kwargs.get('dihedral_match_v', 1) == 1:
        flag_dih_params = True, (1, 3), 1
        int_container_01 = ([], [])
        done_int_01 = [], []
        for neigh_at0 in top_01[0].adj[at0]:
            if (top_pair_ind[0], neigh_at0) in sol.sol_atoms_map:
                neigh_at1 = sol.find_common_atom_pair((top_pair_ind[0], neigh_at0), top_pair_ind[1])
                if neigh_at1 not in (None, Dummy):
                    neigh_at_01 = neigh_at0, neigh_at1
                    for i in range(2):
                        #temp_dih_container_01 = neigh_at_01[i].get_container(top_01[i].DihedralType ,flag_class = True, allow_not_found = True)
                        temp_dih_container_01 = neigh_at_01[i].get_container(container_name, allow_not_found = True)
                        if temp_dih_container_01 is not None:
                            for temp_dih in temp_dih_container_01:
                                if neigh_at_01[i] in temp_dih.atoms[1:3]:
                                    int_container_01[i].append(temp_dih)
        yield int_container_01, done_int_01, flag_dih_params, container_name

def find_interactions_ptp(at0, at1, top_01, top_pair_ind, sol, done_interactions, **kwargs):
    temp_gen_out = gen_int_container_01_list_new(at0, at1, top_01, top_pair_ind, sol, **kwargs)
    for int_container_01, done_int_01, flag_dih_params, container_name in temp_gen_out:
        #temp_gen_inner = _find_int_list_01_to_compare(sol, top_01, top_pair_ind, int_container_01, done_int_01,
        temp_gen_inner = _find_int_list_01_to_compare(sol, top_01, top_pair_ind, int_container_01, done_interactions,
                                                      flag_dih_params, **kwargs)
        # generates lists of interactions from top0 and from top1 based on set of common atoms
        # in principle, the groups are of len == 1, but for e.g. dihedrals this could be > 1
        for (int_list0, int_list1), atoms_01 in temp_gen_inner:
            atoms_01_map = [dict(zip(atoms_01[i], atoms_01[(i+1) % 2])) for i in range(2)] # maps atoms0 <-> atoms1
            # this map is needed for impropers, as one can specify the same improper using different order of atoms
            int_matches, int_ptp = _find_interactions_ptp_atom_set(top_01, list(int_list0), list(int_list1),
                                                                   atoms_01_map, **kwargs)
            yield int_matches, int_ptp, container_name

def _get_sol_int_atoms_simple(sol, top_i, int_i, sol_interaction):
    sol_interaction.atoms = tuple(sol.find_sol_atom((top_i, at)) for at in int_i.atoms)

def _get_sol_int_states(sol, sol_int, flag_dih_v_1=False):
    tops2check = []
    sol_int.states = [None] * len(sol.tops)
    for ND_state in sol_int.ND_states:
        for top_i in sol._sol.columns:
            if top_i in sol_int.ND_states[ND_state]:
                sol_int.states[top_i] = ND_state
            else:
                tops2check.append(top_i)
    if flag_dih_v_1:
        flag_atoms_in_sol = _check_sol_dih_v_1_atoms(sol, sol_int, tops2check=tops2check)
    else:
        flag_atoms_in_sol = _check_sol_int_atoms(sol, sol_int, tops2check=tops2check)
    for ND_state in sol_int.ND_states:
        for top_i in sol._sol.columns:
            if top_i not in sol_int.ND_states[ND_state]:
                if flag_atoms_in_sol[top_i]==False:
                    pass
                    #sol_int.states[top_i] = Dummy

def _check_sol_atoms_if_real(sol_row, list_tops2check):
    for top_set in list_tops2check:
        for temp_top_i in top_set:
            if sol_row[temp_top_i] == Dummy:
                return False
    return True

def _check_sol_int_atoms(sol, sol_int, tops2check=None, atom_range=(None, None)):
    """
    checks if atoms of each topology (corresonding to the sol_int.atoms) are physical (non-Dummy) atoms
    :param sol:
    :param sol_int:
    :return: list of bools
    """
    flag_atoms_in_sol = [True] * len(sol.tops)
    if tops2check is None:
        tops2check = sol._sol.columns
    for top_i in tops2check:
        temp_flag = True
        for sol_at in sol_int.atoms[atom_range[0]:atom_range[1]]:
            sol_row = sol._sol.loc[sol_at.id]
            if sol_row[top_i] == Dummy:
                temp_flag = False
                break
        flag_atoms_in_sol[top_i] = temp_flag
    return flag_atoms_in_sol

def _check_sol_dih_v_1_atoms(sol, sol_int, tops2check=None):
    """
    same as _check_sol_int_atoms (but for dihedrals under dihedral_match_v==1 condition)
    :param sol:
    :param sol_int:
    :return:
    """
    # first_check the middle 2 atoms
    if tops2check is None:
        tops2check = sol._sol.columns
    flag_atoms_in_sol = _check_sol_int_atoms(sol, sol_int, tops2check=tops2check, atom_range=(1,3))
    # now check the neighbours of the middle 2, could be any...
    sol_rows = [sol._sol.loc[sol_at.id] for sol_at in sol_int.atoms[1:3]]
    for top_i in tops2check:
        top_atoms = [sol_row[top_i] for sol_row in sol_rows]
        for top_at in top_atoms:
            if flag_atoms_in_sol[top_i]:
                temp_flag = False
                for at_neigh in sol.tops[top_i].adj[top_at]:
                    if at_neigh != Dummy and at_neigh not in top_atoms:
                        temp_flag = True
                        break
                flag_atoms_in_sol[top_i] *= temp_flag
    return flag_atoms_in_sol


def _check_sol_dih_atoms(sol, sol_int):
    for atom_i in (0, 3):
        sol_row = sol._sol.loc[sol_int.atoms[atom_i].id]
        for top_set in sol_int.ND_states.values():
            for temp_top_i in top_set:
                if sol_row[temp_top_i] == Dummy:
                    return False
    return True

def _get_masses_real_atoms(sol_row, list_tops2check):
    masses = []
    for top_set in list_tops2check:
        for temp_top_i in top_set:
            if sol_row[temp_top_i] != Dummy:
                masses.append(sol_row[temp_top_i].m)
    return masses

############################################
# generating making multi state topologies #
############################################

def generate_multi_state_top(sol, top_state = 0, **kwargs):
    """
    generates multi state topology (and ptp topology) based on the solution
    :param sol: solution (from MCS)
    :param top_state: state to use to generate the topology (0 default)
    :kwargs
        flag_copy_undefined_bl: copies undefined_bl from input topologies
        add_DUM_exclusions: adds exclusions betweeen dummies from different states (True by default)
        flag_EP2excl: generates exclusions and pairse from EP
        flag_bond_constraints: generates constraints to all bonds except the perturbed ones

        also to be passed on to the following functions:
            generate_atom_states, generate_int_states, find_atom_order

    """
    sol.toptp.ff = sol.tops[0].ff
    if kwargs.get('flag_copy_undefined_bl', True):
        try:
            if not hasattr(sol.toptp, 'undefined_bl'):
                sol.toptp.undefined_bl = {}
            sol.toptp.undefined_bl.update(sol.tops[0].undefined_bl)
        except:
            pass
    #sol.toptp.top_state = top_state
    generate_atom_states(sol, **kwargs)
    #generate_atom_top_state(sol, **kwargs)
    generate_int_states(sol, **kwargs)
    if kwargs.get('add_DUM_exclusions', True):
        add_DUM_exclusions(sol, **kwargs)
    sol.toptp.set_top_state(top_state = top_state, **kwargs)
    if kwargs.get('flag_EP2excl'):
        generate_excl_from_EP(sol)
    ordered_atoms = find_atom_order(sol, **kwargs)
    c = 0
    for at in ordered_atoms:
        c+=1
        at.sol_id = at.id
        at.id = c
        #at.gr_id = str(at.id)
    if isinstance(sol.toptp.atoms, list):
        sol.toptp.atoms = ordered_atoms
    else:
        sol.toptp.atoms.clear()
        for at in ordered_atoms:
            sol.toptp.add2container(at)
    for at in ordered_atoms:
        sol.toptp.add2container(at.cg, create=True, db_type=list, replace = -1)
    for i, cg in enumerate(sol.toptp.cg):
        cg.update()
        cg.n = i+1
    # add constraints if there is bond ptp
    if sol.toptp.ptp_int.get('bonds') and kwargs.get('flag_bond_constraints', True):
        excl_bonds = []
        for b in sol.toptp.bonds:
            if len(b.ND_states) != 1:
                excl_bonds.append(b)
        sol.toptp.generate_constraints(excl_bonds = excl_bonds)
    get_residues_v1(sol)
    get_res_atom_names_v1(sol)
    get_sys_title(sol)

def generate_atom_states(sol, ff_dumm = None, **kwargs):
    if ff_dumm is None:
        ff_dumm = sol.toptp.get_DUM_type
    for at in sol.toptp.get_atoms():
        # mass
        at.m_states = [None] * len(sol.tops)
        for m_state in at.ND_m_states:
            for top_i in at.ND_m_states[m_state]:
                at.m_states[top_i] = m_state
        # partial charge
        at.p_ch_states = [None] * len(sol.tops)
        for p_ch_state in at.ND_pch_states:
            for top_i in at.ND_pch_states[p_ch_state]:
                at.p_ch_states[top_i] = p_ch_state
        # atom type
        at.a_type_states = [None] * len(sol.tops)
        for a_type_state in at.ND_a_type_states:
            if a_type_state is None:
                for top_i in at.ND_a_type_states[a_type_state]:
                    at.a_type_states[top_i] = ff_dumm
                    assert at.p_ch_states[top_i] == 0.0
                    assert at.m_states[top_i] is None
            else:
                for top_i in at.ND_a_type_states[a_type_state]:
                    at.a_type_states[top_i] = a_type_state

def add_DUM_exclusions(sol, ff_dumm = None, **kwargs):
    if ff_dumm is None:
        ff_dumm = sol.toptp.get_DUM_type
    atoms_non_DUM_states = [list() for i in range(len(sol.tops))]
    for at in sol.toptp.get_atoms(): 
        st_non_dum = [i for i, at_state in enumerate(at.a_type_states) if at_state != ff_dumm]
        if len(st_non_dum)==1:
            st_non_dum = st_non_dum[0]
            atoms_non_DUM_states[st_non_dum].append(at)
    for top_i_pair in combinations(range(len(sol.tops)), 2):
        for at_0 in atoms_non_DUM_states[top_i_pair[0]]:
            for at_1 in atoms_non_DUM_states[top_i_pair[1]]:
                    temp_interaction = Interaction(ExclusionPairType, atoms=(at_0, at_1))
                    for i in range(len(sol.tops)):
                        temp_interaction.add_state(fnc_type = None, params = (True,))
                    sol.toptp.add2container(temp_interaction, create=True, item_id=frozenset((at_0, at_1)), replace = -1)
                    sol.toptp.add_atom_pair2EP_l(at_0, at_1)

def generate_atom_top_state(sol, other_state = None, **kwargs):
    for at in sol.toptp.get_atoms():
        if kwargs.get('flag_EDS_mass'):
            temp_masses = []
            temp_w = []
            for m in at.ND_m_states:
                if m is not None:
                    temp_masses.append(m)
                    temp_w.append(len(at.ND_m_states[m]))
            at.m = np.average(temp_masses, weights=temp_w)
        elif at.m_states[sol.toptp.top_state] is None:
            if other_state:
                if at.m_states[other_state] is None:
                    txt = 'atom mass not defined for neither of the states:\n{:s}  {:d}  {:d}\n'
                    raise Exception(txt.format(str(at), sol.toptp.top_state, other_state))
                at.m = at.m_states[other_state]
            else:
                if len(at.ND_m_states) != 2:
                    raise Exception('atom mass not defined for top_state and has more than 1 ND_states: ' + str(at))
                else:
                    for m in at.ND_m_states:
                        if m is not None:
                            at.m = m
        else:
            at.m = at.m_states[sol.toptp.top_state]
        at.p_ch = at.p_ch_states[sol.toptp.top_state]
        at.a_type = at.a_type_states[sol.toptp.top_state]

def generate_excl_from_EP(sol, **kwargs):
    for atom_pair in sol.toptp.excl_pair:
        temp_EP = sol.toptp.excl_pair[atom_pair]
        if len(temp_EP.ND_states)!=1:
            assert sol.toptp.ptp_excl_pair!=0
            temp_EP.state = ExclusionPairType(params = (1,))

def generate_int_states(sol, **kwargs):
    #generate_atom_states(sol, **kwargs)
    dihedral_match_v = kwargs.get('dihedral_match_v', 1)
    done_interactions = set()
    dihedrals2fix = []
    for ep_int in sol.toptp.excl_pair.values():
        _get_sol_int_states(sol, ep_int)
    for top_i in sol._sol.columns:
        top = sol.tops[top_i]
        for int_cont in top.get_interaction_containers(EP_cont_exclude = True):
            for int_i in int_cont:
                flag_dih_v_1 = False
                if int_i in sol.toptp.interaction_states_map[top_i]:
                    sol_int = sol.toptp.interaction_states_map[top_i][int_i]
                else:
                    sol_int = sol.toptp.Interaction(int_i.int_type)
                    sol_int.ND_states = {int_i.states[0]: {top_i}}
                    sol_int.int_states = [None] * len(sol.tops)
                    sol_int.int_states[top_i] = int_i
                    _get_sol_int_atoms_simple(sol, top_i, int_i, sol_int)
                    _get_sol_int_states(sol, sol_int) # in this case this is also fine for the dihedrals of v==1 (as they have to be DUM in the other states)
                    sol.toptp.add2container(sol_int, create=True, db_type=list, **kwargs)
                    continue
                    #done_interactions.add(sol_int) # not needed as not matched to any other interaction (all other states DUM)
                if sol_int in done_interactions:
                    continue
                if sol_int.int_states[top_i] is None:
                    continue
                flag_OK = True

                if sol_int.int_states[top_i].int_type == ImproperType:
                    for temp_top_i, temp_int_state in enumerate(sol_int.int_states):
                        if temp_int_state and temp_int_state.states[0] in sol_int.ND_states:
                            _get_sol_int_atoms_simple(sol, temp_top_i, temp_int_state, sol_int)
                            break
                #elif sol_int.int_states[top_i].int_type == top.DihedralType:
                elif sol_int.int_states[top_i].int_type == DihedralType:
                    if dihedral_match_v == 1:
                        flag_dih_v_1 = True
                        _get_sol_int_atoms_simple(sol, top_i, int_i, sol_int)
                        if not _check_sol_dih_atoms(sol, sol_int):
                            flag_OK = False
                            dihedrals2fix.append(sol_int)
                            done_interactions.add(sol_int)
                    else:
                        _get_sol_int_atoms_simple(sol, top_i, int_i, sol_int)
                else:
                    _get_sol_int_atoms_simple(sol, top_i, int_i, sol_int)
                if flag_OK:
                    _get_sol_int_states(sol, sol_int, flag_dih_v_1 = flag_dih_v_1)
                    sol.toptp.add2container(sol_int, create=True, db_type=list, **kwargs)
                    done_interactions.add(sol_int)
    sol.toptp.gen_graph()
    for sol_int in dihedrals2fix:
        N_tops2check = 0
        for top_set in sol_int.ND_states.values():
            N_tops2check += len(top_set)
        sol_int.atoms = list(sol_int.atoms)
        flag_OK = [False, False]
        for i, (atom_i, mid_atom_1_i, mid_atom_2_i) in enumerate(((0,1,2), (3,2,1))):
            atom = sol_int.atoms[atom_i]
            sol_row = sol._sol.loc[atom.id]
            if _check_sol_atoms_if_real(sol_row, sol_int.ND_states.values()):
                flag_OK[i] = True
            else:
                best_masses = _get_masses_real_atoms(sol_row, sol_int.ND_states.values())
                best_masses.sort()
                mid_atom_1, mid_atom_2 = sol_int.atoms[mid_atom_1_i], sol_int.atoms[mid_atom_2_i]
                for neigh_atom in sol.toptp.adj[mid_atom_1]:
                    if neigh_atom not in (mid_atom_2, atom):
                        sol_row = sol._sol.loc[neigh_atom.id]
                        temp_masses = _get_masses_real_atoms(sol_row, sol_int.ND_states.values())
                        temp_masses.sort()
                        flag_better = False
                        if len(temp_masses) > len(best_masses):
                            flag_better = True
                            if len(temp_masses) == N_tops2check:
                                flag_OK[i] = True
                                assert _check_sol_atoms_if_real(sol_row, sol_int.ND_states.values())
                        elif len(temp_masses) == len(best_masses) and temp_masses > best_masses:
                            flag_better = True
                        if flag_better:
                            sol_int.atoms[atom_i] = neigh_atom
                            best_masses = temp_masses
        if not flag_OK[0] * flag_OK[1]:
            do_warn('dihedral is not represented in all states by real atoms:\n' + str(sol_int.atoms))
        sol_int.atoms = tuple(sol_int.atoms)
        _get_sol_int_states(sol, sol_int, flag_dih_v_1 = True)
        sol.toptp.add2container(sol_int, **kwargs)

def _add_next_atom_2_stack(sol, tops_order, atoms2order, stack):
    for top_i in tops_order:
        for at in sol.tops[top_i].get_atoms():
            sol_at = sol.find_sol_atom((top_i, at))
            if sol_at in atoms2order:
                stack.append(sol_at)
                return

def find_atom_order(sol, **kwargs):
    tops_order = kwargs.get('tops_order')
    if not tops_order:
        tops_order = [sol.toptp.top_state]
        tops_order.extend([top_i for top_i in sol._sol.columns if top_i != sol.toptp.top_state])
        # starts with sol.toptp.top_state and other ordered
    atoms2order = set(sol.toptp.get_atoms())
    ordered_atoms = []
    stack = deque()
    while atoms2order:
        while stack:
            sol_at = stack.popleft()
            if sol_at in atoms2order:
                sol_atoms = [sol_at]
                sol_atoms_set = {sol_at}
                cg_atoms_set = None
                row = sol._sol.loc[sol_at.id]
                for top_i in tops_order:
                    at_i = row[top_i]
                    if at_i != Dummy:
                        temp_cg_atoms_set = set()
                        for cg_at in at_i.cg.atoms:
                            temp_sol_at = sol.find_sol_atom((top_i, cg_at))
                            if temp_sol_at in atoms2order:
                                temp_cg_atoms_set.add(temp_sol_at)
                                if temp_sol_at not in sol_atoms_set:
                                    sol_atoms.append(temp_sol_at)
                                    sol_atoms_set.add(temp_sol_at)
                        if cg_atoms_set:
                            cg_atoms_set &= temp_cg_atoms_set
                        else:
                            cg_atoms_set = temp_cg_atoms_set
                temp_stack_set = set(stack)
                cg = sol.toptp.ChargeGroup()
                for sol_at in sol_atoms:
                    if sol_at in cg_atoms_set:
                        cg.add_atom(sol_at)
                        ordered_atoms.append(sol_at)
                        atoms2order.remove(sol_at)
                    elif sol_at not in temp_stack_set:
                        stack.append(sol_at)
        _add_next_atom_2_stack(sol, tops_order, atoms2order, stack)
    return ordered_atoms

def get_residues_v1(sol, **kwargs):
    current_residues_atoms = tuple(dict() for top_i in sol._sol.columns)
    current_res = None
    for sol_at in sol.toptp.get_atoms():
        flag_curr_res = True
        for top_i_curr_res in current_residues_atoms:
            for temp_res_atoms in top_i_curr_res.values():
                if temp_res_atoms:
                    flag_curr_res = False
                    break
            if not flag_curr_res:
                break
        if flag_curr_res:
            current_res = None
            current_residues_atoms = tuple(dict() for top_i in sol._sol.columns)
        if current_res is None:
            current_res = sol.toptp.add_residue(flag_type = int)
        current_res.add_atom(sol_at)
        #current_res.gr_id = current_res.id
        for top_i, top_i_curr_res in enumerate(current_residues_atoms):
            top_at = sol._sol.loc[sol_at.sol_id, top_i]
            if top_at not in (None, Dummy):
                if top_at.res not in current_residues_atoms[top_i]:
                    current_residues_atoms[top_i][top_at.res] = set(top_at.res.atoms)
                current_residues_atoms[top_i][top_at.res] -= {top_at}

def get_state_names(sol, **kwargs):
    names = []
    for t in sol.tops:
        if len(t.residues)==1:
            names.append(list(t.get_residues())[0].name)
    if len(set(names))==len(sol.tops):
        return names
    return ['st_{:d}'.format(i) for i in range(1, len(sol.tops)+1)]

def get_res_atom_names_v1(sol, **kwargs):
    if len(sol.tops) == 2:
        ptp_eds_res_name = 'PTPR'
    else:
        ptp_eds_res_name = 'EDSR'
    for res in sol.toptp.get_residues():
        res_atom_names = set()
        temp_residues_all = set()
        for sol_at in res.atoms:
            atom_names = Counter()
            for top_i, top_at in sol._sol.loc[sol_at.sol_id].iteritems():
                if top_at not in (None, Dummy):
                    temp_residues_all.add(top_at.res)
                    atom_names[top_at.name] += 1
            temp_at_name = atom_names.most_common()[0][0]
            if temp_at_name not in res_atom_names:
                res_atom_names.add(temp_at_name)
                sol_at.name = temp_at_name
            else:
                temp_at_name = res._find_next_code(res_atom_names, pref = 'PT')
                res_atom_names.add(temp_at_name)
                sol_at.name = temp_at_name
        res_names = Counter()
        for temp_res in temp_residues_all:
            res_names[temp_res.name] += 1
        temp_res_name_count = res_names.most_common()[0]
        if temp_res_name_count[1] == 1:
            temp_res_name = ptp_eds_res_name
        else:
            temp_res_name = temp_res_name_count[0]
        res.name = temp_res_name

def get_sys_title(sol, **kwargs):
    sys_title_lines = ['PERTURBATION TOPOLOGY BASED ON N = {:} STATES\n\n'.format(len(sol.tops))]
    for top in sol.tops:
        sys_title = top.get_container('sys_title', allow_not_found = True)
        if sys_title:
            sys_title_lines.append(sys_title[0].lines[0])
            sys_title_lines.append('\n')
    dp = DescriptionPart(lines = sys_title_lines)
    sol.toptp.add2container(dp, create=True, db_type=list,   list_index=0)


from SMArt.incl import combinations, np,  OrderedDict, permutations, bisect_left, Counter, DataDumping
from SMArt.alchemy.incl import AlchemicalSolution, Dummy, TopGraphProperties
from SMArt.alchemy.top_matching_fnc import update_ptp, generate_toptp
from SMArt.md.ana.incl import _RMSD, _RMSD_pairwise

class MCS(DataDumping):
    def __init__(self, *tops, flag_partial_ring=True, max_partial_ring_match=2, **kwargs):
        """
        :param tops:
        :param flag_partial_ring - allow partial ring match, such that partial matches form a ring (e.g. benzene and indene)
        :param max_partial_ring_match - max number of allowed atoms to be matched between chain and ring or 2 rings matched partially
        :param kwargs:
            nb - how far (in neighbours) a new atom pair can be (through bonds, see enumerate function) - default 2 (and only option at the moment)
            flag_top_update - update topology at each step
            score_fnc - general score function
            flag_score_fnc - to specify predefined score functions
            add_RMSD - add RMSD to the score ['simple' or 'pairwise'] - see SMArt.md.incl functions
                RMSD_position - at which position in the score to add the RMSD (by default 1)
            kwargs for AlchemicalSolution
                available_atoms
                non_available_atoms
                common_atoms
                tried_pairs
        """
        kwargs = dict(kwargs)
        ######## steps accounting
        self.solutions = [] # solutions from each step
        self.solution_scores = []
        self.c = 0 # counter in which step (call) we are
        self.c_solutions = 0 # debugging
        self._c_solutions = {} # debugging
        self.c_none_solutions = 0 # debugging
        self.none_solutions = [] # debugging
        self.current = [None] * 4 # debugging
        ######## steps accounting
        self.tops = tuple(tops)
        for top in self.tops:
            assert len(list(top.get_Gs())) == 1
        ####### some params
        self.flag_mem = kwargs.get('flag_mem', True)
        self.flag_top_update = kwargs.get('flag_top_update')
        self.flag_top_prune = kwargs.get('flag_top_prune')
        self.flag_score_fnc = kwargs.get('flag_score_fnc', '')
        self.score_fnc = kwargs.get('score_fnc')
        if self.score_fnc is None:
            self.score_fnc = getattr(self, '_calc_sol_score_' + self.flag_score_fnc)
        flag_add_RMSD = kwargs.get('add_RMSD')
        if flag_add_RMSD:
            assert flag_add_RMSD in ('simple', 'pairwise')
            RMSD_fnc_map = dict(simple=_RMSD, pairwise=_RMSD_pairwise)
            self.calc_RMSD = RMSD_fnc_map[flag_add_RMSD]
            self.RMSD_position = kwargs.get('RMSD_position', 1)
            coords = kwargs.get('coords')
            if coords:
                try:
                    for top_i, top in enumerate(self.tops):
                        for at_i, at in enumerate(top.get_atoms()):
                            at.coord = coords[top_i][at_i]
                except:
                    self.coords = {}
                    for top_i, top in enumerate(self.tops):
                        self.coords[top_i] = {}
                        try:
                            for at_i, at in enumerate(top.get_atoms()):
                                self.coords[top_i][at] = coords[top_i][at_i]
                        except:
                            for at_i, at in enumerate(top.adj):
                                self.coords[top_i][at] = coords[top_i][at_i]
        else:
            self.calc_RMSD = False
        ######################################### initial #########################################
        self.flag_partial_ring = flag_partial_ring
        # add chain ring atom pairs in tried_pairs if max_partial_ring_match==0
        self.max_partial_ring_match = max_partial_ring_match
        self.nb = kwargs.get('nb', 2)
        # get list of neighbours for each atom of each top, separated in lists (1st, 2nd, 3rd... neighbours)
        # the list goes to nb - 1
        # this is used for generating new matched pairs
        temp_list = [self.__get_neigh_levels_dependencies(t, **kwargs) for t in self.tops]
        self.neigh_levels, self.dependencies = np.array(list(zip(*temp_list)), dtype=object)
        self.__find_top_dependencies()
        self.__generate_initial_sol_TGPs(**kwargs)
        self._top_pairs = tuple([top_pair_ind for top_pair_ind in combinations(self.initial_sol._sol.columns, 2)])
        self.pyramid_neighbours = list(self.__get_pyramid_neighbours())
        if self.initial_sol._sol.shape[0] and self.flag_top_update:
            initial_sol = self.initial_sol
            self.make_estimates(initial_sol)
            self.calc_score(initial_sol)
            generate_toptp(sol=initial_sol)
        ######################################### initial #########################################
        ######################################### initial topology props #########################################
        ######################################### initial topology props #########################################
        #self.__get_initial_estimate_matrices() # this is for testing!!!!!!!!!!!!!



    ####################################################################################################################
    #                                        initial topology properties
    ####################################################################################################################
    def __get_neigh_levels_dependencies(self, top, **kwargs):
        neigh_level = {}
        dependencies = {}
        for v0 in top.adj:
            temp_neigh_level = {}
            temp_dependencies = []
            for v, l in top.BFS(v0):
                if l == self.nb:
                    break
                if l>0:
                    if l not in temp_neigh_level:
                        temp_neigh_level[l] = []
                    temp_neigh_level[l].append(v)
                    temp_dependencies.append(v)
            neigh_level[v0] = temp_neigh_level
            dependencies[v0] = temp_dependencies
        return neigh_level, dependencies

    def __get_pyramid_neighbours(self, **kwargs):
        # nb is a parameter that defines how far a new atom pair can be - see enumerate function
        for i in range(1, self.nb // 2 + 1):
            for j in range(i, self.nb - i + 1):
                yield (i, j)
                if i!=j:
                    yield (j, i)

    def __find_top_dependencies(self):# this would be needed for separation of MCS search
        """
        for i, top in enumerate(self.tops):
            for int_cont in  top.get_interaction_containers(self):
                pass
        """
        pass

    def __generate_initial_sol_TGPs(self, **kwargs):
        available_atoms = kwargs.get('available_atoms', [None] * len(self.tops))
        non_available_atoms = kwargs.get('non_available_atoms', [None] * len(self.tops))
        self.TGPs = []
        self.tops_neigh_map = []
        sol_available_atoms = []
        for i, t in enumerate(self.tops):
            temp_TPGs = {}
            temp_kwargs = dict(available_Vs=available_atoms[i], non_available_Vs=non_available_atoms[i], rerun_get_rings=True)
            if not self.flag_partial_ring:
                temp_kwargs['flag_get_ind_ring_parts'] = False
            temp_tgp = TopGraphProperties(t, **temp_kwargs)
            temp_TPGs[(temp_tgp.available_Vs, frozenset(), frozenset())] = temp_tgp
            self.TGPs.append(temp_TPGs)
            sol_available_atoms.append(tuple(temp_tgp.available_Vs))
            temp_map = {}
            #for at in temp_tgp.available_Vs:
            for at in t.adj:
                temp_map[at] = dict(t.BFS(at))
                del(temp_map[at][at])
            self.tops_neigh_map.append(temp_map)
        if self.max_partial_ring_match==0:
            tried_pairs = kwargs.get('tried_pairs', [])
            top_non__ring_at = []
            for top_i, top in enumerate(self.tops):
                ring_at, non_ring_at = [], []
                for at in sol_available_atoms[top_i]:
                    if at in top.rings_map:
                        ring_at.append(at)
                    else:
                        non_ring_at.append(at)
                top_non__ring_at.append((ring_at, non_ring_at))
            for top_pair_idx in combinations(range(len(self.tops)), 2):
                top_idx_0 = top_pair_idx[0]
                top_idx_1 = top_pair_idx[1]
                for i in range(2): # loops over ring and non-ring atoms
                    j = (i + 1) % 2 # picks the opposite of i 
                    for at_0 in top_non__ring_at[top_idx_0][i]:
                        for at_1 in top_non__ring_at[top_idx_1][j]:
                            tried_pairs.append(((top_idx_0, at_0), (top_idx_1, at_1)))
            kwargs['tried_pairs'] = tried_pairs

        kwargs['available_atoms'] = sol_available_atoms
        initial_sol = kwargs.get('initial_sol')
        if initial_sol:
            self.initial_sol = initial_sol
        else:
            self.initial_sol = AlchemicalSolution(self.tops, **kwargs)
            if self.flag_top_update:
                self.initial_sol.toptp = generate_toptp(tops = self.tops, **kwargs)

    def get_TGP_top_available_atoms(self, top_i, available_atoms, in_solution_atoms, in_solution_dummy_match_atoms):
        assert isinstance(available_atoms, frozenset) and isinstance(in_solution_atoms, frozenset) and isinstance(in_solution_dummy_match_atoms, frozenset)
        if self.flag_mem:
            k = (available_atoms, in_solution_atoms, in_solution_dummy_match_atoms)
            if k not in self.TGPs[top_i]:
                self.TGPs[top_i][k] = TopGraphProperties(self.tops[top_i], available_Vs=available_atoms,
            in_solution_Vs = in_solution_atoms, in_solution_dummy_match_Vs = in_solution_dummy_match_atoms)
                self.TGPs[top_i][k].get_separate_groups()
            return self.TGPs[top_i][k]
        else:
            tgp = TopGraphProperties(self.tops[top_i], available_Vs=available_atoms, in_solution_Vs = in_solution_atoms,
                                     in_solution_dummy_match_Vs = in_solution_dummy_match_atoms)
            tgp.get_separate_groups()
            return tgp

    ####################################################################################################################
    #                                                   END                                                            #
    #                                        initial topology properties
    ####################################################################################################################

    ####################################################################################################################
    #                                       enumeration part
    #                                       enumeration part
    ####################################################################################################################

    ####################################################################################################################
    #                                    brute force approach
    ####################################################################################################################
    def enumerate_brute_force(self, **kwargs):
        sol = self.initial_sol.copy()
        if sol._sol.shape[0] == 0:
            for temp_pair in self.get_initial_pairs_brute_force(**kwargs):
                new_sol = sol.add_match2sol(temp_pair) ################################ this could already be atom_atom
                self.update_solutions(new_sol)
                ########## debugging
                new_sol.solution_history = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
                ########## debugging
                self.enumerate_brute_force_call(new_sol, **kwargs)
                sol.tried_pairs.append(temp_pair)
        else:
            self.enumerate_brute_force_call(sol, **kwargs)

    def enumerate_brute_force_call(self, sol, **kwargs):
        for temp_pair in self.get_new_pairs_brute_force(sol, **kwargs):
            self.current[0] = sol ################################### debugging
            new_sol = sol.add_match2sol(temp_pair)
            ########## debugging
            new_sol.solution_history = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
            ########## debugging
            self.update_solutions(new_sol)
            self.enumerate_brute_force_call(new_sol, **kwargs)
            sol.tried_pairs.append(temp_pair)

    # actually this might be not good!!! it should go one atom at the time, new line, one atom...
    def get_initial_pairs_brute_force(self, **kwargs):
        sol = self.initial_sol
        for top_pair_ind in combinations(sol._sol.columns, 2):
            for at1 in sol.available_atoms[top_pair_ind[0]]:
                for at2 in sol.available_atoms[top_pair_ind[1]]:
                    temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, at2))
                    if self.__check_atom_atom_pair(temp_pair, sol, **kwargs):
                        yield temp_pair
                    # row atom pair?
                temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, Dummy))
                yield temp_pair
#                break
#            break

    def get_initial_atom_brute_force(self, **kwargs):
        sol = self.initial_sol
        for i, av_atoms in enumerate(sol.available_atoms):
            for at in av_atoms:
                return (i, at), (None, None)

    def get_new_pairs_brute_force(self, sol, **kwargs):
        # rows with None in columns (tops) that are not in sol
        flag_add2row = False
        for row_i, row in sol._sol.iterrows():
            if (row.values == None).any():
                flag_add2row = True
                row_identifier = (None, row_i)
                for top_i in sol._sol.columns:
                    if row[top_i] is None:
                        if sol.flag_top_not_in_sol[top_i]:
                            for temp_atom in sol.available_atoms[top_i]:
                                temp_pair = (row_identifier, (top_i, temp_atom))
                                if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                                    yield temp_pair
                        temp_pair = row_identifier, (top_i, Dummy)
                        if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                            yield temp_pair
#            else:
#                break
#        return

        # new atom-atom pairs from atom-atom pairs in sol -> new rows as well as replacing None in existing rows
        for top_pair_ind in combinations(sol._sol.columns, 2):
            for neigh_pair_level in self.pyramid_neighbours: # loop over the pyramid
            # 1st neighbours, 1st & 2nd, 2nd & 2nd, 2nd & 3rd, 3rd & 3rd, ... till level + level == nb
                top_pair_neigh_levels = self.neigh_levels[list(top_pair_ind)]
                for temp_sol_pair in sol._sol.values[:,top_pair_ind]:
                    if None in temp_sol_pair or Dummy in temp_sol_pair: # skip as None and Dummy don't have neighbours
                        continue
                    if sum([neigh_pair_level[i] in top_pair_neigh_levels[i][temp_sol_pair[i]] for i in range(2)]) != 2:
# skip as in a small molecule (or a small subproblem), not each atom has to have nth level of neighbours
                        continue
                    for temp_at0 in top_pair_neigh_levels[0][temp_sol_pair[0]][neigh_pair_level[0]]:
                        temp_top_i_atom0 = (top_pair_ind[0], temp_at0) # make into top_index_atom
                        if temp_top_i_atom0 in sol.sol_atoms_map:
                            temp_top_i_atom0 = (None, sol.sol_atoms_map[temp_top_i_atom0])
                        elif temp_at0 not in sol.available_atoms[top_pair_ind[0]]: # check if in available
                            continue
                        for temp_at1 in top_pair_neigh_levels[1][temp_sol_pair[1]][neigh_pair_level[1]]:
                            temp_top_i_atom1 = (top_pair_ind[1], temp_at1)
                            if temp_top_i_atom1 in sol.sol_atoms_map:# already a row
                                temp_top_i_atom1 = (None, sol.sol_atoms_map[temp_top_i_atom1])
                                temp_pair = (temp_top_i_atom1, temp_top_i_atom0)
                                if temp_top_i_atom0[0] is None:# both rows!
                                    if not kwargs.get('flag_row_row_match', False):continue
                                    temp_pair = tuple(sorted(temp_pair))
                            elif temp_at1 not in sol.available_atoms[top_pair_ind[1]]: # check if in available
                                continue
                            else:
                                temp_pair = (temp_top_i_atom0, temp_top_i_atom1)
                            if self.check_match_pair(temp_pair, sol, **kwargs):
                                yield temp_pair

        # new row! - first neighbours from atoms in sol + (all atoms from tops not in sol) or (Dummy for each other top)
        for top_i, top_i_sol_atoms in enumerate(sol._sol.values.T):
            top_new_atoms = []
            for top_i_sol_at in top_i_sol_atoms:
                if top_i_sol_at not in (None, Dummy):
                    for new_at in self.neigh_levels[top_i][top_i_sol_at][1]:
                        if new_at in sol.available_atoms[top_i] and new_at not in top_new_atoms:
                            top_new_atoms.append(new_at)
                            for top_j in  sol._sol.columns:
                                if top_i == top_j:
                                    continue
                                if sol.flag_top_not_in_sol[top_j]:
                                    for temp_atom in sol.available_atoms[top_j]:
                                        temp_pair = (top_i, new_at), (top_j, temp_atom)
                                        if self.check_match_pair(temp_pair, sol, **kwargs):
                                            yield temp_pair
                                if 1:
                                    temp_pair = (top_i, new_at), (top_j, Dummy)
                                    if self.check_match_pair(temp_pair, sol, **kwargs):
                                        yield temp_pair
#                            yield (top_i, new_at), (None, None)

    ####################################################################################################################
    #                                          END
    #                                    brute force approach
    ####################################################################################################################

    ####################################################################################################################
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    ####################################################################################################################
    def __check_monotony_full_estiamte(self, sol):################### DEBUGGING
        sol.c_SH_sol = 0
        full_esitmate_full_sol = sol.full_estimate
        sol_SH = sol
        while hasattr(sol_SH, 'SH'):
            sol_SH = sol_SH.SH[0]
            for i, temp_value in enumerate(full_esitmate_full_sol):
                assert temp_value <= sol_SH.full_estimate[i]
            sol.c_SH_sol += 1


    def update_solutions(self, new_solution, n_solutions = -1, **kwargs):
        self.c += 1 ################################## debugging
        if None not in new_solution._sol.values:
            for top_available_atoms in new_solution.available_atoms:
                if top_available_atoms:
                    return
            try:
                flag_find_solution, check_res = self.__find_rings2check_in_sol(new_solution)
                assert check_res[0]
                assert flag_find_solution # actually does not have to be, but for all test cases should be like this
                # so remove once tests done!!!!!!!!
                # here makes more sense - as here is enumeration that does checks on the fly
            except:
                self.make_estimates(new_solution)
                flag_find_solution, check_res = self.__find_rings2check_in_sol(new_solution)
                if not flag_find_solution or not check_res[0]:
                    assert new_solution.full_estimate == 0
                else:
                    assert new_solution.full_estimate != 0
                if new_solution.full_estimate == 0:
                    return
            if not flag_find_solution or not check_res[0]:return ########### not sure if this would be needed
            # but I wrote before that flag_find could actually be 0 even with smart enumeration...
            # possibly if the last match is ring atom to dummy or something like that
            # so maybe keep it...

            ################################################## DEBUGGING
            # check monotony of full_estimates
            self.__check_monotony_full_estiamte(new_solution)
            # check monotony of full_estimates

            self.c_solutions += 1
            if new_solution._sol.shape[0] not in self._c_solutions:
                self._c_solutions[new_solution._sol.shape[0]] = 1
            else:
                self._c_solutions[new_solution._sol.shape[0]] += 1

            if self.solutions:
                ind = bisect_left(self.solution_scores, new_solution.score)
            else:
                ind = 0
            self.solution_scores.insert(ind, new_solution.score)
            self.solutions.insert(ind, new_solution)
            if n_solutions != -1:
                self.solutions = self.solutions[:n_solutions]
            st = kwargs.get('start_time')
            if kwargs.get('verbose') and st and new_solution._sol.shape[0] == len(self.tops[0].adj):
                st[0] += 1
                print(st[0], (st[1].now() - st[2]).seconds)
        else:
#            self.none_solutions.append(new_solution)
            self.c_none_solutions += 1 ################################## debugging

    @staticmethod
    def __merge_top_atom_pair(top_pair_ind, atom_pair):
        return (top_pair_ind[0], atom_pair[0]), (top_pair_ind[1], atom_pair[1])

    ############ checking functions ############
    def __check_bonds2break_make(self, bond2break, bond2make):
        if self.nb==2:
            return False
        ############# missing ################################ fix

    def __find_check_bonds2break_make(self, matched_atom_atom_pair, sol):
        bonds2break, bonds2make = [], []
        top1, at1 = matched_atom_atom_pair[0]
        top2, at2 = matched_atom_atom_pair[1]
        if at1 in (None, Dummy) or at2 in (None, Dummy):
            return True
        for row_i, row in sol._sol.iterrows():
            if row[top1] not in (None, Dummy) and row[top2] not in (None, Dummy):
                if row[top1] in self.tops[top1].adj[at1] and row[top2] not in self.tops[top2].adj[at2]:
                    if not self.__check_bonds2break_make((top1, (at1, row[top1])), (top2, (at2, row[top2]))):
                        return False
                    bonds2make.append((top2, (at2, row[top2])))
                    bonds2break.append((top1, (at1, row[top1])))
                if row[top1] not in self.tops[top1].adj[at1] and row[top2] in self.tops[top2].adj[at2]:
                    if not self.__check_bonds2break_make((top2, (at2, row[top2])), (top1, (at1, row[top1]))):
                        return False
                    bonds2make.append((top1, (at1, row[top1])))
                    bonds2break.append((top2, (at2, row[top2])))
        return True
#        return bonds2break, bonds2make

    def __check_atom_atom_pair(self, matched_pair, sol, **kwargs):
        # check if atoms from both topologies still in available atoms
        # not actually needed, done in the generation of the pair
        for temp_top_ind_atom in matched_pair:
            #assert temp_top_ind_atom[1]!=Dummy ################ maybe this makes sense to ensure with enumeration
            if temp_top_ind_atom[1]!=Dummy:
                assert temp_top_ind_atom[1] in sol.available_atoms[temp_top_ind_atom[0]]
                if temp_top_ind_atom[1] not in sol.available_atoms[temp_top_ind_atom[0]]:
                    return False
        # check if this pair was tried already
        if matched_pair in sol.tried_pairs:
            return False
        # check bonds
        if not self.__find_check_bonds2break_make(matched_pair, sol):
            return False
        # checks passed
        return True

    def __check_row_atom_pair(self, matched_pair, sol, **kwargs):
        # check if atom (second from the pair) still in available atoms
        # not actually needed, done in the generation of the pair
        if matched_pair[1][1]!=Dummy and matched_pair[1][1] not in sol.available_atoms[matched_pair[1][0]]:
            assert 0
            return False
        # check if the ith position (topology) from the proposed atom (top_i_atom) in the proposed row is None
        if sol._sol.loc[matched_pair[0][1], matched_pair[1][0]] != None:
            return False
        # check if this pair was tried already
        if matched_pair in sol.tried_pairs:
            return False
        # check bonds
        temp_row = sol._sol.loc[matched_pair[0][1], :]
        for top_i in sol._sol.columns:
            if temp_row[top_i] not in (None, Dummy):
                matched_atom_atom_pair = (top_i, temp_row[top_i]), matched_pair[1]
                if not self.__find_check_bonds2break_make(matched_atom_atom_pair, sol):
                    return False
        # checks passed
        return True

    @staticmethod
    def __check_if_in_tried(matched_pair, sol, **kwargs):
        return matched_pair in sol.tried_pairs

    def __check_row_row_pair(self, matched_pair, sol, **kwargs):
        # no need to check availability
        # first check that this is not the same row
        if matched_pair[0][1] == matched_pair[1][1]:
            return False
        # check if this pair was tried already
        if matched_pair in sol.tried_pairs:
            return False
        # check that states (topologies) don't have assigned atoms in both rows
        # technically, we do this by checking if these 2 rows have at least 1 None or 2 Dummies in each column
        flag_none = (sol._sol.loc[(matched_pair[0][1], matched_pair[1][1]),:].values == None).any(axis = 0)
        flag_dummy = sol._sol.loc[(matched_pair[0][1], matched_pair[1][1]),:].values == Dummy
        flag_dummy = flag_dummy[0] * flag_dummy[1]
        temp_flags = flag_none + flag_dummy
        if not temp_flags.all():
            return False
        return True

    def check_match_pair(self, matched_pair, sol, **kwargs):
        """
        :param matched_pair:
        :param sol:
        :param kwargs:
        :return:
        """
        # temp_pair is of this format (top_index, atom_key), (top_index, atom_key) -> sorted by top_index
        # or (None, row_index), (top_index, atom_key) -> always None, row_index first
        # or (None, row_index), (None, row_index) -> row index not sorted (not sorted???)
        if matched_pair[0][0] is not None:
            return self.__check_atom_atom_pair(matched_pair, sol, **kwargs)
        elif matched_pair[1][0] is not None:
            return self.__check_row_atom_pair(matched_pair, sol, **kwargs)
        else:
            #return False ############# in principle this is fine - see first comment in the function!
            return self.__check_row_row_pair(matched_pair, sol, **kwargs)
        ############################### missing topology check!!!
        # maybe topology check could be done within the sol class
        #self.top_check_match_pair(self, temp_pair, sol, **kwargs)

    ############ checking functions ############

    ####################################################################################################################
    #                                          END
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    ####################################################################################################################

    ####################################################################################################################
    #                                    stepwise approach - supports sorting
    ####################################################################################################################
    def enumerate_stepwise(self, **kwargs):
        sol = self.initial_sol.copy()
        if sol._sol.shape[0] == 0:
            for temp_pair in self.get_initial_pairs_stepwise(**kwargs):
                new_sol = sol.add_match2sol(temp_pair) ################################ this could already be atom_atom
                self.update_solutions(new_sol)
                ########## debugging
                new_sol.solution_history = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
                ########## debugging
                self.enumerate_stepwise_call(new_sol, **kwargs)
                sol.tried_pairs.append(temp_pair)
        else:
            self.enumerate_stepwise_call(sol, **kwargs)

    def enumerate_stepwise_call(self, sol, **kwargs):
        for temp_pair in self.get_new_pairs_stepwise(sol, **kwargs):
            self.current[0] = sol ################################### debugging
            new_sol = sol.add_match2sol(temp_pair)
            ########## debugging
            new_sol.solution_history = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
            ########## debugging
            self.update_solutions(new_sol)
            self.enumerate_stepwise_call(new_sol, **kwargs)
            sol.tried_pairs.append(temp_pair)

    def enumerate_stepwise_call_partial_sol_invalid(self, sol, top_pair, rows, independent_ring_parts, **kwargs):
        pass

    def get_initial_pairs_stepwise(self, **kwargs):
        sol = self.initial_sol
        for top_pair_ind in combinations(sol._sol.columns, 2):
            for at1 in sol.available_atoms[top_pair_ind[0]]:
                for at2 in sol.available_atoms[top_pair_ind[1]]:
                    temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, at2))
                    if self.__check_atom_atom_pair(temp_pair, sol, **kwargs):
                        yield temp_pair
                temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, Dummy))
                yield temp_pair
                return

    def get_new_pairs_stepwise(self, sol, **kwargs):
        row_atom_pairs = set()
        atom_atom_pairs = set()
        first_row_atom_pair = None
        first_atom_atom_pair = None
        # this is randomized here for test
        """
        BP_pairs = list(self.get_new_pairs_brute_force(sol, **kwargs))
        for temp_pair_i in np.random.choice(range(len(BP_pairs)), len(BP_pairs), False):
            temp_pair = BP_pairs[temp_pair_i]
        """
        for temp_pair in self.get_new_pairs_brute_force(sol, **kwargs):
            if temp_pair[0][0] is None:
                if first_row_atom_pair is None:
                    first_row_atom_pair = temp_pair
                    row_atom_pairs.add(temp_pair)
                elif temp_pair[0][1] == first_row_atom_pair[0][1] and temp_pair[1][0] == first_row_atom_pair[1][0]:
                        row_atom_pairs.add(temp_pair)
            else:
                if first_atom_atom_pair is None:
                    first_atom_atom_pair = temp_pair
                    atom_atom_pairs.add(temp_pair)
                elif first_atom_atom_pair[0] in temp_pair:
                    temp_index = temp_pair.index(first_atom_atom_pair[0])
                    other_atom = temp_pair[(temp_index+1) % 2]
                    if first_atom_atom_pair[1][0] == other_atom[0]:
                        atom_atom_pairs.add(temp_pair)
        if row_atom_pairs:
            for temp_pair in row_atom_pairs:
                yield temp_pair
            """
            row_atom_pairs = list(row_atom_pairs)
            for temp_pair_i in np.random.choice(range(len(row_atom_pairs)), len(row_atom_pairs), replace = False):
                yield row_atom_pairs[temp_pair_i]
            """
            return
        for temp_pair in atom_atom_pairs:
            yield temp_pair
        """
        if atom_atom_pairs:
            atom_atom_pairs = list(atom_atom_pairs)
            for temp_pair_i in np.random.choice(range(len(atom_atom_pairs)), len(atom_atom_pairs), replace = False):
                yield atom_atom_pairs[temp_pair_i]
        """
        return
    # this would have been part of it...
    """
        ## this is a smarter way, but first test the above brute-force-like approach
        flag_yielded = False
        for row_i, row in sol._sol.iterrows():
            if (row.values == None).any():
                row_identifier = (None, row_i)
                # topologies that have a true atom (not None or Dummy) in the row
                tops_flag_atoms = [top_i for top_i in sol._sol.columns if row[top_i] not in (None, Dummy)]
                for top_i in sol._sol.columns:
                    if row[top_i] is None:
                        if sol.flag_top_not_in_sol[top_i]:
                            for temp_atom in sol.available_atoms[top_i]:
                                temp_pair = (row_identifier, (top_i, temp_atom))
                                if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                                    yield temp_pair
                                    flag_yielded = True
                        else:
                            for temp_pair in self.get_new_row_atom_pairs(sol, row, row_identifier, top_i,
                                                                         tops_flag_atoms, **kwargs):
                                yield temp_pair
                                flag_yielded = True
                        temp_pair = row_identifier, (top_i, Dummy)
                        if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                            yield temp_pair
                            flag_yielded = True
                        if flag_yielded:
                            return

    def get_new_row_atom_pairs(self, sol, row, row_identifier, top_i, tops_flag_atoms, **kwargs):
        # top_i is the column where None is
        # looping through columns with true atoms
        for top_j in tops_flag_atoms:
         # looping through levels and lth neighbours of the atom in top_j column of the row
            for temp_l, temp_atom_j in self.neigh_levels[top_j][row[top_j]].items():
                temp_atom_j_pos = np.where(sol._sol[:,top_j] == temp_atom_j)
                # check if temp_atom in sol and if it is matched to a true atom in top_i column
                if temp_atom_j_pos[0] and sol._sol[temp_atom_j_pos[0][0],top_i] not in (None, Dummy):
                    temp_atom_top_i = sol._sol[temp_atom_j_pos[0][0],top_i]
                    for neigh_pair_level in self.pyramid_neighbours:
                        if neigh_pair_level[0]!=temp_l:continue # makes only sense to use the same level as temp_l
                        for lth_neighbour_temp_atom_top_i in self.neigh_levels[top_i][temp_atom_top_i][neigh_pair_level[1]].values():
                            temp_pair = row_identifier, (top_i, lth_neighbour_temp_atom_top_i)
                            if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                                yield temp_pair
    """

    ####################################################################################################################
    #                                          END
    #                                    stepwise approach
    ####################################################################################################################

    ####################################################################################################################
    #                                    stepwise approach with sorting
    ####################################################################################################################

    def __find_ind_ring_part_from_atoms_multy(self, sol, top_i, r, atoms):
        temp_set = set()
        for at in atoms:
            if len(r.independent_parts_map[at]) == 1 and len(r.independent_parts_map[at] - atoms)!=0:
                temp_set.add(r.independent_parts_map[at])
        return temp_set

    def __find_common_ind_ring_part_from_atoms(self, r, atoms):
        temp_iter = iter(atoms)
        first_atom = next(temp_iter)
        temp_set = [ind_ring_part for ind_ring_part in r.independent_parts_map[first_atom]]
        for at in temp_iter:
            temp_set = [ind_ring_part for ind_ring_part in temp_set if at in ind_ring_part]
        return temp_set

    def __find_ind_ring_part_from_atoms(self, sol, top_i, atoms):
        for at in atoms:
            if at not in sol.tops[top_i].rings_map:
                return list()
        temp_set = list()
        for r in sol.tops[top_i].rings_map[atoms[0]]:
            temp_set.extend([ind_ring_part for ind_ring_part in r.independent_parts_map[atoms[0]]])
        for at in atoms[1:]:
            temp_set = [ind_ring_part for ind_ring_part in temp_set if at in ind_ring_part]
        return temp_set

    def __find_rings2check_in_sol_brute_force_2_sides(self, sol, **kwargs): ############# DEBUGGING
        sol_rev = sol.copy()
        sol_rev.tops = sol_rev.tops[::-1]
        self.rev_mcs = getattr(self, 'rev_mcs', MCS(*sol_rev.tops))
        sol_rev.available_atoms = sol_rev.available_atoms[::-1]
        sol_rev.flag_top_not_in_sol = sol_rev.flag_top_not_in_sol[::-1]
        for i, top_i in enumerate(sol._sol.columns):
            top_rev = sol._sol.columns[::-1][i]
            sol_rev._sol.loc[:, top_rev] = sol._sol.loc[:, top_i]
        sol_rev._get_sol_atoms_map()
        self.rev_mcs.make_estimates(sol_rev)
        return self.rev_mcs.__find_rings2check_in_sol(sol_rev, **kwargs)

    def __find_rings2check_in_sol(self, sol, **kwargs): ## assumes connected sol
        min_at_in_sol_allowed = 1
        flag_partial_valid_solution = True
        rings2check = {}
        tops_in_sol = tuple(top_i for top_i in sol._sol.columns if not sol.flag_top_not_in_sol[top_i])
        # find rings with 2 or more atoms in the solution (save the atoms per ring)
        for top_i in tops_in_sol:
            temp_ring_atoms_in_sol = {}
            #for atom_ii, atom_i in enumerate(sol._sol.values[:, top_i]): # loop over atoms of top_i
            for atom_i in sol._sol[top_i]: # loop over atoms of top_i
                if atom_i not in (None, Dummy) and atom_i in sol.tops[top_i].rings_map:
                    for temp_r in sol.tops[top_i].rings_map[atom_i]:
                        if temp_r not in temp_ring_atoms_in_sol:
                            temp_ring_atoms_in_sol[temp_r] = []
                        temp_ring_atoms_in_sol[temp_r].append(atom_i)
            r2del = []
            for temp_r in temp_ring_atoms_in_sol:
                if len(temp_ring_atoms_in_sol[temp_r]) <= self.max_partial_ring_match:
                    r2del.append(temp_r)
                else:
                    #temp_ring_atoms_in_sol[temp_r] = frozenset(temp_ring_atoms_in_sol[temp_r])
                    pass
            for temp_r in r2del:
                del(temp_ring_atoms_in_sol[temp_r])
            rings2check[top_i] = temp_ring_atoms_in_sol
        # rings2check - for each top_i, dict  with keys = rings and values = ring atoms in sol
        tops_ind_ring_parts_pairs = []
        tops_allowed_atoms_pairs = {}
        tops_disallowed_atoms_pairs = {}
        for top_pair_ind in combinations(tops_in_sol, 2):
            flag_partial_valid_solution_top_pair = True
            for r2check_0 in rings2check[top_pair_ind[0]]:
                ind_ring_part_1_skip = []
                atoms_closed_ind_ring_parts_0 = set()
                flag_find_solution = True # if this partial is valid
                atoms_top_0 = {}
                atoms_top_1 = {}
                r2check_1 = None
                for atom_0 in rings2check[top_pair_ind[0]][r2check_0]:
                    atom_1 = sol.find_common_atom_pair((top_pair_ind[0], atom_0), top_pair_ind[1])
                    #temp_row = sol.sol_atoms_map[top_pair_ind[0], atom_0]
                    #atom_1 = sol._sol.values[temp_row, top_pair_ind[1]]
                    if atom_1 not in (None, Dummy):
                        if atom_1 in sol.tops[top_pair_ind[1]].rings_map:
                            #if not atoms_top_0:
                            if r2check_1 is None:
                                r2check_1 = list(sol.tops[top_pair_ind[1]].rings_map[atom_1])
                            else:
                                r2check_1 = [temp_r for temp_r in r2check_1 if atom_1 in temp_r.adj]
                        else:
                            flag_find_solution = False # this means 1 ring matched to non-ring atom (used below if number of matched atoms > max_partial_ring_match)
                        atoms_top_0[atom_0] = atom_1
                        atoms_top_1[atom_1] = atom_0
                if kwargs.get('verbose'):
                    print(atoms_top_0)
                    print(atoms_top_1)
                    print(r2check_0)
                    print(r2check_1)
                if r2check_1 is None:
                    r2check_1 = list()
                if len(r2check_1) != 1: # this means 1 ring matched to 0 or 2 (or more) rings...
                    flag_find_solution = False
                if len(atoms_top_0) > self.max_partial_ring_match: # more matched atoms in a ring than allowed for partial sol
                    temp_flag_partial_valid_solution = False
                    if flag_find_solution == False:
                        return False, None
                    r2check_1 = r2check_1.pop()
                    for atom_1 in rings2check[top_pair_ind[1]][r2check_1]:
                        if atom_1 not in atoms_top_1:
                            atom_0 = sol.find_common_atom_pair((top_pair_ind[1], atom_1), top_pair_ind[0])
                            #temp_row = sol.sol_atoms_map[top_pair_ind[1], atom_1]
                            #atom_0 = sol._sol.values[temp_row, top_pair_ind[0]]
                            if atom_0 not in (None, Dummy): # this means 1 ring (with more matches than allowed for partial sol) matched to non-ring atom
                                return False, None
                    # find all ind ring parts that are compatible with both atoms_top_0 and 1
                    if kwargs.get('verbose'):print('JOS SMO TU')
                    atoms_set_top_0 = frozenset(atoms_top_0)
                    atoms_set_top_1 = frozenset(atoms_top_1)
                    allowed_atoms_0, allowed_atoms_1 = set(), set()
                    for ind_ring_part_0 in r2check_0.independent_parts:
                        at_in_sol_0 = ind_ring_part_0 & atoms_set_top_0#atoms in sol (top_0) and in one of ind_ring_parts
                        if at_in_sol_0:
                            #if len(at_in_sol_0) >= min_at_in_sol_allowed:
                            if 1:
                                allowed_atoms_0 = allowed_atoms_0 | ind_ring_part_0
                            if len(at_in_sol_0) > 2:#this ensures that 1 distinct ind_ring_part is defined (can't be more)
                                at_in_sol_1 = frozenset(atoms_top_0[at_0] for at_0 in at_in_sol_0)
                                ind_ring_part_1 = self.__find_common_ind_ring_part_from_atoms(r2check_1, at_in_sol_1)
                                if len(ind_ring_part_1) != 1:
                                    # in this case, ind_ring_part has to match to exactly 1 ind_ring_part of the other top
                                    return False, None
                                ind_ring_part_1 = ind_ring_part_1[0]
                            else: # this ind_part could be matched to more ind parts
                                ind_ring_part_0 = self.__find_common_ind_ring_part_from_atoms(r2check_0, at_in_sol_0)
                                if len(ind_ring_part_0) != 1:# this subset of atoms defines more ind_ring_parts
                                    continue
                                ind_ring_part_0 = ind_ring_part_0[0]
                                at_in_sol_1 = frozenset(atoms_top_0[at_0] for at_0 in at_in_sol_0)
                                ind_ring_part_1 = self.__find_common_ind_ring_part_from_atoms(r2check_1, at_in_sol_1)
                                if len(ind_ring_part_1) != 1:
                                    continue
                                ind_ring_part_1 = ind_ring_part_1[0]
                            # both cases when len(at_in_sol_0) > 2 or len(at_in_sol_0) == 2
                            if len(ind_ring_part_0) != len(ind_ring_part_1):
                                # ind ring parts of diff len, will never find a full match!
                                return False, None
                            if at_in_sol_1 != ind_ring_part_1 & atoms_set_top_1:
                                # this means that ind_ring_part_1 matches more than ind_ring_part_0
                                return False, None
                            if kwargs.get('verbose'):print('JOS SMO TU 2')
                            ring_ind_part_in_sol_0 = tuple(at_in_sol_0)
                            ring_ind_part_in_sol_1 = tuple(atoms_top_0[at_0] for at_0 in at_in_sol_0)
                            ring_ind_part_in_sol_0_1 = ring_ind_part_in_sol_0, ring_ind_part_in_sol_1
                            atoms_rest_top_0_1 = ind_ring_part_0 - at_in_sol_0, ind_ring_part_1 - at_in_sol_1
                            for i in range(2):
                                tgp = sol.TGPs_tops_in_sol[top_pair_ind][i]
                                if tgp.in_solution_dummy_match_Vs &  atoms_rest_top_0_1[i]:
                                    # this means that some of the rest atoms is matched to Dummy...
                                    ################################## ISSUE #########################
                                    return False, None
                            if kwargs.get('verbose'):print('JOS SMO TU 3')
                            assert atoms_rest_top_0_1 == (ind_ring_part_0 - atoms_set_top_0, ind_ring_part_1 - atoms_set_top_1)
                            if len(atoms_rest_top_0_1[0]) != 0:
                                tops_ind_ring_parts_pairs.append((top_pair_ind, ring_ind_part_in_sol_0_1, atoms_rest_top_0_1))
                            else:
                                atoms_closed_ind_ring_parts_0 = atoms_closed_ind_ring_parts_0 | ind_ring_part_0
                            ind_ring_part_1_skip.append(ind_ring_part_1)
                            #assert len(at_in_sol_1) > min_at_in_sol_allowed
                            allowed_atoms_1 = allowed_atoms_1 | ind_ring_part_1
                            if kwargs.get('verbose'):print('JOS SMO TU 4 END OF THE LOOP')
                    for ind_ring_part_1 in r2check_1.independent_parts:
                        if ind_ring_part_1 in ind_ring_part_1_skip:
                            continue
                        at_in_sol_1 = ind_ring_part_1 & atoms_set_top_1
                        #if len(at_in_sol_1) > min_at_in_sol_allowed:
                        if at_in_sol_1:
                            allowed_atoms_1 = allowed_atoms_1 | ind_ring_part_1
                            if len(at_in_sol_1) > 2:
                                # in this case, ind_ring_part would match more than 1 ind_ring_part of the other top
                                return False, None
                    if len(atoms_closed_ind_ring_parts_0) == len(atoms_top_0):
                        temp_flag_partial_valid_solution = True
                    if top_pair_ind not in tops_allowed_atoms_pairs:
                        tops_allowed_atoms_pairs[top_pair_ind] = [frozenset(), frozenset()]
                    tops_allowed_atoms_pairs[top_pair_ind][0]  = tops_allowed_atoms_pairs[top_pair_ind][0] | allowed_atoms_0
                    tops_allowed_atoms_pairs[top_pair_ind][1]  = tops_allowed_atoms_pairs[top_pair_ind][1] | allowed_atoms_1
                else:
                    # only max_partial_ring_match atoms matched (or less) in a ring!
                    temp_flag_partial_valid_solution = True
                    """
                    if len(atoms_top_0) == 2 and len(r2check_1)==1:
                        r2check_1 = r2check_1.pop()
                        allowed_atoms_0, allowed_atoms_1 = set(), set()
                        ind_ring_parts_0 = self.__find_common_ind_ring_part_from_atoms(r2check_0, atoms_top_0)
                        for ind_ring_part_0 in ind_ring_parts_0:
                            allowed_atoms_0 = allowed_atoms_0 | ind_ring_part_0
                        ind_ring_parts_1 = self.__find_common_ind_ring_part_from_atoms(r2check_1, atoms_top_1)
                        for ind_ring_part_1 in ind_ring_parts_1:
                            allowed_atoms_1 = allowed_atoms_1 | ind_ring_part_1
                        disallowed_atoms_0 = frozenset(r2check_0.adj) - allowed_atoms_0
                        disallowed_atoms_1 = frozenset(r2check_1.adj) - allowed_atoms_1
                        if top_pair_ind not in tops_disallowed_atoms_pairs:
                            tops_disallowed_atoms_pairs[top_pair_ind] = [frozenset(), frozenset()]
                        tops_disallowed_atoms_pairs[top_pair_ind][0] = \
                            tops_disallowed_atoms_pairs[top_pair_ind][0] | disallowed_atoms_0
                        tops_disallowed_atoms_pairs[top_pair_ind][1] = \
                            tops_disallowed_atoms_pairs[top_pair_ind][1] | disallowed_atoms_1
                    """
                flag_partial_valid_solution_top_pair *= temp_flag_partial_valid_solution

            # this is probably good!!!!!!!!
            #"""
            for r2check_1 in rings2check[top_pair_ind[1]]:
                flag_find_solution = True
                atoms_top_0 = {}
                atoms_top_1 = {}
                r2check_0 = None
                for atom_1 in rings2check[top_pair_ind[1]][r2check_1]:
                    atom_0 = sol.find_common_atom_pair((top_pair_ind[1], atom_1), top_pair_ind[0])
                    #temp_row = sol.sol_atoms_map[top_pair_ind[1], atom_1]
                    #atom_0 = sol._sol.values[temp_row, top_pair_ind[0]]
                    if atom_0 not in (None, Dummy):
                        if atom_0 in sol.tops[top_pair_ind[0]].rings_map:
                            #if not atoms_top_1:
                            if r2check_0 is None:
                                r2check_0 = list(sol.tops[top_pair_ind[0]].rings_map[atom_0])
                            else:
                                r2check_0 = [temp_r for temp_r in r2check_0 if atom_0 in temp_r.adj]
                        else:
                            flag_find_solution = False  # this means 1 ring matched to non-ring atom (used below if number of matched atoms > max_partial_ring_match)
                        atoms_top_0[atom_0] = atom_1
                        atoms_top_1[atom_1] = atom_0
                if r2check_0 is None:
                    r2check_0 = list()
                if len(r2check_0) != 1: # this means 1 ring matched to 0 or 2 (or more) rings...
                    flag_find_solution = False
                if len(atoms_top_0) > self.max_partial_ring_match: # more matched atoms in a ring than allowed for partial sol
                    #temp_flag_partial_valid_solution = False
                    if flag_find_solution == False:
                        return False, None
                    r2check_0 = r2check_0.pop()
                    for atom_0 in rings2check[top_pair_ind[0]][r2check_0]:
                        if atom_0 not in atoms_top_0:
                            atom_1 = sol.find_common_atom_pair((top_pair_ind[0], atom_0), top_pair_ind[1])
                            #temp_row = sol.sol_atoms_map[top_pair_ind[0], atom_0]
                            #atom_1 = sol._sol.values[temp_row, top_pair_ind[1]]
                            if atom_1 not in (None, Dummy): # this means 1 ring (with more matches than allowed for partial sol) matched to non-ring atom
                                return False, None
                    #flag_partial_valid_solution_top_pair *= temp_flag_partial_valid_solution
                    # in this case this is not needed, as this is already checked above
                    # this needs to be check only when all against all atoms are from ring1 and ring2 (no ring1 atoms matched to non-ring2 atoms)
            #"""

            if top_pair_ind in tops_allowed_atoms_pairs:
                if not flag_partial_valid_solution_top_pair:
                    assert tops_allowed_atoms_pairs[top_pair_ind][0] and tops_allowed_atoms_pairs[top_pair_ind][1]
                else:
                    del(tops_allowed_atoms_pairs[top_pair_ind])
            if top_pair_ind in tops_disallowed_atoms_pairs:
                if not tops_disallowed_atoms_pairs[top_pair_ind][0] and not tops_disallowed_atoms_pairs[top_pair_ind][1]:
                    del(tops_disallowed_atoms_pairs[top_pair_ind])
            flag_partial_valid_solution *= flag_partial_valid_solution_top_pair
        res = (flag_partial_valid_solution,
               tops_ind_ring_parts_pairs, tops_allowed_atoms_pairs, tops_disallowed_atoms_pairs)
        # tops_ind_ring_parts_pairs pairs of partial ind ring parts matches - if exists, only this
        # tops_allowed_atoms_pairs - if flag_partial_valid_solution is False, these are allowed atoms that can be added
        return True, res

    def __find_partial_ind_ring_parts_in_sol(self, sol, **kwargs): ## assumes connected sol
        flag_find_solution = True
        ind_ring_parts2check = {}
        tops_in_sol = tuple(top_i for top_i in sol._sol.columns if not sol.flag_top_not_in_sol[top_i])
        for top_i in tops_in_sol:
            temp_ring_atoms_in_sol = []
            #for atom_ii, atom_i in enumerate(sol._sol.values[:, top_i]): # loop over atoms of top_i
            for atom_i in sol._sol[top_i]: # loop over atoms of top_i
                if atom_i not in (None, Dummy) and atom_i in sol.tops[top_i].rings_map:
                    temp_ring_atoms_in_sol.append(atom_i)
            temp_ring_atoms_in_sol = frozenset(temp_ring_atoms_in_sol)
            temp_ind_ring_parts2check = {}
            for atom_i in temp_ring_atoms_in_sol:
                for r in sol.tops[top_i].rings_map[atom_i]:
                    for ind_ring_part in r.independent_parts_map[atom_i]:
                        temp_ind_ring_part_in_sol = ind_ring_part & temp_ring_atoms_in_sol
                        if len(temp_ind_ring_part_in_sol) > 2:
                            temp_ind_ring_parts2check[ind_ring_part] = temp_ind_ring_part_in_sol
            ind_ring_parts2check[top_i] = temp_ind_ring_parts2check
        tops_ind_ring_parts_pairs = []
        for top_pair_ind in combinations(tops_in_sol, 2):
            ind_ring_part_1_skip = []
            for ind_ring_part_0 in ind_ring_parts2check[top_pair_ind[0]]:
                atoms_top_0 = []
                atoms_top_1 = []
                for atom_0 in ind_ring_parts2check[top_pair_ind[0]][ind_ring_part_0]:
                    atom_1 = sol.find_common_atom_pair((top_pair_ind[0], atom_0), top_pair_ind[1])
                    #temp_row = sol.sol_atoms_map[top_pair_ind[0], atom_0]
                    #atom_1 = sol._sol.values[temp_row, top_pair_ind[1]]
                    if atom_1 not in (None, Dummy):
                        atoms_top_0.append(atom_0)
                        atoms_top_1.append(atom_1)
                if len(atoms_top_1) > 2:
                    ind_ring_part_1 = self.__find_ind_ring_part_from_atoms(sol, top_pair_ind[1], atoms_top_1)
                    if len(ind_ring_part_1) != 1: # 1 independent ring part does not map to 1 ind ring parts in other top...
                        flag_find_solution = False
                        return flag_find_solution, None
                    ind_ring_part_1 = ind_ring_part_1[0]
                    if len(ind_ring_part_0) != len(ind_ring_part_1):
                        # ind ring parts of diff len, will never find a full match!
                        flag_find_solution = False
                        return flag_find_solution, None
                    for atom_1 in ind_ring_parts2check[top_pair_ind[1]][ind_ring_part_1]:
                        if atom_1 not in atoms_top_1:
                            temp_atom_0 = sol.find_common_atom_pair((top_pair_ind[1], atom_1), top_pair_ind[0])
                            #temp_row = sol.sol_atoms_map[top_pair_ind[1], atom_1]
                            #temp_atom_0 = sol._sol.values[temp_row, top_pair_ind[0]]
                            if temp_atom_0 not in (None, Dummy):
                                # 1 independent ring part does not map to 1 ind ring parts in other top...
                                flag_find_solution = False
                                return flag_find_solution, None
                    atoms_rest_top_0_1 = ind_ring_part_0 - frozenset(atoms_top_0), ind_ring_part_1 - frozenset(atoms_top_1)
                    if len(atoms_rest_top_0_1[0]) != 0:
                        tops_ind_ring_parts_pairs.append((top_pair_ind, (atoms_top_0, atoms_top_1), atoms_rest_top_0_1))
                    ind_ring_part_1_skip.append(ind_ring_part_1)
            for ind_ring_part_1 in ind_ring_parts2check[top_pair_ind[1]]:
                if ind_ring_part_1 in ind_ring_part_1_skip:
                    continue
                atoms_top_0 = []
                for atom_1 in ind_ring_parts2check[top_pair_ind[1]][ind_ring_part_1]:
                    atom_0 = sol.find_common_atom_pair((top_pair_ind[1], atom_1), top_pair_ind[0])
                    #temp_row = sol.sol_atoms_map[top_pair_ind[1], atom_1]
                    #atom_0 = sol._sol.values[temp_row, top_pair_ind[0]]
                    if atom_0 not in (None, Dummy):
                        atoms_top_0.append(atom_0)
                if len(atoms_top_0) > 2:
                    ind_ring_part_0 = self.__find_ind_ring_part_from_atoms(sol, top_pair_ind[0], atoms_top_0)
                    assert len(ind_ring_part_1) != 1
                    # 1 independent ring part does not map to 1 ind ring parts in other top...
                    flag_find_solution = False
                    return flag_find_solution, None
        return flag_find_solution, tops_ind_ring_parts_pairs

    def enumerate_stepwise_sorted(self, **kwargs):
        sol = self.initial_sol.copy()
        if sol._sol.shape[0] == 0:
            for temp_pair in self.get_sorted_initial_pairs_stepwise_2(sol, **kwargs):
                if kwargs.get('verbose'):
                    print('1', self.c, temp_pair)
                    if self.c>1:exit()
                new_sol = sol.add_match2sol(temp_pair)
                new_sol = self.update_solution(new_sol, temp_pair, sol, **kwargs)
                if new_sol:
                    ########## debugging
                    new_sol.SH = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
                    ########## debugging
                    self.update_solutions(new_sol, **kwargs)
                    assert new_sol.full_estimate <= sol.full_estimate
                    self.enumerate_stepwise_sorted_call(new_sol, **kwargs)
                sol.tried_pairs.append(temp_pair)
        else:
            self.enumerate_stepwise_sorted_call(sol, **kwargs)

    def enumerate_stepwise_sorted_call(self, sol, **kwargs):
        for temp_pair in self.get_sorted_new_pairs_stepwise(sol, **kwargs):
            self.current[0] = sol ################################### debugging
            if kwargs.get('verbose'):
                print('2', self.c, temp_pair)
                if self.c>10:exit()
            new_sol = sol.add_match2sol(temp_pair)
            self.current[1:3] = new_sol, temp_pair
            new_sol = self.update_solution(new_sol, temp_pair, sol, **kwargs)
            self.current[3] = new_sol
            if new_sol:
                ########## debugging
                new_sol.SH = sol, sol.copy(), new_sol.copy(), temp_pair, self.c
                ########## debugging
                self.update_solutions(new_sol, **kwargs)
                #assert new_sol.full_estimate <= sol.full_estimate
                self.enumerate_stepwise_sorted_call(new_sol, **kwargs)
            sol.tried_pairs.append(temp_pair)

    def update_solution(self, sol, temp_pair, old_sol, **kwargs): # this can go in update_solutions
        self.make_estimates(sol)
        if not sol.flag_find_solution:
            # if some of the ring pairs does not have a valid solution for closure...
            return False

        assert temp_pair[0][1] != Dummy

        # topology update (interactions)
        if self.flag_top_update:
            old_sol.update_copy_sol_toptp(sol)
            update_ptp(sol, temp_pair, old_sol, **kwargs)

        # propagate dummies
        if not kwargs.get('flag_skip_update_solution'):
            ############################# IMPROVEMENT POSSIBLE ##########################
            # if a temp pair does not have Dummy, it still could be subject to dummies propagation!!!!!!!!!!
            # this should be when an atom from the pair has a first or second neighbour matched to dummy
            dummies_matches = []
            if temp_pair[1][1] == Dummy:
                top_j = temp_pair[1][0]
                if temp_pair[0][0] == None:
                    row = temp_pair[0][1]
                    top_atom_list_2_propagate = [(top_i, sol._sol.loc[row, top_i])
                                        for top_i in sol._sol.columns if sol._sol.loc[row, top_i] not in (None, Dummy)]
                else:
                    top_atom_list_2_propagate = (temp_pair[0],)
                for top_i, atom_i in top_atom_list_2_propagate:
                    atoms_i = []
                    if atom_i in sol.tops[top_i].rings_map:
                        for row in sol._sol.values:
                            if row[top_i] not in (None, Dummy) and row[top_j] not in (None, Dummy):
                                atoms_i.append(row[top_i])
                        atoms_i = frozenset(atoms_i)
                    else:# same looping, but with break as only 1 pair is enough to do the dummy propagation
                        for row in sol._sol.values:
                            if row[top_i] not in (None, Dummy) and row[top_j] not in (None, Dummy):
                                atoms_i.append(row[top_i])
                                break
                    if atoms_i:
                        if top_i < top_j:
                            top_pair_ind = (top_i, top_j)
                            tgp = sol.TGPs_tops_in_sol[top_pair_ind][0]
                        else:
                            top_pair_ind = (top_j, top_i)
                            tgp = sol.TGPs_tops_in_sol[top_pair_ind][1]
                        for at_i in tgp.dummy_branches[atom_i]:
                            temp_pair = (top_i, at_i), (top_j, Dummy)
                            temp_pair = sol.transform_top_index_atom_pair(temp_pair)
                            if temp_pair not in dummies_matches:
                                if self.check_match_pair(temp_pair, sol, **kwargs):
                                    dummies_matches.append(temp_pair)
                                else:
                                    return False
            for temp_pair in dummies_matches:
                new_sol = sol.add_match2sol(temp_pair)
                if self.flag_top_update:
                    new_sol.toptp  = sol.toptp
                    update_ptp(new_sol, temp_pair, None, **kwargs)
                sol = new_sol
            if dummies_matches:
                self.make_estimates(sol)
                if not sol.flag_find_solution:
                    # if some of the ring pairs does not have a valid solution for closure...
                    return False
        if self.prune(sol, **kwargs):
            return False
        return sol

    def _prune_top_EDS(self, sol, **kwargs):
        if sol.toptp.ptp_int or sol.toptp.ptp_make_break_int:
            return True
        elif  kwargs.get('flag_prune_EDS_match_mass', True) and sol.toptp.ptp_atom[0] != 0.:
            return True
        else:
            return False

    def _prune_top_bond(self, sol, **kwargs):
        if sol.toptp.ptp_int.get('bonds') or sol.toptp.ptp_make_break_int:
            return True

    @staticmethod
    def _calc_sol_score_(sol, **kwargs): # the lower the better
        num_matches = sum(sol.full_estimate)
        if hasattr(sol, 'toptp'):
            top_score = sum(sol.toptp.ptp_atom) + sum(sol.toptp.ptp_int.values()) + sol.toptp.ptp_excl_pair
            sol.score = (-num_matches, top_score)
        else:
            sol.score = -num_matches

    @staticmethod
    def _calc_sol_score_EDS(sol, **kwargs): # the lower the better
        sol.score = (-sum(sol.full_estimate), sol.toptp.ptp_atom[1], sol.toptp.ptp_atom[2])

    def __get_c_space(self, sol):
        shape = list(sol._sol.shape[::-1])
        shape.append(3)
        c_space = np.zeros(shape)
        for row_i, row in sol._sol.iterrows():
            for top_i, at in enumerate(row):
                if at is None:
                    at_coord = np.nan
                else:
                    try:
                        at_coord = at.coord
                    except:
                        at_coord = self.coords[top_i][at]
                c_space[top_i, row_i,:] = at_coord
        return c_space
        
    def __add_RMSD_score(self, sol, **kwargs):
        c_space = self.__get_c_space(sol)
        RMSD = self.calc_RMSD(c_space, flag_make_c_space=False)
        sol._score_without_RMSD = sol.score
        try:
            sol.score = list(sol.score)
        except:
            sol.score = [sol.score]
        if  0 <= self.RMSD_position < len(sol.score):
            sol.score.insert(self.RMSD_position, RMSD)
        else:
            sol.score.append(RMSD)
        #sol.score = tuple(sol.score)

    def calc_score(self, sol, **kwargs):
        self.score_fnc(sol, **kwargs)
        if self.calc_RMSD:
            self.__add_RMSD_score(sol, **kwargs)
        #score_fnc = getattr(self, '_calc_sol_score_' + self.flag_score_fnc)
        #score_fnc(sol)

    def prune(self, sol, **kwargs):
        self.calc_score(sol, **kwargs) # actually the same as the line below...
        if not kwargs.get('flag_skip_prune'):
            if self.solutions:
                if hasattr(sol, '_score_without_RMSD'):
                    score_attr = '_score_without_RMSD'
                else:
                    score_attr = 'score'
                if getattr(self.solutions[0], score_attr) < getattr(sol, score_attr):
                    return True
            if self.flag_top_prune:
                top_prune_fnc = getattr(self, '_prune_top_' + self.flag_top_prune)
                if top_prune_fnc(sol, **kwargs):
                    return True

    @staticmethod
    def __get_best_estimate_pos(pair_estimate_matrix):
        best_estimate_pos = np.unravel_index(pair_estimate_matrix.argmax(), pair_estimate_matrix.shape)
        return pair_estimate_matrix[best_estimate_pos], best_estimate_pos

    @staticmethod
    def __get_best_estimate_row_pos(pair_estimate_row_col):
        best_estimate_pos = pair_estimate_row_col.argmax()
        return pair_estimate_row_col[best_estimate_pos], best_estimate_pos

    @staticmethod
    def __get_overall_best_estimate_pos(top_pair_best_estimates):
        best_estimate_pos = np.argmax(top_pair_best_estimates)
        return top_pair_best_estimates[best_estimate_pos], best_estimate_pos

    @staticmethod
    def __get_best_estimate_row_sorted(best_estimate_row):
        #best_estimate_row_sorted_index = (np.array(best_estimate_row) * -1).argsort()
        best_estimate_row_sorted_index = (best_estimate_row * -1).argsort() # assumes that best_est_row is np.array
        return best_estimate_row_sorted_index

    def get_sorted_initial_pairs_stepwise(self, sol, **kwargs):
        ######## things to improve:
        # for each atom, make a list (N-1) tops of best estimate of this atom with that top
        # sort this list of atoms
        # i.e. find atom that has best (max) list of best estimates
        self.make_estimates(sol)
        top_pair_best_estimates = []
        top_pair_best_estimates_pos = []
        for top_pair_ind in self._top_pairs:
            noDummy_estimate_matrix = sol.pair_estimates[top_pair_ind][:-1, :-1]
            temp_best_estimate, temp_best_estimate_pos = self.__get_best_estimate_pos(noDummy_estimate_matrix)
            top_pair_best_estimates.append(temp_best_estimate)
            top_pair_best_estimates_pos.append(temp_best_estimate_pos)
        temp_best_estimate, temp_best_estimate_pos = self.__get_overall_best_estimate_pos(top_pair_best_estimates)
        top_pair_ind = self._top_pairs[temp_best_estimate_pos]
        top_pair_best_estimate_pos = top_pair_best_estimates_pos[temp_best_estimate_pos]
        best_estimate_row = sol.pair_estimates[top_pair_ind][top_pair_best_estimate_pos[0], :]
        best_estimate_row_sorted_index = self.__get_best_estimate_row_sorted(best_estimate_row)
        assert best_estimate_row_sorted_index[0] == top_pair_best_estimate_pos[1]
        at1 = sol.available_atoms[top_pair_ind[0]][top_pair_best_estimate_pos[0]]
        c_dummy = 0
        for top2_available_atom_index in best_estimate_row_sorted_index:
            try:
                at2 = sol.available_atoms[top_pair_ind[1]][top2_available_atom_index]
            except IndexError:
                at2 = Dummy
                c_dummy += 1
            temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, at2))
            yield temp_pair
        assert c_dummy==1

    @staticmethod
    def __get_initial_best_atom_estimate(sol):
        return np.zeros(len(sol._sol.columns) - 1)

    @staticmethod
    def __check_best_atom_estimate_dummy_flag(temp_best_atom_estimate, temp_best_atom_estimate_dummy_flag,
                                              temp_best_atom_estimate_sorted_index):
        for i in temp_best_atom_estimate_sorted_index:
            if temp_best_atom_estimate[i] == temp_best_atom_estimate[temp_best_atom_estimate_sorted_index[0]]:
                if temp_best_atom_estimate_dummy_flag[i] == False:
                    return True, i
            else:
                return False, None
        return False, None

    def __compare_best_atom_estimate(self, curr_max, curr_max_dummy_count,
                                     temp_best_atom_estimate, temp_best_atom_estimate_dummy_flag):
        temp_best_atom_estimate_sorted_index = self.__get_best_estimate_row_sorted(temp_best_atom_estimate)
        # first we check that at least one of the best matches is with a non-Dummy atom
        flag_check_dummy_flag, best_top_j = self.__check_best_atom_estimate_dummy_flag(temp_best_atom_estimate,
                                            temp_best_atom_estimate_dummy_flag, temp_best_atom_estimate_sorted_index)
        # best_top_j is the top index of the best estimate with a non-Dummy atom
        # note that this is true index if < top_i, but has needs top_j + 1 to become the correct index if top_j >= top_i
        if flag_check_dummy_flag:
            # comparing each element at the time (the estimates have to be ordered)
            curr_temp_zip = zip(curr_max, temp_best_atom_estimate[temp_best_atom_estimate_sorted_index])
            for curr_max_estimate, temp_estimate in curr_temp_zip:
                if temp_estimate > curr_max_estimate:
                    return True, (temp_best_atom_estimate_sorted_index, sum(temp_best_atom_estimate_dummy_flag), best_top_j)
                elif temp_estimate < curr_max_estimate:
                    return False, None
            if curr_max_dummy_count > sum(temp_best_atom_estimate_dummy_flag):
                return True, (temp_best_atom_estimate_sorted_index, sum(temp_best_atom_estimate_dummy_flag), best_top_j)
        return False, None

    def get_sorted_initial_pairs_stepwise_2(self, sol, **kwargs):
        # uses the same estimates as above, but takes into account not only pair with the best estimate, but
        # looks for the atom with best list of estimates with all other tops
        #   for each atom, find best estimates with other tops
        #   sort these list of best estimates (default as tuples of nums would be sorted; or by sum)
        self.make_estimates(sol)
        best_atom_estimate = self.__get_initial_best_atom_estimate(sol) # list of best estimates with different tops
        curr_max_dummy_count = len(sol._sol.columns)
        curr_max_top_pair_atom_index = None
        for top_i in sol._sol.columns:
            for top_i_atom_i  in range(len(sol.available_atoms[top_i])):
                temp_best_atom_estimate = []
                temp_best_atom_estimate_dummy_flag = [] # stores True/False if the best estimate is with a Dummy
                for top_j in sol._sol.columns:
                    if top_i == top_j:
                        continue
                    if top_i < top_j:
                        top_pair_ind = (top_i, top_j)
                        temp_row_col = sol.pair_estimates[top_pair_ind][top_i_atom_i, :]
                    else:
                        top_pair_ind = (top_j, top_i)
                        temp_row_col = sol.pair_estimates[top_pair_ind][:, top_i_atom_i]
                    temp_best_estimate, temp_best_estimate_pos = self.__get_best_estimate_row_pos(temp_row_col)
                    temp_best_atom_estimate.append(temp_best_estimate)
                    temp_best_atom_estimate_dummy_flag.append(temp_best_estimate_pos == len(temp_row_col) - 1)
                temp_best_atom_estimate = np.array(temp_best_atom_estimate)
                flag_compare2curr_max, additional_data = self.__compare_best_atom_estimate(best_atom_estimate,
                                    curr_max_dummy_count, temp_best_atom_estimate, temp_best_atom_estimate_dummy_flag)
                if flag_compare2curr_max:
                    temp_best_atom_estimate_sorted_index, temp_best_dummy_count, best_top_j = additional_data
                    best_atom_estimate = temp_best_atom_estimate[temp_best_atom_estimate_sorted_index]
                    curr_max_dummy_count = temp_best_dummy_count
                    curr_max_top_pair_atom_index = top_i, top_i_atom_i, best_top_j
        if curr_max_top_pair_atom_index:
            top_i, top_i_atom_i, best_top_j = curr_max_top_pair_atom_index
            if best_top_j >= top_i:
                best_top_j += 1 # as len(best_atom_estimate) is len(tops) - 1, and best_top_j == top_i makes no sense
                top_pair_ind = (top_i, best_top_j)
                best_estimate_row_col = sol.pair_estimates[top_pair_ind][top_i_atom_i, :]
                at1 = sol.available_atoms[top_i][top_i_atom_i]
                flag_at_flip = 1
            else:
                top_pair_ind = (best_top_j, top_i)
                best_estimate_row_col = sol.pair_estimates[top_pair_ind][:, top_i_atom_i]
                at1 = sol.available_atoms[top_i][top_i_atom_i]
                flag_at_flip = -1
            best_estimate_row_col_sorted_index = self.__get_best_estimate_row_sorted(best_estimate_row_col)
            c_dummy = 0
            for top_other_available_atom_index in best_estimate_row_col_sorted_index:
                try:
                    at2 = sol.available_atoms[best_top_j][top_other_available_atom_index]
                except IndexError:
                    at2 = Dummy
                    c_dummy += 1
                temp_pair = self.__merge_top_atom_pair(top_pair_ind, (at1, at2)[::flag_at_flip])
                if self.check_match_pair(temp_pair, sol, **kwargs):
                    yield temp_pair
            assert c_dummy==1

    def __get_unsorted_new_pairs_stepwise_1_part(self, sol, **kwargs):
# filling the None places in the current sol
# generate pairs of atoms in sol and all other atoms from tops not in sol,
# also add Dummy pairs where None is in a pair with one of the Vs_of_interest
        for row_i, row in sol._sol.iterrows():
            if (row.values == None).any():
                row_identifier = (None, row_i)
                row_tops_in_sol, row_tops_none = [], []
                for top_i in sol._sol.columns:
                    if row[top_i] is None:
                        row_tops_none.append(top_i)
                    elif row[top_i] != Dummy:
                        row_tops_in_sol.append(top_i)
                for top_i in row_tops_none:
                    if sol.flag_top_not_in_sol[top_i]:
                        for temp_atom in sol.available_atoms[top_i]:
                            temp_pair = (row_identifier, (top_i, temp_atom))
                            if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                                yield temp_pair
                        temp_pair = row_identifier, (top_i, Dummy)
                        if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                            yield temp_pair
                    else:
                        temp_pair = None
                        for top_i_in_sol in row_tops_in_sol:
                            if top_i < top_i_in_sol:
                                #atoms_top_in_sol = top_pair_atoms[(top_i, top_i_in_sol)][1]
                                tgp = sol.TGPs_tops_in_sol[(top_i, top_i_in_sol)][1]
                            else:
                                tgp = sol.TGPs_tops_in_sol[(top_i_in_sol, top_i)][0]
                                #atoms_top_in_sol = top_pair_atoms[(top_i_in_sol, top_i)][0]
                            #if row[top_i_in_sol] in atoms_top_in_sol:
                            # this is somehow wrong with atoms_top_in_sol...
                            if row[top_i_in_sol] in tgp.Vs_of_interest:
                                temp_pair = row_identifier, (top_i, Dummy)
                                break
                        if temp_pair and self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                            yield temp_pair

    def __get_unsorted_new_pairs_stepwise_1_part_old(self, sol, **kwargs):
# filling the None places in the current sol, generate pairs of atoms in sol and all other atoms from tops not in sol
        for row_i, row in sol._sol.iterrows():
            if (row.values == None).any():
                row_identifier = (None, row_i)
                for top_i in sol._sol.columns:
                    if row[top_i] is None:
                        if sol.flag_top_not_in_sol[top_i]:
                            for temp_atom in sol.available_atoms[top_i]:
                                temp_pair = (row_identifier, (top_i, temp_atom))
                                if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                                    yield temp_pair
                        temp_pair = row_identifier, (top_i, Dummy)
                        #######################################
                        ################################################################### FIX!!!!!!!!!!!!!!!
                        # this is wrong if no atoms in the row that are first neighbours to atoms that are matched to top_i!!!!!!!
                        if self.__check_row_atom_pair(temp_pair, sol, **kwargs):
                            yield temp_pair

    def __get_unsorted_new_pairs_stepwise_2_part(self, sol, **kwargs):
        # new atom-atom pairs from atom-atom pairs in sol -> new rows as well as replacing None in existing rows
        top_pairs_ring_allowed_atoms = kwargs.get('top_pairs_ring_allowed_atoms')
        top_pairs_ring_disallowed_atoms = kwargs.get('top_pairs_ring_disallowed_atoms')
        tops_in_sol = tuple(top_i for top_i in sol._sol.columns if not sol.flag_top_not_in_sol[top_i])
        for top_pair_ind in combinations(tops_in_sol, 2):
            #top_pair_atoms[top_pair_ind] = [set(), set()]
            if top_pairs_ring_allowed_atoms and top_pair_ind not in top_pairs_ring_allowed_atoms:
                assert top_pair_ind[::-1] not in top_pairs_ring_allowed_atoms
                continue # only allowed top pairs to be considered
            if top_pairs_ring_disallowed_atoms and top_pair_ind in top_pairs_ring_disallowed_atoms:
                temp_ring_disallowed_atoms = top_pairs_ring_disallowed_atoms[top_pair_ind]
            else:
                temp_ring_disallowed_atoms = None
            for neigh_pair_level in self.pyramid_neighbours: # loop over the pyramid
            # 1st neighbours, 1st & 2nd, 2nd & 2nd, 2nd & 3rd, 3rd & 3rd, ... till level + level == nb
                top_pair_neigh_levels = self.neigh_levels[list(top_pair_ind)]
                for temp_sol_pair in sol._sol.values[:,top_pair_ind]:
                    if None in temp_sol_pair or Dummy in temp_sol_pair: # skip as None and Dummy don't have neighbours
                        continue
                    top_i_j_top_atom = ([], [])
                    temp_top_pair_atoms = ([], [])
                    for top_i_j in range(2):
                        if neigh_pair_level[top_i_j] in top_pair_neigh_levels[top_i_j][temp_sol_pair[top_i_j]]:
                            for temp_at in top_pair_neigh_levels[top_i_j][temp_sol_pair[top_i_j]][neigh_pair_level[top_i_j]]:
                                if top_pairs_ring_allowed_atoms and \
                                        temp_at not in top_pairs_ring_allowed_atoms[top_pair_ind][top_i_j]:
                                    continue # only allowed atoms to be considered
                                if temp_ring_disallowed_atoms and temp_at in temp_ring_disallowed_atoms[top_i_j]:
                                    continue # this is a disallowed atom
                                temp_top_i_j_atom_i_j = (top_pair_ind[top_i_j], temp_at) # make into top_index_atom
                                if temp_top_i_j_atom_i_j in sol.sol_atoms_map:
                                    temp_top_i_j_atom_i_j = (None, sol.sol_atoms_map[temp_top_i_j_atom_i_j]) # this is row identifier
                                elif temp_at not in sol.available_atoms[top_pair_ind[top_i_j]]: # check if in available
                                    continue
                                if temp_at not in temp_top_pair_atoms[top_i_j]:
                                    top_i_j_top_atom[top_i_j].append(temp_top_i_j_atom_i_j)
                                    temp_top_pair_atoms[top_i_j].append(temp_at)
                        #top_pair_atoms[top_pair_ind][top_i_j] |= set(temp_top_pair_atoms[top_i_j])
                    for temp_idx_atom0, temp_top_i_atom0 in enumerate(top_i_j_top_atom[0]):
                        for temp_idx_atom1, temp_top_i_atom1 in enumerate(top_i_j_top_atom[1]):
                            if temp_top_i_atom1[0] == None:
                                temp_pair = (temp_top_i_atom1, temp_top_i_atom0)
                                if temp_top_i_atom0[0] is None:# both rows!
                                    temp_pair = tuple(sorted(temp_pair))
                            else:
                                temp_pair = (temp_top_i_atom0, temp_top_i_atom1)
                            if self.check_match_pair(temp_pair, sol, **kwargs):
                                temp_at0 = temp_top_pair_atoms[0][temp_idx_atom0]
                                temp_at1 =  temp_top_pair_atoms[1][temp_idx_atom1]
                                self.calculate_estimate_atom_pair_in_sol(sol, top_pair_ind, temp_at0, temp_at1)
                                yield temp_pair

    def __get_unsorted_new_pairs_stepwise_2_part_old(self, sol, **kwargs):
        # new atom-atom pairs from atom-atom pairs in sol -> new rows as well as replacing None in existing rows
        top_pairs_ring_allowed_atoms = kwargs.get('top_pairs_ring_allowed_atoms')
        top_pairs_ring_disallowed_atoms = kwargs.get('top_pairs_ring_disallowed_atoms')
        tops_in_sol = tuple(top_i for top_i in sol._sol.columns if not sol.flag_top_not_in_sol[top_i])
        for top_pair_ind in combinations(tops_in_sol, 2):
            if top_pairs_ring_allowed_atoms and top_pair_ind not in top_pairs_ring_allowed_atoms:
                assert top_pair_ind[::-1] not in top_pairs_ring_allowed_atoms
                continue # only allowed top pairs to be considered
            if top_pairs_ring_disallowed_atoms and top_pair_ind in top_pairs_ring_disallowed_atoms:
                temp_ring_disallowed_atoms = top_pairs_ring_disallowed_atoms[top_pair_ind]
            else:
                temp_ring_disallowed_atoms = None
            for neigh_pair_level in self.pyramid_neighbours: # loop over the pyramid
            # 1st neighbours, 1st & 2nd, 2nd & 2nd, 2nd & 3rd, 3rd & 3rd, ... till level + level == nb
                top_pair_neigh_levels = self.neigh_levels[list(top_pair_ind)]
                for temp_sol_pair in sol._sol.values[:,top_pair_ind]:
                    if None in temp_sol_pair or Dummy in temp_sol_pair: # skip as None and Dummy don't have neighbours
                        continue
                    if sum([neigh_pair_level[i] in top_pair_neigh_levels[i][temp_sol_pair[i]] for i in range(2)]) != 2:
# skip as in a small molecule (or a small subproblem), not each atom has to have nth level of neighbours
                        continue
                    for temp_at0 in top_pair_neigh_levels[0][temp_sol_pair[0]][neigh_pair_level[0]]:
                        if top_pairs_ring_allowed_atoms and temp_at0 not in top_pairs_ring_allowed_atoms[top_pair_ind][0]:
                            continue # only allowed atoms to be considered
                        if temp_ring_disallowed_atoms and temp_at0 in temp_ring_disallowed_atoms[0]:
                            continue # this is a disallowed atom
                        temp_top_i_atom0 = (top_pair_ind[0], temp_at0) # make into top_index_atom
                        if temp_top_i_atom0 in sol.sol_atoms_map:
                            temp_top_i_atom0 = (None, sol.sol_atoms_map[temp_top_i_atom0]) # this is row identifier
                        elif temp_at0 not in sol.available_atoms[top_pair_ind[0]]: # check if in available
                            continue
                        for temp_at1 in top_pair_neigh_levels[1][temp_sol_pair[1]][neigh_pair_level[1]]:
                            if top_pairs_ring_allowed_atoms and temp_at1 not in top_pairs_ring_allowed_atoms[top_pair_ind][1]:
                                continue # only allowed atoms to be considered
                            if temp_ring_disallowed_atoms and temp_at1 in temp_ring_disallowed_atoms[1]:
                                continue # this is a disallowed atom
                            temp_top_i_atom1 = (top_pair_ind[1], temp_at1)
                            if temp_top_i_atom1 in sol.sol_atoms_map:# already a row
                                temp_top_i_atom1 = (None, sol.sol_atoms_map[temp_top_i_atom1])
                                temp_pair = (temp_top_i_atom1, temp_top_i_atom0)
                                if temp_top_i_atom0[0] is None:# both rows!
                                    temp_pair = tuple(sorted(temp_pair))
                            elif temp_at1 not in sol.available_atoms[top_pair_ind[1]]: # check if in available
                                continue
                            else:
                                temp_pair = (temp_top_i_atom0, temp_top_i_atom1)
                            if self.check_match_pair(temp_pair, sol, **kwargs):
                                self.calculate_estimate_atom_pair_in_sol(sol, top_pair_ind, temp_at0, temp_at1)
                                yield temp_pair

    def __check_atoms_in_sol_avail(self, sol, top_pair_ind, atoms_0_1):
        for i in range(2):
            for temp_avail_atom in atoms_0_1[i]:
                if temp_avail_atom not in sol.available_atoms[top_pair_ind[i]]:
                    if (top_pair_ind[i], temp_avail_atom) not in sol.sol_atoms_map:
                        return False
                    row = sol.sol_atoms_map[(top_pair_ind[i], temp_avail_atom)]
                    #if sol._sol.values[row, top_pair_ind[(i + 1) % 2]] != None:
                    if sol._sol.loc[row, top_pair_ind[(i + 1) % 2]] != None:
                        return False
        return True

    def __get_unsorted_new_pairs_stepwise_2_part_ring_ind_parts(self, sol, tops_ind_ring_parts_pairs,
                                                                flag_valid_closure, **kwargs):
        for tops_ind_ring_parts_pair in tops_ind_ring_parts_pairs:
            temp_flag_valid_closure = False
            top_pair_ind, ring_ind_part_in_sol_0_1, ring_ind_part_avail_0_1 = tops_ind_ring_parts_pair
            flag_atoms_in_sol_avail = self.__check_atoms_in_sol_avail(sol, top_pair_ind, ring_ind_part_avail_0_1)
            if flag_atoms_in_sol_avail == False:
                flag_valid_closure[0] = False
                return
            for temp_sol_pair in zip(*ring_ind_part_in_sol_0_1):
                for temp_at0 in self.tops[top_pair_ind[0]].adj[temp_sol_pair[0]]:
                    if temp_at0 not in ring_ind_part_avail_0_1[0]:
                        continue
                    temp_top_i_atom0 = (top_pair_ind[0], temp_at0) # make into top_index_atom
                    if temp_top_i_atom0 in sol.sol_atoms_map:
                        temp_top_i_atom0 = (None, sol.sol_atoms_map[temp_top_i_atom0]) # this is row identifier
                    for temp_at1 in self.tops[top_pair_ind[1]].adj[temp_sol_pair[1]]:
                        if temp_at1 not in ring_ind_part_avail_0_1[1]:
                            continue
                        temp_top_i_atom1 = (top_pair_ind[1], temp_at1)
                        if temp_top_i_atom1 in sol.sol_atoms_map:# already a row
                            temp_top_i_atom1 = (None, sol.sol_atoms_map[temp_top_i_atom1])
                            temp_pair = (temp_top_i_atom1, temp_top_i_atom0)
                            if temp_top_i_atom0[0] is None:# both rows!
                                temp_pair = tuple(sorted(temp_pair))
                        else:
                            temp_pair = (temp_top_i_atom0, temp_top_i_atom1)
                        if self.check_match_pair(temp_pair, sol, **kwargs):
                            RIP_in_sol_0_1_sets = tuple(frozenset(temp_list) for temp_list in ring_ind_part_in_sol_0_1)
                            temp_ind_ring_estimate = \
                                self.calculate_estimate_atom_pair_in_sol_ring_ind(sol, top_pair_ind,
                                temp_at0, temp_at1, RIP_in_sol_0_1_sets, ring_ind_part_avail_0_1)
                            if temp_ind_ring_estimate != len(ring_ind_part_avail_0_1[0]):
                                continue
                            self.calculate_estimate_atom_pair_in_sol(sol, top_pair_ind, temp_at0, temp_at1)
                            temp_flag_valid_closure = True
                            yield temp_pair
            if not temp_flag_valid_closure:
                flag_valid_closure[0] = False
                return

    def __get_unsorted_new_pairs_stepwise_3_part_top_top_gen(self, sol, **kwargs):
        pass

    def __get_unsorted_new_pairs_stepwise_3_part(self, sol, **kwargs):
        # new row! - first neighbours from atoms in sol + (all atoms from tops not in sol) or (Dummy for each other top)
        for top_i, top_i_sol_atoms in enumerate(sol._sol.values.T):
            if sol.flag_top_not_in_sol[top_i]:
                continue
            top_new_atoms = []
            for top_i_sol_at in top_i_sol_atoms:
                if top_i_sol_at not in (None, Dummy):
                    for new_at in self.neigh_levels[top_i][top_i_sol_at][1]:
                        if new_at in sol.available_atoms[top_i] and new_at not in top_new_atoms:
                            top_new_atoms.append(new_at)
                            for top_j in  sol._sol.columns:
                                if top_i == top_j:
                                    continue
                                temp_pair_group = []
                                if top_i < top_j:
                                    fac = 1
                                    top_pair_ind = (top_i, top_j)
                                else:
                                    top_pair_ind = (top_j, top_i)
                                    fac = -1
                                if sol.flag_top_not_in_sol[top_j]:
                                    for temp_atom in sol.available_atoms[top_j]:
                                        temp_pair = ((top_i, new_at), (top_j, temp_atom))[::fac]
                                        if self.check_match_pair(temp_pair, sol, **kwargs):
                                            temp_pair_group.append(temp_pair)
                                            #yield temp_pair
                                    temp_pair = (top_i, new_at), (top_j, Dummy)
                                else:
                                    temp_pair = None
                                    tgp = sol.TGPs_tops_in_sol[top_pair_ind][0]
                                    if tgp._in_solution_Vs:
                                        temp_pair = (top_i, new_at), (top_j, Dummy)
                                if temp_pair and self.check_match_pair(temp_pair, sol, **kwargs):
                                    temp_pair_group.append(temp_pair)
                                        #yield temp_pair
                                yield temp_pair_group
#                            yield (top_i, new_at), (None, None)

    @staticmethod
    def __get_top_atom_pair_from_temp_pair(temp_pair):
        top_i, atom_i = temp_pair[0]
        top_j, atom_j = temp_pair[1]
        if top_i > top_j:
            top_i, atom_i, top_j, atom_j = top_j, atom_j, top_i, atom_i
        return top_i, atom_i, top_j, atom_j

    def __get_estimate_top_atom_pair_from_temp_pair(self, temp_pair, sol):
        top_i, atom_i, top_j, atom_j = self.__get_top_atom_pair_from_temp_pair(temp_pair)
        if atom_i == Dummy:
            atom_ii = -1
        else:
            atom_ii = sol.pair_estimates_available_atoms_map[(top_i, top_j)][0][atom_i]
        if atom_j == Dummy:
            atom_jj = -1
        else:
            atom_jj = sol.pair_estimates_available_atoms_map[(top_i, top_j)][1][atom_j]
        return top_i, atom_ii, top_j, atom_jj

    def __get_row_atom_estimate(self, temp_pair, sol, **kwargs):
        num_non_dummies = [0,0]
        non_neigh_estimate = [-1]
        top_j, atom_j = temp_pair[1]
        temp_estimates = []
        tops2enumerate = kwargs.get('tops2enumerate')
        if tops2enumerate:
            temp_gen = enumerate(sol._sol.loc[temp_pair[0][1], tops2enumerate])
        else:
            temp_gen = enumerate(sol._sol.loc[temp_pair[0][1], :])
        for top_i, atom_i in temp_gen:
            if tops2enumerate:
                top_i = tops2enumerate[top_i]
            if atom_i is not None:
                if top_i < top_j:
                    if atom_i == Dummy:
                        atom_ii = -1
                    else:
                        num_non_dummies[0]+=1
                        if atom_i not in sol.pair_estimates_available_atoms_map[(top_i, top_j)][0]:
                            return non_neigh_estimate, num_non_dummies
                        atom_ii = sol.pair_estimates_available_atoms_map[(top_i, top_j)][0][atom_i]
                    if atom_j == Dummy:
                        atom_jj = -1
                    else:
                        if atom_j not in sol.pair_estimates_available_atoms_map[(top_i, top_j)][1]:
                            return non_neigh_estimate, num_non_dummies
                        atom_jj = sol.pair_estimates_available_atoms_map[(top_i, top_j)][1][atom_j]
                    temp_estimates.append(sol.pair_estimates[(top_i, top_j)][atom_ii, atom_jj])
                else:
                    if atom_i == Dummy:
                        atom_ii = -1
                    else:
                        num_non_dummies[1]+=1
                        if atom_i not in sol.pair_estimates_available_atoms_map[(top_j, top_i)][1]:
                            return non_neigh_estimate, num_non_dummies
                        atom_ii = sol.pair_estimates_available_atoms_map[(top_j, top_i)][1][atom_i]
                    if atom_j == Dummy:
                        atom_jj = -1
                    else:
                        if atom_j not in sol.pair_estimates_available_atoms_map[(top_j, top_i)][0]:
                            return non_neigh_estimate, num_non_dummies
                        atom_jj = sol.pair_estimates_available_atoms_map[(top_j, top_i)][0][atom_j]
                    temp_estimates.append(sol.pair_estimates[(top_j, top_i)][atom_jj, atom_ii])
        return sorted(temp_estimates, reverse=True), num_non_dummies

    def __get_row_row_atom_estimate(self, temp_pair, fake_pair, row, sol, **kwargs):
        tops2enumerate = []
        for top_i in sol._sol.columns:
            temp_atom = sol._sol.loc[row, top_i]
            if temp_atom == None:
                tops2enumerate.append(top_i)
        return self.__get_row_atom_estimate(fake_pair, sol, tops2enumerate = tops2enumerate, **kwargs)

    def __compare_best_atom_estimate_add2row(self, current_best, temp_pair_estimate, num_non_dummies, dummy_flag):
        # current best is (num_non_dummies, estimates, dummy_flag)
        if current_best is None:return True
        if num_non_dummies > current_best[0]:
            return True
        if num_non_dummies < current_best[0]:
            return False
        for i, single_estimate in enumerate(temp_pair_estimate):
            try:
                if single_estimate > current_best[1][i]:
                    return True
                elif single_estimate < current_best[1][i]:
                    return False
            except:
                if single_estimate > 0:
                    return True
        if current_best[2] and not dummy_flag:
            return True
        return False

    def __generate_row_atom_pairs_from_row_row_pair(self, sol, temp_pair):
        row_tops_in_sol_0_1 = [], []
        row_tops_none_0_1 = [], []
        for top_i in sol._sol.columns:
            for i in range(2):
                if sol._sol.loc[temp_pair[i][1], top_i] not in (None, Dummy):
                    row_tops_in_sol_0_1[i].append(top_i)
                if sol._sol.loc[temp_pair[i][1], top_i] == None:
                    row_tops_none_0_1[i].append(top_i)
        for i in range(2):
            other_i = (i+1) % 2
            for top_i_none in row_tops_none_0_1[i]:
                if top_i_none in row_tops_in_sol_0_1[other_i]:
                    temp_row_top = temp_pair[i][1], top_i_none
                    row_other = temp_pair[other_i][1]
                    temp_atom = sol._sol.loc[row_other, top_i_none]
                    yield temp_row_top, temp_atom, row_other

    def get_sorted_new_pairs_stepwise(self, sol, **kwargs):
        #self.make_estimates(sol)
#        flag_find_solution, tops_ind_ring_parts_pairs = self.__find_partial_ind_ring_parts_in_sol(sol)
        #flag_find_solution, check_res = self.__find_rings2check_in_sol(sol)
        flag_find_solution, check_res = sol.flag_find_solution, sol.check_res
        if not flag_find_solution:
            # if some of the ring pairs does not have a valid solution for closure...
            return
        flag_find_solution = [flag_find_solution]
        flag_partial_valid_solution, tops_ind_ring_parts_pairs, top_pairs_allowed_atoms, top_pairs_ring_disallowed_atoms \
            = check_res
        best_pair = None
        best_atom_estimate = None
        additional_add2row = []
        row_top_all_atoms = {}
        new_row_top_atom_all_atoms = {}
        if tops_ind_ring_parts_pairs:
            temp_pair_generator = \
                self.__get_unsorted_new_pairs_stepwise_2_part_ring_ind_parts(sol, tops_ind_ring_parts_pairs,
                                                                             flag_find_solution, **kwargs)
        elif flag_partial_valid_solution == False:
            new_kwargs = {'top_pairs_ring_allowed_atoms':top_pairs_allowed_atoms}
            new_kwargs.update(kwargs)
            temp_pair_generator = self.__get_unsorted_new_pairs_stepwise_2_part(sol, **new_kwargs)
            #temp_pair_generator = self.__get_unsorted_new_pairs_stepwise_2_part_old(sol, **new_kwargs)
        else:
            kwargs['top_pairs_ring_disallowed_atoms'] = top_pairs_ring_disallowed_atoms
            temp_pair_generator = self.__get_unsorted_new_pairs_stepwise_2_part(sol, **kwargs)
            #temp_pair_generator = self.__get_unsorted_new_pairs_stepwise_2_part_old(sol, **kwargs)
        for temp_pair in temp_pair_generator:
            if temp_pair[0][0] is None:
                if temp_pair not in additional_add2row: # as it might yield duplicates
                    additional_add2row.append(temp_pair)
                    if temp_pair[1][0] is None:
                        for temp_row_top, temp_atom, row_other in \
                                self.__generate_row_atom_pairs_from_row_row_pair(sol, temp_pair):
                            if temp_row_top not in row_top_all_atoms:
                                row_top_all_atoms[temp_row_top] = list()
                            if temp_atom not in row_top_all_atoms[temp_row_top]:
                                row_top_all_atoms[temp_row_top].append(temp_atom)
                    else:
                        temp_row_top = temp_pair[0][1], temp_pair[1][0]
                        if temp_row_top not in row_top_all_atoms:
                            row_top_all_atoms[temp_row_top] = list()
                        if temp_pair[1][1] not in row_top_all_atoms[temp_row_top]:
                            row_top_all_atoms[temp_row_top].append(temp_pair[1][1])
            else:
                top_atom_top = temp_pair[0][0], temp_pair[0][1], temp_pair[1][0]
                if top_atom_top not in new_row_top_atom_all_atoms:
                    new_row_top_atom_all_atoms[top_atom_top] = set()
                new_row_top_atom_all_atoms[top_atom_top].add(temp_pair[1][1])
        if not flag_find_solution[0]:
            # if some of the ind ring part pairs does not have a valid solution for closure...
            return

        for top_pair_ind in tuple(sol.TGPs_tops_in_sol):
            self.calculate_estimate_atom_dummy_pairs_in_sol(sol, top_pair_ind) # sorts out estimates to dummy atoms
        #if not tops_ind_ring_parts_pairs:
        if flag_partial_valid_solution:
            for temp_pair in self.__get_unsorted_new_pairs_stepwise_1_part(sol, **kwargs):
                temp_row_top = temp_pair[0][1], temp_pair[1][0]
                if temp_row_top not in row_top_all_atoms:
                    row_top_all_atoms[temp_row_top] = list()
#                elif isinstance(row_top_all_atoms[temp_row_top], set):
#                    row_top_all_atoms[temp_row_top] = list(row_top_all_atoms[temp_row_top])
                assert temp_pair[1][1] not in row_top_all_atoms[temp_row_top]
                row_top_all_atoms[temp_row_top].append(temp_pair[1][1])
                temp_estimate, num_non_dummies = self.__get_row_atom_estimate(temp_pair, sol, **kwargs)
                dummy_flag = temp_pair[1][1] == Dummy
                if self.__compare_best_atom_estimate_add2row(best_atom_estimate, temp_estimate, num_non_dummies, dummy_flag):
                    best_atom_estimate = num_non_dummies, temp_estimate, dummy_flag
                    best_pair = temp_pair
        for temp_pair in additional_add2row:
            if temp_pair[1][0] != None:
                temp_estimate, num_non_dummies = self.__get_row_atom_estimate(temp_pair, sol, **kwargs)
                dummy_flag = temp_pair[1][1] == Dummy
                if self.__compare_best_atom_estimate_add2row(best_atom_estimate, temp_estimate, num_non_dummies, dummy_flag):
                    best_atom_estimate = num_non_dummies, temp_estimate, dummy_flag
                    best_pair = temp_pair
            else:
                for temp_row_top, temp_atom, row in self.__generate_row_atom_pairs_from_row_row_pair(sol, temp_pair):
                    temp_fake_pair = (None, temp_row_top[0]), (temp_row_top[1], temp_atom)
                    temp_estimate, num_non_dummies = self.__get_row_row_atom_estimate(temp_pair, temp_fake_pair, row, sol, **kwargs)
#                    temp_estimate = self.__get_row_atom_estimate(temp_pair, sol, **kwargs)
                    dummy_flag = False
                    assert temp_pair[1][1] != Dummy
                    if self.__compare_best_atom_estimate_add2row(best_atom_estimate, temp_estimate, num_non_dummies, dummy_flag):
                        best_atom_estimate = num_non_dummies, temp_estimate, dummy_flag
#                        best_pair = temp_pair
                        best_pair = temp_fake_pair
        if best_pair:
            pairs_from_best_row_top = {}
            best_row_top = best_pair[0][1], best_pair[1][0]
            check_row_row = sol.transform_top_index_atom(best_pair[1])
            ############## DEBUGGING ######################################
            if check_row_row[0] == None:
                true_best_pair = tuple(sorted((best_pair[0], check_row_row)))
            else:
                true_best_pair = best_pair
            if true_best_pair[1][0] == None:
                pass
            ############## DEBUGGING ######################################
            for temp_at_j in row_top_all_atoms[best_row_top]:
                check_row_row = sol.transform_top_index_atom((best_pair[1][0], temp_at_j))
                if check_row_row[0] == None:
                    temp_pair = tuple(sorted((best_pair[0], check_row_row)))
                    temp_fake_pair = best_pair[0], (best_pair[1][0], temp_at_j)
                    row = check_row_row[1]
                    pairs_from_best_row_top[temp_pair], num_non_dummies =\
                        self.__get_row_row_atom_estimate(temp_pair, temp_fake_pair, row, sol, **kwargs)
                else:
                    temp_pair = (best_pair[0], check_row_row)
                    pairs_from_best_row_top[temp_pair], num_non_dummies = self.__get_row_atom_estimate(temp_pair, sol, **kwargs)
                pairs_from_best_row_top[temp_pair].append(temp_pair[1][1] != Dummy)
            for temp_pair, _ in sorted(pairs_from_best_row_top.items(), key=lambda item: item[1], reverse=True):
                yield temp_pair
            return
        best_atom_estimate = -1
        best_atom_estimate_dummy_flag = True
        best_atom_estimate_group = None
        for temp_top_atom_top in new_row_top_atom_all_atoms:
            for temp_atom in new_row_top_atom_all_atoms[temp_top_atom_top]:
                temp_pair = (temp_top_atom_top[:2], (temp_top_atom_top[2], temp_atom))
                top_i, atom_ii, top_j, atom_jj = self.__get_estimate_top_atom_pair_from_temp_pair(temp_pair, sol)
                temp_estimate = sol.pair_estimates[(top_i, top_j)][atom_ii, atom_jj]
                if temp_estimate > best_atom_estimate:
                    best_pair = temp_pair
                    best_atom_estimate = temp_estimate
                    best_atom_estimate_dummy_flag = False
                    best_atom_estimate_group = temp_top_atom_top
        best_top_atom_top_estimate_group = None
        if best_atom_estimate_group:
            best_top_atom_top_estimate_group = best_atom_estimate_group
            temp_pair_group = []
            for temp_atom in new_row_top_atom_all_atoms[best_atom_estimate_group]:
                temp_pair = (best_atom_estimate_group[:2], (best_atom_estimate_group[2], temp_atom))
                temp_pair_group.append(temp_pair)
            best_atom_estimate_group = temp_pair_group
        #if not tops_ind_ring_parts_pairs:
        if flag_partial_valid_solution:
            for temp_pair_group in self.__get_unsorted_new_pairs_stepwise_3_part(sol, **kwargs):
                if temp_pair_group:
                    temp_pair = temp_pair_group[0]
                    top_atom_top = temp_pair[0][0], temp_pair[0][1], temp_pair[1][0]
                    if top_atom_top in new_row_top_atom_all_atoms:
                        if top_atom_top == best_top_atom_top_estimate_group:
                            best_atom_estimate_group.append(temp_pair)
                        continue
                for temp_pair in temp_pair_group:
                    top_i, atom_ii, top_j, atom_jj = self.__get_estimate_top_atom_pair_from_temp_pair(temp_pair, sol)
                    temp_estimate = sol.pair_estimates[(top_i, top_j)][atom_ii, atom_jj]
                    if temp_estimate > best_atom_estimate:
                        best_pair = temp_pair
                        best_atom_estimate = temp_estimate
                        best_atom_estimate_dummy_flag = temp_pair[1][1] == Dummy
                        best_atom_estimate_group = temp_pair_group
                        best_top_atom_top_estimate_group = top_atom_top
                    elif temp_estimate == best_atom_estimate and best_atom_estimate_dummy_flag and temp_pair[1][1] != Dummy:
                        best_pair = temp_pair
                        best_atom_estimate_dummy_flag = False
                        best_atom_estimate_group = temp_pair_group
                        best_top_atom_top_estimate_group = top_atom_top
        if best_pair:
            pairs_from_best_pair = {}
            for temp_pair in best_atom_estimate_group:
                top_i, atom_ii, top_j, atom_jj = self.__get_estimate_top_atom_pair_from_temp_pair(temp_pair, sol)
                temp_estimate = sol.pair_estimates[(top_i, top_j)][atom_ii, atom_jj]
                pairs_from_best_pair[temp_pair] = (temp_estimate, temp_pair[1][1] != Dummy)
            for temp_pair, _ in sorted(pairs_from_best_pair.items(), key=lambda item: item[1], reverse=True):
                yield temp_pair

    ####################################################################################################################
    #                                                   END
    #                                    stepwise approach with sorting
    ####################################################################################################################

    ####################################################################################################################
    #                                        match estimates for atoms and dummies
    ####################################################################################################################
    def __get_initial_estimate_matrices(self, **kwargs): # just for testing!!!!!!
        sol = self.initial_sol
        self.initial_available_atoms = [tuple(temp_avail_atoms) for temp_avail_atoms in sol.available_atoms]
        initial_available_atoms_frozenset = [frozenset(temp_avail_atoms) for temp_avail_atoms in sol.available_atoms]
        self.initial_pair_estimates = {}
        self.initial_pair_estimates_available_atoms = {}
        for top_pair_ind in combinations(self.initial_sol._sol.columns, 2):
            available_atoms = tuple(self.initial_available_atoms[top_i] for top_i in top_pair_ind)
            available_atoms_frozenset = tuple(initial_available_atoms_frozenset[top_i] for top_i in top_pair_ind)
            temp_pair_estimates = self.get_all_vs_all_estimate_matrix_top_pair(top_pair_ind, available_atoms,
                                                                               available_atoms_frozenset)
            self.initial_pair_estimates[top_pair_ind] = temp_pair_estimates
            self.initial_pair_estimates_available_atoms[top_pair_ind] = available_atoms

    def make_estimates(self, sol, **kwargs):
        sol.pair_estimates = {}
        sol.pair_estimates_available_atoms = {}
        sol.pair_estimates_available_atoms_map = {}
        sol.TGPs_tops_in_sol = {}
        for top_pair_ind in combinations(self.initial_sol._sol.columns, 2):
            # both topologies not in solution
            if sol.flag_top_not_in_sol[top_pair_ind[0]] and sol.flag_top_not_in_sol[top_pair_ind[1]]:
                available_atoms = tuple(sol.available_atoms[top_i] for top_i in top_pair_ind)
                available_atoms_frozenset = tuple(frozenset(sol.available_atoms[top_i]) for top_i in top_pair_ind)
                temp_pair_estimates = self.get_all_vs_all_estimate_matrix_top_pair(top_pair_ind, available_atoms,
                                                                                   available_atoms_frozenset)
                sol.pair_estimates[top_pair_ind] = temp_pair_estimates
                sol.pair_estimates_available_atoms[top_pair_ind] = available_atoms
            # both topologies already in solution
            elif not sol.flag_top_not_in_sol[top_pair_ind[0]] and not sol.flag_top_not_in_sol[top_pair_ind[1]]:
                # find available and in sol atoms
                # available atoms are all not in sol + these that are at the time matched with None
                available_atoms = tuple(list(sol.available_atoms[top_i]) for top_i in top_pair_ind)
                in_solution_atoms = ([], [])
                in_solution_dummy_match = ([], [])
                for temp_atom_pair in sol._sol.values[:, top_pair_ind]:
                    """
                    if temp_atom_pair[0] not in (None, Dummy) and temp_atom_pair[1] not in (None, Dummy):
                        in_solution_atoms[0].append(temp_atom_pair[0])
                        in_solution_atoms[1].append(temp_atom_pair[1])
                        continue
                    """
                    for i in range(2):
                        ii = i
                        jj = (i + 1) % 2
                        if temp_atom_pair[ii] not in (None, Dummy):
                            if temp_atom_pair[jj] is None:
                                available_atoms[ii].append(temp_atom_pair[ii])
                            else:
                                if temp_atom_pair[jj] == Dummy:
                                    in_solution_dummy_match[ii].append(temp_atom_pair[ii])
                                    #in_solution_atoms[ii].append(temp_atom_pair[ii])
                                else:
                                    in_solution_atoms[ii].append(temp_atom_pair[ii])
                available_atoms = (frozenset(available_atoms[0]), frozenset(available_atoms[1]))
                in_solution_atoms = (frozenset(in_solution_atoms[0]), frozenset(in_solution_atoms[1]))
                in_solution_dummy_match = tuple(frozenset(in_solution_dummy_match[i]) for i in range(2))
                sol.pair_estimates[top_pair_ind], sol.pair_estimates_available_atoms[top_pair_ind],\
                    sol.TGPs_tops_in_sol[top_pair_ind] = self.create_estimate_top_pair_in_sol\
                    (top_pair_ind, available_atoms, in_solution_atoms, in_solution_dummy_match)
                # here, sol.pair_estimates_available_atoms is actually Vs_of_interest
            # one of the topologies in the solution and one not
            else:
                # finding which one is in the solution and which one not
                if sol.flag_top_not_in_sol[top_pair_ind[0]]:
                    available_atoms_top_in_sol = list(sol.available_atoms[top_pair_ind[1]])
                    available_atoms = (tuple(sol.available_atoms[top_pair_ind[0]]), available_atoms_top_in_sol)
                    top_in_sol_i = 1
                    top_not_in_sol_i = 0
                else:
                    available_atoms_top_in_sol = list(sol.available_atoms[top_pair_ind[0]])
                    available_atoms = (available_atoms_top_in_sol, tuple(sol.available_atoms[top_pair_ind[1]]))
                    top_in_sol_i = 0
                    top_not_in_sol_i = 1
                # find which atoms are in the solution that are still None for the other topology
                for temp_atom_pair in sol._sol.values[:, top_pair_ind]:
                    if temp_atom_pair[top_in_sol_i] not in (None, Dummy) and temp_atom_pair[top_not_in_sol_i] is None:
                        available_atoms_top_in_sol.append(temp_atom_pair[top_in_sol_i])
                available_atoms_frozenset = tuple(frozenset(available_atoms[i]) for i in range(2))
                temp_pair_estimates = self.get_all_vs_all_estimate_matrix_top_pair(top_pair_ind, available_atoms,
                                                                                   available_atoms_frozenset)
                sol.pair_estimates[top_pair_ind] = temp_pair_estimates
                sol.pair_estimates_available_atoms[top_pair_ind] = available_atoms
            available_atoms_map = tuple(OrderedDict(zip(av_atoms_tup, range(len(av_atoms_tup))))
                                   for av_atoms_tup in sol.pair_estimates_available_atoms[top_pair_ind])
            sol.pair_estimates_available_atoms_map[top_pair_ind] = available_atoms_map
        self._make_full_sol_estimate(sol, **kwargs)

    def _make_full_sol_estimate(self, sol, **kwargs):
        full_estimate = []
        flag_find_solution, check_res = self.__find_rings2check_in_sol(sol)
        sol.flag_find_solution, sol.check_res = flag_find_solution, check_res
        if not flag_find_solution:
            sol.full_estimate = 0
            return 0
        for top_pair_ind in self._top_pairs:
            if top_pair_ind not in sol.TGPs_tops_in_sol:
                temp_estimate = sol.pair_estimates[top_pair_ind].max()
                full_estimate.append(temp_estimate)
            else:
                TGPs = sol.TGPs_tops_in_sol[top_pair_ind]
                if not TGPs[0]._in_solution_Vs:
                    temp_estimate = sol.pair_estimates[top_pair_ind].max()
                    full_estimate.append(temp_estimate)
                    continue
                temp_top_pair_estimate = 0
                done_Vs_in_sol = [[], []]
                done_shared_Vs_in_sol = [[], []]
                for ii in range(2):
                    jj = (ii + 1) % 2
                    for common_Vs_in_sol_ii in TGPs[ii].shared_Vs_in_sol:
                        if common_Vs_in_sol_ii in done_shared_Vs_in_sol[ii]:
                            continue
                        common_Vs_in_sol_jj = []
                        at_ii_jj_map = {}
                        for at_ii in common_Vs_in_sol_ii:
                            at_jj = sol.find_common_atom_pair((top_pair_ind[ii], at_ii), top_pair_ind[jj])
                            #row = sol.sol_atoms_map[(top_pair_ind[ii], at_ii)]
                            #at_jj = sol._sol.loc[row, top_pair_ind[jj]]
                            at_ii_jj_map[at_ii] = at_jj
                            assert at_jj not in (None, Dummy)
                            common_Vs_in_sol_jj.append(at_jj)
                        common_Vs_in_sol_jj = frozenset(common_Vs_in_sol_jj)
                        if common_Vs_in_sol_jj in TGPs[jj].groups4estimates:
                            common_Vs_in_sol_ii_jj = common_Vs_in_sol_ii, common_Vs_in_sol_jj
                            temp_estimates = []
                            for at_ii in common_Vs_in_sol_ii:
                                comb_branches_ii_jj = []
                                comb_branches_ii_jj_neigh_count_map = []
                                at_jj = at_ii_jj_map[at_ii]
                                for i in range(2):
                                    ij = (ii, jj)[i]
                                    atom_ij = (at_ii, at_jj)[i]
                                    temp_b = TGPs[ij].combined_branches[atom_ij]
                                    all_parents = set(temp_b.parents)
                                    parents = all_parents & TGPs[ij].groups4estimates[common_Vs_in_sol_ii_jj[i]]
                                    if all_parents != parents:
                                        temp_b = \
                                            TGPs[ij]._get_combined_neigh_count_Nv_allowed_parents(atom_ij, parents)
                                    comb_branches_ii_jj.append(temp_b)
                                    comb_branches_ii_jj_neigh_count_map.append(self.__get_neigh_count_from_neigh_map
                                                                    (top_pair_ind[ij], atom_ij, temp_b.V_in_branch))
                                temp_estimate = self.calculate_combined_branch_pair_estimate(*comb_branches_ii_jj) - 1
                                temp_estimates.append(temp_estimate)
                                temp_estimate = temp_b.calculate_neigh_count_pair_estimate(comb_branches_ii_jj_neigh_count_map)
                                temp_estimates.append(temp_estimate)
                            # roots parents together branches estimate
                            """
                            old_comb_branches_ii_jj = list(comb_branches_ii_jj)
                            available_Vs_ij = [comb_branches_ii_jj[i].V_in_branch for i in range(2)]
                            comb_branches_ii_jj = []
                            for i in range(2):
                                ij = (ii, jj)[i]
                                roots = common_Vs_in_sol_ii_jj[i]
                                parents = TGPs[ij].groups4estimates[common_Vs_in_sol_ii_jj[i]]
                                available_Vs = available_Vs_ij[i]
                                temp_b = TGPs[ij]._g_get_branch_roots_parents\
                                    (TGPs[ij].g, roots, parents, available_Vs = available_Vs)
                                comb_branches_ii_jj.append(temp_b)
                                old_comb_branches_ii_jj[i].comb_ring_comm_parents_b = comb_branches_ii_jj[i]
                                if kwargs.get('verbose'):
                                    print('HM\t\t', roots, parents, '\t\t',old_comb_branches_ii_jj[i].parents, old_comb_branches_ii_jj[i].roots)
                            temp_estimate = comb_branches_ii_jj[0].calculate_pair_estimate(comb_branches_ii_jj[1])
                            temp_estimates.append(temp_estimate)
                            if kwargs.get('verbose'):print(temp_estimate)
                            """
                            # roots parents together branches estimate
                            temp_top_pair_estimate += min(temp_estimates)
                            done_shared_Vs_in_sol[ii].append(common_Vs_in_sol_ii)
                            done_shared_Vs_in_sol[jj].append(common_Vs_in_sol_jj)
                        for at_ii in common_Vs_in_sol_ii:
                            at_jj = at_ii_jj_map[at_ii]
                            if at_ii not in done_Vs_in_sol[ii]:
                                assert at_jj not in done_Vs_in_sol[jj]
                                done_Vs_in_sol[ii].append(at_ii)
                                done_Vs_in_sol[jj].append(at_jj)
                            else:
                                continue
                            comb_alone_branches_ii_jj = []
                            for i in range(2):
                                ij = (ii, jj)[i]
                                temp_k = frozenset({(at_ii, at_jj)[i]})
                                if temp_k in TGPs[ij].groups4estimates:
                                    parents = TGPs[ij].groups4estimates[temp_k]
                                    temp_b_alone = TGPs[ij]._get_combined_neigh_count_Nv_allowed_parents\
                                        ((at_ii, at_jj)[i], parents)
                                else:
                                    temp_b_alone = None
                                comb_alone_branches_ii_jj.append(temp_b_alone)
                            if None not in comb_alone_branches_ii_jj:
                                temp_estimate_alone = \
                                    self.calculate_combined_branch_pair_estimate(*comb_alone_branches_ii_jj) - 1
                            else:
                                temp_estimate_alone = 0
                            temp_top_pair_estimate += temp_estimate_alone
                for at_0 in TGPs[0]._in_solution_Vs:
                    if at_0 not in done_Vs_in_sol[0]:
                        at_1 = sol.find_common_atom_pair((top_pair_ind[0], at_0), top_pair_ind[1])
                        #row = sol.sol_atoms_map[(top_pair_ind[0], at_0)]
                        #at_1 = sol._sol.loc[row, top_pair_ind[1]]
                        assert at_1 not in (None, Dummy)
                        temp_b_0_1 = [TGPs[i].combined_branches[(at_0, at_1)[i]] for i in range(2)]
                        temp_estimate = self.calculate_combined_branch_pair_estimate(*temp_b_0_1)
                        temp_top_pair_estimate += temp_estimate - 1
                        assert at_1 not in done_Vs_in_sol[1]
                        done_Vs_in_sol[1].append(at_1)
                assert len(TGPs[0]._in_solution_Vs) == len(TGPs[1]._in_solution_Vs)
                assert len(done_Vs_in_sol[1]) == len(TGPs[0]._in_solution_Vs)
                full_estimate.append(temp_top_pair_estimate + len(TGPs[0]._in_solution_Vs))
        sol.full_estimate = full_estimate
        return full_estimate

    def __get_neigh_count_from_neigh_map(self, top_i, at, available_at):
        temp_neigh_count = Counter()
        for temp_at in available_at:
            temp_neigh_count[self.tops_neigh_map[top_i][at][temp_at]] += 1
        neigh_count =  np.zeros(len(temp_neigh_count), dtype = int)
        for i in range(len(temp_neigh_count)):
            neigh_count[i] = temp_neigh_count[i + 1]
        return neigh_count

    def get_all_vs_all_estimate_matrix_top_pair(self, top_pair_ind, available_atoms, available_atoms_frozenset=None):
        if available_atoms_frozenset is None:
            available_atoms_frozenset = [frozenset(temp_available_atoms) for temp_available_atoms in available_atoms]
        temp_shape = (len(available_atoms[0]) + 1, len(available_atoms[1]) + 1)
        temp_pair_estimates = np.zeros(temp_shape)
        TGPs = [self.get_TGP_top_available_atoms
                (top_pair_ind[i], available_atoms_frozenset[i], frozenset(), frozenset()) for i in range(2)]
        for at1_i, at1 in enumerate(available_atoms[0]):
            for at2_i, at2 in enumerate(available_atoms[1]):
                # at1, at2 are atom pairs from 2 different topologies
                b1, b2 = TGPs[0].combined_branches[at1], TGPs[1].combined_branches[at2]
                temp_pair_estimates[at1_i, at2_i] = self.calculate_combined_branch_pair_estimate(b1, b2)
        for at1_i, at1 in enumerate(available_atoms[0]):
            temp_pair_estimates[at1_i][-1] = min(len(available_atoms[0]) - 1, TGPs[0].dummy_estimates[at1])
        for at2_i, at2 in enumerate(available_atoms[1]):
            temp_pair_estimates[-1][at2_i] = min(len(available_atoms[1]) - 1, TGPs[1].dummy_estimates[at2])
        self.reduce_dummy_estimates_1(temp_pair_estimates) # based on best_atom_atom estimate
        return temp_pair_estimates

    def create_estimate_top_pair_in_sol(self, top_pair_ind, available_atoms_frozenset, in_solution_atoms_frozenset,
                                        in_solution_dummy_match_atoms_frozenset):
        TGPs = [self.get_TGP_top_available_atoms(top_pair_ind[i], available_atoms_frozenset[i],
                in_solution_atoms_frozenset[i], in_solution_dummy_match_atoms_frozenset[i]) for i in range(2)]
        if TGPs[0]._in_solution_Vs:
            temp_shape = (len(TGPs[0].Vs_of_interest) + 1, len(TGPs[1].Vs_of_interest) + 1)
            temp_pair_estimates = np.zeros(temp_shape)
            for i in range(2):
                if i == 0:
                    fac = 1
                else:
                    fac = -1
                for ati_ii, at_i in enumerate(TGPs[i].Vs_of_interest):
                    pos = (ati_ii, -1)[::fac]
                    temp_pair_estimates[pos] = TGPs[i].dummy_estimates[at_i]
        else:
            available_atoms = tuple(tuple(TGPs[i].available_Vs) for i in range(2))
            available_atoms_frozenset = tuple(TGPs[i].available_Vs for i in range(2))
            temp_pair_estimates = self.get_all_vs_all_estimate_matrix_top_pair(top_pair_ind, available_atoms,
                                                                               available_atoms_frozenset)
        return temp_pair_estimates, (TGPs[0].Vs_of_interest, TGPs[1].Vs_of_interest), TGPs

    def __calculate_estimate_atom_pair_in_sol(self, sol, top_pair_ind, atom_i, atom_j, TGPs):
        b1, b2 = TGPs[0].combined_branches[atom_i], TGPs[1].combined_branches[atom_j]
        return self.calculate_combined_branch_pair_estimate(b1, b2)

    def calculate_estimate_atom_pair_in_sol(self, sol, top_pair_ind, atom_i, atom_j):
        TGPs = sol.TGPs_tops_in_sol[top_pair_ind]
        at_ii = sol.pair_estimates_available_atoms_map[top_pair_ind][0][atom_i]
        at_jj = sol.pair_estimates_available_atoms_map[top_pair_ind][1][atom_j]
        sol.pair_estimates[top_pair_ind][at_ii, at_jj] = \
            self.__calculate_estimate_atom_pair_in_sol(sol, top_pair_ind, atom_i, atom_j, TGPs)
        return sol.pair_estimates[top_pair_ind][at_ii, at_jj]

    @staticmethod
    def reduce_dummy_estimates_1(pair_estimates):
        try:
            best_atom_atom_estimate = pair_estimates[:-1,:-1].max()
            for i in range(pair_estimates.shape[0] - 1):
                pair_estimates[i, -1] = min(best_atom_atom_estimate, pair_estimates[i, -1])
            for i in range(pair_estimates.shape[1] - 1):
                pair_estimates[-1, i] = min(best_atom_atom_estimate, pair_estimates[-1, i])
            pair_estimates[-1, -1] = best_atom_atom_estimate
        except:
            pair_estimates[:] = 0

    @staticmethod
    def reduce_dummy_estimates_2(pair_estimates):
        try:
            for i in range(pair_estimates.shape[0] - 1):
                pair_estimates[i, -1] = min(pair_estimates[i, :-1].max(), pair_estimates[i, -1])
            for i in range(pair_estimates.shape[1] - 1):
                pair_estimates[-1, i] = min(pair_estimates[:-1, i].max(), pair_estimates[-1, i])
            pair_estimates[-1, -1] = max(pair_estimates[-1, :-1].max(), pair_estimates[:-1, -1].max())
        except:
            pair_estimates[:] = 0
            pass

    def calculate_estimate_atom_dummy_pairs_in_sol(self, sol, top_pair_ind):
        TGPs = sol.TGPs_tops_in_sol[top_pair_ind]
        self.reduce_dummy_estimates_2(sol.pair_estimates[top_pair_ind])
        #del(sol.TGPs_tops_in_sol[top_pair_ind])

    def calculate_estimate_atom_pair_in_sol_ring_ind(self, sol, top_pair_ind, atom_i, atom_j, ind_ring_parts_in_sol,
                                                     ind_ring_parts_avail):
        TGPs = [self.get_TGP_top_available_atoms
          (top_pair_ind[i], ind_ring_parts_avail[i], ind_ring_parts_in_sol[i], frozenset()) for i in range(2)]
        return self.__calculate_estimate_atom_pair_in_sol(sol, top_pair_ind, atom_i, atom_j, TGPs)

    @staticmethod
    def calculate_combined_branch_pair_estimate(b1, b2):
        if b1.ring_linear_combined_branch is None and b2.ring_linear_combined_branch is not None:
            return b2.ring_linear_combined_branch.calculate_pair_estimate(b1) + 1
        elif b2.ring_linear_combined_branch is None and b1.ring_linear_combined_branch is not None:
            return b1.ring_linear_combined_branch.calculate_pair_estimate(b2) + 1
        else:
            return b1.calculate_pair_estimate(b2) + 1

    ####################################################################################################################
    #                                                   END
    #                                        match estimates for atoms and dummies
    ####################################################################################################################

    ####################################################################################################################
    #                                           generate atom order of PTP topology
    ####################################################################################################################

    def find_atom_order(self, sol):
        for at in sol.toptp.atoms.values():
            at.sol_id = at.id
        atoms2sort = set(sol._sol.index)
        while atoms2sort:
            pass

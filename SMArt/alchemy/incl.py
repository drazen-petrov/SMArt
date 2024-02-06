from SMArt.incl import copy, np, pd, combinations, product, defaultdict, frozendict, OrderedDict, do_warn
from SMArt.graph import Graph
# bond breaking / making (if True, max dist for bond breaking needed)
# breaking other bonded interactions (angles, dihedrals, impropers (needed if bond breaks = True)
    # interesting for matching dihedrals in gromacs with different multiplicity or rings with one bond (2 atoms) partial solutions (see below)
# ring solutions - all or nothing; all or one atom; all or one bond (two atoms); partial (needs bond breaking)
# exclusions / pairs

from SMArt.md.incl import Dummy

"""
class Dummy:
    m = None
    p_ch = None
    a_type = None
    coord = np.array([np.nan] * 3)
"""


class Branch:
    def __init__(self, roots, parents, neigh_count, **kwargs):
        self.roots = roots
        self.parents = parents
        self.neigh_count = neigh_count
        self.Nv = kwargs.get('Nv', neigh_count.sum())
        # V_in_branch could be tuple, or frozendict with levels
        self.V_in_branch = kwargs.get('V_in_branch', ())

    def combine_neigh_count(self, other_branch, flag_other_neigh_count = False):
        if flag_other_neigh_count:
            other_neigh_count = other_branch
        else:
            other_neigh_count = other_branch.neigh_count
        n1, n2 = len(self.neigh_count), len(other_neigh_count)
        if n1 > n2:
            return np.vstack([self.neigh_count, np.append(other_neigh_count, [0]*(n1 - n2))])
        else:
            return np.vstack([np.append(self.neigh_count, [0]*(n2 - n1)), other_neigh_count])

    def sum_neigh_count(self, other_branch, flag_other_neigh_count = False):
        return self.combine_neigh_count(other_branch, flag_other_neigh_count).sum(axis = 0)

    def calculate_pair_estimate(self, other_branch):
        stacked_neigh_count = self.combine_neigh_count(other_branch)
        max_match = self.calculate_stacked_neigh_count_pair_estimate(stacked_neigh_count)
        return min(min(self.Nv, other_branch.Nv), max_match)

    @staticmethod
    def calculate_branch_pair_estimate(branches):
        return branches[0].calculate_pair_estimate(branches[1])

    @classmethod
    def calculate_neigh_count_pair_estimate(cls, neigh_count_pair):
        stacked_neigh_count = cls.combine_neigh_count_list(neigh_count_pair)
        return cls.calculate_stacked_neigh_count_pair_estimate(stacked_neigh_count)

    @staticmethod
    def calculate_stacked_neigh_count_pair_estimate(stacked_neigh_count):
        max_match = np.sum(np.min(stacked_neigh_count, axis = 0))
        return max_match

    @staticmethod
    def combine_neigh_count_list(neigh_count_list):
        if neigh_count_list:
            len_lists = [len(ncl) for ncl in neigh_count_list]
            N_max = max(len_lists)
            return np.vstack([np.append(ncl, [0] * (N_max - len_lists[i])) for i, ncl in enumerate(neigh_count_list)])
        return np.array([])

    @staticmethod
    def sum_neigh_count_list(neigh_count_list):
        if neigh_count_list:
            return Branch.combine_neigh_count_list(neigh_count_list).sum(axis = 0)
        return np.array([])

    @staticmethod
    def sum_neigh_count_list_branches(branches):
        temp_neigh_count = np.array([])
        for b in branches:
            temp_neigh_count = b.sum_neigh_count(temp_neigh_count, flag_other_neigh_count = True)
        return temp_neigh_count

    @staticmethod
    def combine_mutually_exclusive_neigh_count_list(neigh_count_list):
        """
        calculates combined neigh_count from a list of neigh_count, such that it takes max at each position
        :param neigh_count_list:
        :return:
        """
        temp_combined_neigh_count_list = Branch.combine_neigh_count_list(neigh_count_list)
        return temp_combined_neigh_count_list.max(axis = 0)

    @staticmethod
    def combine_mutually_exclusive_neigh_count_list_with_max_Nv(neigh_count_list):
        """
        calculates combined neigh_count from a list of neigh_count, such that it takes max at each position
        returns also the max Nv from the neigh count lists
        :param neigh_count_list:
        :return:
        """
        temp_combined_neigh_count_list = Branch.combine_neigh_count_list(neigh_count_list)
        return temp_combined_neigh_count_list.max(axis = 0), temp_combined_neigh_count_list.sum(axis = 1).max()

class TopGraphProperties:
    """
    graph properties: rings, different branches, dummy-match estimates etc...

    attributes within g:
        rings, rings_map
        within each ring:
            breaks, breaks_map, break_parents, independent_parts, independent_parts_map, independent_parts_graph

    attributes within the instance:
        g, available_Vs, non_available_Vs
        first_non_ring_neighbours
        branches, ring_branches - true branches
        combined_FromRing_branches - for each ring / atom in r, a combined branch to all parents that are not in r
        each of ring_branches
            combined_linear_branch - max of neighbouring combined_FromRing_branches within the r
            fake_branches - it's a full branch having a single parent within the r
        break_branches - for each ring
        combined_FromRing_branches - branch that is a combination off branches going out of the r
        combined_branches, linear_combined_branches - for each atom, combination of all available branches
                if v in ring, linear_combined_branches is a bit different as it only takes first neighbours within r
        dummy_estimates
    """
    def __init__(self, g, available_Vs = None, non_available_Vs = None, in_solution_Vs = None,
                 in_solution_dummy_match_Vs = None, remove_V_in_branch=False, **kwargs):
        """
        graph properties: rings, different branches, dummy-match estimates etc...
        :param g:
        :param available_Vs:
        :param non_available_Vs:
        :param in_solution:
        :param kwargs:
            non__av_insol_Vs
            flag_get_ind_ring_parts - if False, no ind_ring_parts (in case partial ring match not allowed)
        """
        if 'non__av_insol_Vs' in kwargs:
            self.available_Vs, self.non_available_Vs, self.in_solution_Vs, self.in_solution_dummy_match_Vs\
                = kwargs['non__av_insol_Vs']
        else:
            self.available_Vs, self.non_available_Vs, self.in_solution_Vs, self.in_solution_dummy_match_Vs = \
                self._g_non__available__insol_Vs(g, available_Vs, non_available_Vs, in_solution_Vs, in_solution_dummy_match_Vs)
        self._in_solution_Vs = self.in_solution_Vs - self.in_solution_dummy_match_Vs
        self.g = g
        self.__get_rings(**kwargs) # gets the rings and rings_map (independent of available_Vs)
        #self.__get_independent_ring_parts(flag_get_ind_ring_branches = False)
        if 'ring_independent_parts_availability_map' in kwargs:
            self.ring_independent_parts_availability_map = kwargs['ring_independent_parts_availability_map']
        else:
            self.__get_independent_ring_parts_availability_map()
        if 1: # this is a generalized way
            self.get_branches()
            self.get_combined_neigh_count_Nv()
            if not self.in_solution_Vs:
                self.get_dummy_match_estimates()
            else:
                self.get_dummy_match_estimates_in_solution()
        else:##################################################### fix
            # this should be if we include all vertices as available
            self.__get_initial_branches()
            self.__get_independent_ring_parts(branches = self.branches, ring_branches = self.ring_branches,
                                              combined_FromRing_branches = self.combined_FromRing_branches)
            self.__initial_combined_neigh_count()
            self.__get_initial_dummy_match_estimate()
        if remove_V_in_branch:
            self.remove_V_in_branch()

    @staticmethod
    def _g_non__available__insol_Vs(g, available_Vs = None, non_available_Vs = None, in_solution_Vs = None,
                                    in_solution_dummy_match_Vs = None):
        if available_Vs:
            available_Vs = frozenset(available_Vs)
            non_available_Vs = frozenset(g.adj) - frozenset(available_Vs)
        elif non_available_Vs:
            available_Vs = frozenset(g.adj) - frozenset(non_available_Vs)
            non_available_Vs = frozenset(non_available_Vs)
        else:
            available_Vs = frozenset(g.adj)
            non_available_Vs = frozenset()
        if in_solution_Vs:
            in_solution_Vs = frozenset(in_solution_Vs)
        else:
            in_solution_Vs = frozenset()
        if in_solution_dummy_match_Vs:
            in_solution_dummy_match_Vs = frozenset(in_solution_dummy_match_Vs)
            in_solution_Vs  = in_solution_Vs | in_solution_dummy_match_Vs
        else:
            in_solution_dummy_match_Vs = frozenset()
        if in_solution_Vs:
            available_Vs = available_Vs - in_solution_Vs
            non_available_Vs = non_available_Vs | in_solution_Vs
        return available_Vs, non_available_Vs, in_solution_Vs, in_solution_dummy_match_Vs

    @classmethod
    def _g_check_branch(cls, g, b, flag_use_V_in_branch = False):
        if flag_use_V_in_branch:
            if not b.V_in_branch:
                do_warn("can't check branch without V_in_branch")
                return
            visited = set(g.adj) - set(b.V_in_branch)
            visited = visited | set(b.roots)
        else:
            visited = set(b.roots)
        temp_neigh_count, _ = cls.get_neigh_count_from_generator(g.BFS(b.parents, visited=visited, flag_v_list = 1))
        assert (temp_neigh_count == b.neigh_count).all()
        if b.V_in_branch:
            assert len(b.V_in_branch) == b.Nv

    def check_branches(self):
        if self.non_available_Vs:
            flag_use_V_in_branch = True
        else:
            flag_use_V_in_branch = False
        for v in self.g.adj:
            for b in self.branches[v].values():
                self._g_check_branch(self.g, b, flag_use_V_in_branch)
            for b in self.ring_branches[v].values():
                self._g_check_branch(self.g, b, flag_use_V_in_branch)
                for fake_b in b.fake_branches.values():
                    self._g_check_branch(self.g, fake_b, flag_use_V_in_branch)
        for r in self.combined_FromRing_branches:
            for b in self.combined_FromRing_branches[r].values():
                self._g_check_branch(self.g, b, flag_use_V_in_branch)
        for r in self.break_branches:
            for ring_ind_part_branches in self.break_branches[r].values():
                for b in ring_ind_part_branches:
                    self._g_check_branch(self.g, b, flag_use_V_in_branch)

    @classmethod
    def _g_get_rings(cls, g, flag_get_ind_ring_parts=True, rerun_get_rings=False):
        if not hasattr(g, 'rings_map') or not hasattr(g, 'rings') or rerun_get_rings:
            g.rings = tuple(g.find_rings_graph(flag_root_at = 1))
            g.rings_map = {}
            for r_graph in g.rings:
                for v in r_graph.adj:
                    if v not in g.rings_map:
                        g.rings_map[v] = []
                    g.rings_map[v].append(r_graph)
            cls._g_get_independent_ring_parts(g, flag_get_ind_ring_branches=False, flag_get_ind_ring_parts=flag_get_ind_ring_parts)

    def __get_rings(self, flag_get_ind_ring_parts=True, rerun_get_rings=False, **kwargs):
        self._g_get_rings(self.g, flag_get_ind_ring_parts=flag_get_ind_ring_parts, rerun_get_rings=rerun_get_rings)

    @classmethod
    def _g_get_independent_ring_parts(cls, g, flag_get_ind_ring_branches = True, branches = None, ring_branches = None,
                                      combined_FromRing_branches = None, flag_get_ind_ring_parts=True):
        break_branches = {}
        for r in g.rings:
            r.breaks = {}
            r.break_parents = {}
            r.independent_parts = []
            if flag_get_ind_ring_parts: # all is one independent part if partial match is not allowed...
                for v in r.adj:
                    if len(r.adj[v]) > 2:
                        for v2 in r.adj[v]:
                            if len(r.adj[v2]) > 2 and (v2, v) not in r.breaks:
                                temp_pair = (v, v2)
                                r.breaks[temp_pair] = []
                                r.break_parents[temp_pair] = []
                                first_neigh = set(r.BFS_l(temp_pair, 1, flag_v_list=True))
                                ind_ring_part_i2del, flag_already_searched_ind_ring_parts = None, None
                                while first_neigh:
                                    temp_first_neigh = first_neigh.pop()
                                    new_ind_ring_part = set(temp_pair)
                                    temp_branch = frozenset(r._BFS(temp_first_neigh, visited=new_ind_ring_part))
                                    ind_ring_part_parents = first_neigh & temp_branch
                                    if len(temp_branch) != len(r.adj) - 2:
                                        new_ind_ring_part = frozenset(new_ind_ring_part)
                                        ind_ring_part_parents.add(temp_first_neigh)
                                        r.break_parents[temp_pair].append(ind_ring_part_parents)
                                        r.breaks[temp_pair].append(tuple(temp_branch))
                                        if ind_ring_part_i2del is None and flag_already_searched_ind_ring_parts is None:
                                            flag_already_searched_ind_ring_parts = True
                                            for temp_ind_ring_part_i, temp_ind_ring_part in enumerate(r.independent_parts):
                                                if v in temp_ind_ring_part and v2 in temp_ind_ring_part:
                                                    ind_ring_part_i2del = temp_ind_ring_part_i
                                                    r.independent_parts.pop(ind_ring_part_i2del)
                                                    break
                                        if ind_ring_part_i2del is None:
                                            r.independent_parts.append(new_ind_ring_part)
                                        else:
                                            r.independent_parts.append(temp_ind_ring_part & new_ind_ring_part)
                                    else:
                                        r.breaks.pop(temp_pair)
                                        r.break_parents.pop(temp_pair)
                                    first_neigh = first_neigh - frozenset(ind_ring_part_parents)
            if not r.independent_parts:
                r.independent_parts = [frozenset(r.adj)]
            r.breaks_map = {}
            for temp_break in r.breaks:
                for v in temp_break:
                    r.breaks_map[v] = temp_break
            r.independent_parts_map = {}
            for v in r.adj:
                r.independent_parts_map[v] = []
            independent_parts_adj = {}
            for ind_ring_part in r.independent_parts:
                independent_parts_adj[ind_ring_part] = []
                for v in ind_ring_part:
                    r.independent_parts_map[v].append(ind_ring_part)
            for ind_ring_part_pair in combinations(r.independent_parts, 2):
                if ind_ring_part_pair[0] & ind_ring_part_pair[1]:
                    for i in range(2):
                        independent_parts_adj[ind_ring_part_pair[i]].append(ind_ring_part_pair[(i + 1) % 2])
            r.independent_parts_graph = Graph(independent_parts_adj)
            if flag_get_ind_ring_branches: # this part depends on available_Vs!!!!!!!!!!!!
                break_branches[r] = {}
                for temp_pair in r.break_parents:
                    for ind_ring_part_parents in r.break_parents[temp_pair]:
                        temp_ind_ring_part_branch = cls.__initial_r_get_ind_ring_branch_fromRB(r, temp_pair,
                    tuple(ind_ring_part_parents), branches, ring_branches, combined_FromRing_branches)
                        break_branches[r][temp_pair].append(temp_ind_ring_part_branch)
        return break_branches # not needed in principle

    def __get_independent_ring_parts(self, **kwargs):
        self.break_branches = self._g_get_independent_ring_parts(self.g, **kwargs)

    @staticmethod
    def _g_get_independent_ring_parts_availability_map(g, non_available_Vs, in_solution_Vs):
        ring_ind_part_non_available_Vs = non_available_Vs - in_solution_Vs
        ring_independent_parts_availability_map = {}
        for r in g.rings:
            flag_all_non_available = False
            roots_ind_part = []
            visited_ind_part = set()
            ring_independent_parts_availability_map[r] = {}
            for ind_part in r.independent_parts:
                if ind_part & ring_ind_part_non_available_Vs:
                    ring_independent_parts_availability_map[r][ind_part] = False
                    ind_part_in_solution_Vs = ind_part & in_solution_Vs
                    if ind_part_in_solution_Vs:
                        if ind_part_in_solution_Vs - set(r.breaks_map):
                            flag_all_non_available = True
                            break
                        roots_ind_part.append(ind_part)
                else:
                    ring_independent_parts_availability_map[r][ind_part] = True
                    if ind_part & in_solution_Vs:
                        visited_ind_part.add(ind_part)
            if flag_all_non_available:
                for ind_part in r.independent_parts:
                    ring_independent_parts_availability_map[r][ind_part] = False
            else:
                for ind_part in r.independent_parts_graph._BFS(roots_ind_part, flag_v_list=True, visited=visited_ind_part):
                    ring_independent_parts_availability_map[r][ind_part] = False
        return ring_independent_parts_availability_map


    def __get_independent_ring_parts_availability_map(self):
        self.ring_independent_parts_availability_map = self._g_get_independent_ring_parts_availability_map(self.g,
                                                                    self.non_available_Vs, self._in_solution_Vs)

    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #

    @classmethod
    def update_neigh_count_VinB(cls, v, l, collect_res, V_in_branch):
        V_in_branch.append(v)
        cls.update_neigh_count(v, l, collect_res)

    @staticmethod
    def update_neigh_count(v, l, collect_res):
        if l==len(collect_res):
            collect_res.append(1)
        else:
            collect_res[l] += 1

    @staticmethod
    def get_neigh_count_from_generator(BFS_generator, **kwargs):
        collect_res = []
        V_in_branch = []
        for v, l in BFS_generator:
            V_in_branch.append(v)
            if l==len(collect_res):
                collect_res.append(1)
            else:
                collect_res[l] += 1
        return np.array(collect_res, dtype=int), V_in_branch

    @staticmethod
    def _get_neigh_count_from_generator_VinB_list(BFS_generator):
        collect_res = []
        V_in_branch = []
        for v, l in BFS_generator:
            V_in_branch.append(v)
            if l==len(collect_res):
                collect_res.append(1)
            else:
                collect_res[l] += 1
        return np.array(collect_res, dtype=int), V_in_branch

    @staticmethod
    def _get_neigh_count_from_generator_VinB_frozendict(BFS_generator):
        collect_res = []
        V_in_branch = []
        for v, l in BFS_generator:
            V_in_branch.append((v,l))
            if l==len(collect_res):
                collect_res.append(1)
            else:
                collect_res[l] += 1
        return np.array(collect_res, dtype=int), frozendict(V_in_branch)

    @classmethod
    def _g_get_neigh_count_Vpair(cls, g, v0, v1):
        visited = set([v0])
        return cls.get_neigh_count_from_generator(g.BFS(v1, visited=visited))

    def get_neigh_count_Vpair(self, v0, v1):
        visited = set([v0])
        return self.get_neigh_count_from_generator(self.g.BFS(v1, visited=visited))

    @classmethod
    def _g_get_neigh_count_roots_parents(cls, g, roots, parents):
        visited = set(roots)
        temp_generator = g.BFS(parents, visited=visited, flag_v_list = True)
        return cls.get_neigh_count_from_generator(temp_generator)

    def get_neigh_count_roots_parents(self, roots, parents):
        visited = set(roots)
        temp_generator = self.g.BFS(parents, visited=visited, flag_v_list = True)
        return self.get_neigh_count_from_generator(temp_generator)

    @classmethod
    def __g_non__available_to_visited_v0(cls, g, v0, available_Vs = None, non_available_Vs = None):
        return cls.__g_non__available_to_visited_roots(g, [v0], available_Vs, non_available_Vs)

    @classmethod
    def __g_non__available_to_visited_roots(cls, g, roots, available_Vs = None, non_available_Vs = None):
        non_available_Vs = cls.__g_non__available_Vs(g, available_Vs, non_available_Vs)
        return cls.__roots2visited(roots, non_available_Vs)

    @staticmethod
    def __g_non__available_Vs(g, available_Vs = None, non_available_Vs = None):
        if non_available_Vs:
            return frozenset(non_available_Vs)
        elif available_Vs:
            non_available_Vs = frozenset(g.adj) - frozenset(available_Vs)
            return non_available_Vs
        else:
            return frozenset()

    @staticmethod
    def __roots2visited(roots, non_available_Vs):
        if non_available_Vs:
            visited = set(roots) | non_available_Vs
        else:
            visited = set(roots)
        return visited

    def __available_to_visited_v0(self, v0):
        return self.__roots2visited([v0], self.non_available_Vs)

    def __available_to_visited_roots(self, roots, available_Vs):
        return self.__roots2visited(roots, self.non_available_Vs)

    @classmethod
    def _g_get_branch_roots_parents(cls, g, roots, parents, available_Vs = None, non_available_Vs = None):
        visited = cls.__g_non__available_to_visited_roots(g, roots, available_Vs, non_available_Vs)
        temp_generator = g.BFS(parents, visited=visited, flag_v_list = True)
        neigh_count, V_in_branch = cls.get_neigh_count_from_generator(temp_generator)
        return Branch(roots, parents, neigh_count, V_in_branch = tuple(V_in_branch))

    @classmethod
    def _r_get_branch_roots_parents_fromBR(cls, r, roots, parents, branches, ring_branches, combined_FromRing_branches):
        BFS_generator = r.BFS(parents, visited = set(roots), flag_v_list = True)
        V_in_branch = []
        neigh_count_list = []
        for v, l in BFS_generator:
            temp_combined_FromRing_b = combined_FromRing_branches[r][v]
            V_in_branch.extend(temp_combined_FromRing_b.V_in_branch)
            temp_neigh_count = np.append([0] * l, temp_combined_FromRing_b.neigh_count)
            neigh_count_list.append(temp_neigh_count)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        return Branch(roots, parents, b_neigh_count, V_in_branch = V_in_branch)

    @staticmethod
    def __initial_r_get_sorted_ring_linear_combined_branches(r, ring_branches):
        temp_list = [ring_branches[v][r].combined_linear_branch for v in r.adj]
        temp_list.sort(key = lambda x:x.tot_Nv, reverse=True)
        return temp_list

# most general cases - not using already known branches
    def __get_branch_Vpair(self, v0, v1): # does not allow for available_Vs argument - all Vs available
        roots = (v0,)
        parents = (v1,)
        neigh_count, V_in_branch = self.get_neigh_count_Vpair(v0, v1)
        return Branch(roots, parents, neigh_count, V_in_branch = V_in_branch)

    def get_branch_Vpair(self, v0, v1):
        roots = (v0,)
        parents = (v1,)
        visited = self.__roots2visited([v0], self.non_available_Vs)
        temp_generator = self.g.BFS(v1, visited=visited)
        neigh_count, V_in_branch = self.get_neigh_count_from_generator(temp_generator)
        return Branch(roots, parents, neigh_count, V_in_branch = V_in_branch)

    def __get_branch_roots_parents(self, roots, parents): # does not allow for available_Vs argument - all Vs available
        visited = set(roots)
        temp_generator = self.g.BFS(parents, visited=visited, flag_v_list = True)
        neigh_count, V_in_branch = self.get_neigh_count_from_generator(temp_generator)
        return Branch(roots, parents, neigh_count, V_in_branch = V_in_branch)

    def get_branch_roots_parents(self, roots, parents):
        visited = self.__roots2visited(roots, self.non_available_Vs)
        temp_generator = self.g.BFS(parents, visited=visited, flag_v_list = True)
        neigh_count, V_in_branch = self.get_neigh_count_from_generator(temp_generator)
        return Branch(roots, parents, neigh_count, V_in_branch = V_in_branch)

    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions #

    ####################################################################################################################
    # __initaial all available atoms # __initaial all available atoms # __initaial all available atoms #################
    # __initaial all available atoms # __initaial all available atoms # __initaial all available atoms #################
    ####################################################################################################################

    # getting branches - assuming only bonded pairs and non-fake branches
    # get branch recursively using known branches
    def __get_initial_first_non_ring_neighbour(self):
        self.first_non_ring_neighbours = {}
        for v0 in self.g.adj:
            first_neighbours = self.g.adj[v0]
            if v0 in self.g.rings_map:
                first_neighbours = set(self.g.adj[v0])
                for r in self.g.rings_map[v0]:
                    first_neighbours = first_neighbours - set(r.adj[v0])
            self.first_non_ring_neighbours[v0] = tuple(first_neighbours)

    def __get_initial_branches(self): # does not allow for available_Vs argument - all Vs available
        self.__get_initial_first_non_ring_neighbour()
        for v0 in self.g.adj:
            if v0 in self.g.rings_map:
                for r in self.g.rings_map[v0]:
                    self.__get_ring_branches_recursive(v0, r)
            for v1 in self.first_non_ring_neighbours[v0]:
                self.__get_true_branch_Vpair_recursive(v0, v1)
        self.__initial_get_ring_combined_linear_branches_fromBR()

    def __get_true_branch_Vpair_recursive(self, v0, v1):
        if v1 in self.branches[v0]:
            return self.branches[v0][v1]
        V_in_branch = [v1]
        neigh_count_list = [[1]]
        if v1 in self.g.rings_map:
            for r in self.g.rings_map[v1]:
                assert v0 not in r.adj
                temp_branch = self.__get_ring_branches_recursive(v1, r)
                V_in_branch.extend(temp_branch.V_in_branch)
                temp_neigh_count = np.append([0], temp_branch.neigh_count)
                neigh_count_list.append(temp_neigh_count)
        for v_next in self.first_non_ring_neighbours[v1]:
            if v_next!=v0:
                temp_branch = self.__get_true_branch_Vpair_recursive(v1, v_next)
                V_in_branch.extend(temp_branch.V_in_branch)
                temp_neigh_count = np.append([0], temp_branch.neigh_count)
                neigh_count_list.append(temp_neigh_count)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        linear_branch = Branch((v0,), (v1, ), b_neigh_count, V_in_branch = V_in_branch)
        self.branches[v0][v1] = linear_branch
        return linear_branch

    def __get_ring_branches_recursive(self, v0, r):
        if r in self.ring_branches[v0]:
            return self.ring_branches[v0][r]
        V_in_branch = []
        neigh_count_list = []
        for v, l in r.BFS(r.adj[v0], visited = set([v0]), flag_v_list = True):
            temp_comb_branch = self.__initial_get_ring_combined_FromRing_branch_recursive(r, v)
            V_in_branch.extend(temp_comb_branch.V_in_branch)
            temp_neigh_count = np.append([0] * l, temp_comb_branch.neigh_count)
            neigh_count_list.append(temp_neigh_count)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        ring_branch = Branch((v0,), tuple(r.adj[v0]), b_neigh_count, V_in_branch = V_in_branch)
        self.ring_branches[v0][r] = ring_branch
        # do fake linear branches from each first neighbour within the ring
        ring_branch.fake_branches = {}
        for v1 in r.adj[v0]:
            neigh_count_list = []
            for v, l  in r.BFS(v1, visited = set([v0])):
                temp_comb_branch = self.__initial_get_ring_combined_FromRing_branch_recursive(r, v)
                temp_neigh_count = np.append([0] * l, temp_comb_branch.neigh_count)
                neigh_count_list.append(temp_neigh_count)
            b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
            within_ring_branch = Branch((v0,), (v1, ), b_neigh_count, V_in_branch = V_in_branch)
            ring_branch.fake_branches[v1] = within_ring_branch
        return ring_branch

    def __initial_get_FromRing_branches_recursive(self, r, v):
        for ring_outside_r in self.g.rings_map[v]:
            if ring_outside_r != r:
                yield self.__get_ring_branches_recursive(v, ring_outside_r)
        for v_outside_r in self.first_non_ring_neighbours[v]:
            yield self.__get_true_branch_Vpair_recursive(v, v_outside_r)

    def __initial_get_ring_combined_FromRing_branch_recursive(self, r, v):
        if v in self.combined_FromRing_branches[r]:
            return self.combined_FromRing_branches[r][v]
        neigh_count_list = []
        V_in_branch = [v]
        for temp_branch in self.__initial_get_FromRing_branches_recursive(r, v):
            neigh_count_list.append(temp_branch.neigh_count)
            V_in_branch.extend(temp_branch.V_in_branch)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        b_neigh_count = np.append([1], b_neigh_count)
        V_in_branch = tuple(V_in_branch)
        self.combined_FromRing_branches[r][v] = Branch(tuple(r.adj[v]), (v,), b_neigh_count, V_in_branch = V_in_branch)
        return self.combined_FromRing_branches[r][v]

    @classmethod
    def __initial_r_get_ind_ring_branch_fromRB(cls, r, temp_pair, ind_ring_part_parents, branches, ring_branches,
                                      combined_FromRing_branches):
        temp_b = cls._r_get_branch_roots_parents_fromBR(r, temp_pair, ind_ring_part_parents, branches, ring_branches, combined_FromRing_branches)
        temp_b.tot_Nv = temp_b.Nv
        for v in temp_pair:
            FromRing_b = combined_FromRing_branches[r][v]
            temp_b.tot_Nv += FromRing_b.Nv
        return temp_b

    def __initial_get_ring_combined_linear_branches_fromBR(self):
        self.__initial_g_get_ring_combined_linear_branches_fromBR(self.g, self.ring_branches, self.combined_FromRing_branches)

    @classmethod
    def __initial_g_get_ring_combined_linear_branches_fromBR(cls, g, ring_branches, combined_FromRing_branches):
        for r in g.rings:
            for v in r.adj:
                neigh_count_list = [combined_FromRing_branches[r][v1].neigh_count for v1 in r.adj[v]]
                Nv_max = 0
                for v1 in r.adj[v]:
                    if Nv_max < combined_FromRing_branches[r][v1].Nv:
                        Nv_max = combined_FromRing_branches[r][v1].Nv
                combined_neigh_count = Branch.combine_mutually_exclusive_neigh_count_list(neigh_count_list)
                combined_linear_branch = Branch((v,), ring_branches[v][r].parents, combined_neigh_count, Nv = Nv_max)
                combined_linear_branch.tot_Nv = combined_linear_branch.Nv + combined_FromRing_branches[r][v].Nv# - 1 (this minus one is actually not wrong, but we use this to get to the total number of atoms, including the roots...)
                ring_branches[v][r].combined_linear_branch = combined_linear_branch

    def __initial_combined_neigh_count(self):
        #self.ring_combined_branches, self.linear_combined_branches =
        self.combined_branches = self.__initial_g_get_combined_neigh_count(self.g, self.branches, self.ring_branches)

    @staticmethod
    def __initial_g_get_combined_neigh_count(g, branches, ring_branches):
        combined_branches = {}
#        ring_combined_branches = {}
#        linear_combined_branches = {}
        for v in g.adj:
            Nv_sum = 0
            neigh_count_list = []
            for b in branches[v].values():
                neigh_count_list.append(b.neigh_count)
                Nv_sum+=b.Nv
            temp_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
            if ring_branches[v]:
                neigh_count_list = [temp_neigh_count]
                Nv_sum_normal = Nv_sum
                for b in ring_branches[v].values():
                    neigh_count_list.append(b.neigh_count)
                    Nv_sum_normal+=b.Nv
                comb_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
                ring_combined_branch = Branch((v,), tuple(g.adj[v]), comb_neigh_count, Nv = Nv_sum_normal)
                combined_branches[v] = ring_combined_branch
                """
                ring_combined_branches[v] = ring_combined_branch
                neigh_count_list = [temp_neigh_count]
                Nv_sum_linear = Nv_sum
                for b in ring_branches[v].values():
                    neigh_count_list.append(b.combined_linear_branch.neigh_count)
                    Nv_sum_linear+=b.combined_linear_branch.Nv
                comb_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
                linear_combined_branch = Branch((v,), tuple(g.adj[v]), comb_neigh_count, Nv = Nv_sum_linear)
                linear_combined_branches[v] = linear_combined_branch
                """
            else:
                linear_combined_branch = Branch((v,), tuple(g.adj[v]), temp_neigh_count, Nv = Nv_sum)
                combined_branches[v] = linear_combined_branch
                """
                linear_combined_branches[v] = linear_combined_branch
                ring_combined_branches[v] = None
                """
        return combined_branches
#        return ring_combined_branches, linear_combined_branches

    def __get_initial_dummy_match_estimate(self):
        self.dummy_estimates = {}
        sorted_RLCBs = {}
        for r in self.g.rings:
            sorted_RLCBs[r] = self.__initial_r_get_sorted_ring_linear_combined_branches(r, self.ring_branches)
        for v in self.g.adj:
            max_Nv_neigh_branches = 0
            for r in self.ring_branches[v]:
                Vs_independent = []
                for temp_ind_ring_part in r.independent_parts:
                    if v in temp_ind_ring_part: # finding independent ring parts with this atom
                        for v1 in temp_ind_ring_part: # looping over each atom of the independent part
                            if v1 in r.breaks_map:
                                if v not in r.breaks_map[v1]:
                                    for break_branch in self.break_branches[r][r.breaks_map[v1]]:
                                        if break_branch.parents[0] not in temp_ind_ring_part:
                                            if max_Nv_neigh_branches < break_branch.tot_Nv:
                                                max_Nv_neigh_branches = break_branch.tot_Nv
                            #else:Vs_independent.append(v1)
                            if v1 != v:
                                Vs_independent.append(v1)
                temp_max_Nv = self._r_get_max_Nv_ring_linear_combined_branch(r, v, Vs_independent, sorted_RLCBs[r],
                                                                             self.combined_FromRing_branches)
                if max_Nv_neigh_branches < temp_max_Nv:
                    max_Nv_neigh_branches = temp_max_Nv
            for b in self.branches[v].values():
                if b.Nv > max_Nv_neigh_branches:
                    max_Nv_neigh_branches = b.Nv
            self.dummy_estimates[v] = max_Nv_neigh_branches

    ####################################################################################################################
    #                                                   END                                                            #
    # __initaial all available atoms # __initaial all available atoms # __initaial all available atoms #################
    ####################################################################################################################

    ####################################################################################################################
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #################
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #################
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #################
    ####################################################################################################################

    ####################################################################################################################
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    ####################################################################################################################

    def get_first_non_ring_neighbour(self):
        self.first_non_ring_neighbours = {}
        for v0 in self.g.adj:
            first_neighbours = set(self.g.adj[v0]) - self.non_available_Vs
            if v0 in self.g.rings_map:
                for r in self.g.rings_map[v0]:
                    first_neighbours = first_neighbours - set(r.adj[v0])
            self.first_non_ring_neighbours[v0] = tuple(first_neighbours)

    def combine_neigh_count_Nv_general(self, neigh_count_list, Nv_sum, temp_branch):
        neigh_count_list.append(temp_branch.neigh_count)
        Nv_sum += temp_branch.Nv
        return Nv_sum

    def combine_neigh_count_Nv_VinB_general(self, V_in_branch, neigh_count_list, Nv_sum, temp_branch):
        V_in_branch.extend(temp_branch.V_in_branch)
        return self.combine_neigh_count_Nv_general(neigh_count_list, Nv_sum, temp_branch)

    def update_neigh_count_Nv_general(self, neigh_count_list, Nv_sum, temp_branch, l = 1):
        temp_neigh_count = np.append([0] * l, temp_branch.neigh_count)
        neigh_count_list.append(temp_neigh_count)
        Nv_sum += temp_branch.Nv
        return Nv_sum

    def update_neigh_count_Nv_VinB_general(self, V_in_branch, neigh_count_list, Nv_sum, temp_branch, l = 1):
        V_in_branch.extend(temp_branch.V_in_branch)
        return self.update_neigh_count_Nv_general(neigh_count_list, Nv_sum, temp_branch, l=l)

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################
    def remove_V_in_branch(self):
        V_in_branch_empty = frozenset()
        for branches in self.branches.values():
            for temp_branch in branches.values():
                temp_branch.V_in_branch = V_in_branch_empty
        for branches in self.ring_branches.values():
            for temp_branch in branches.values():
                temp_branch.V_in_branch = V_in_branch_empty
                temp_branch.combined_linear_branch.V_in_branch = V_in_branch_empty
                for temp_fake_branch in temp_branch.fake_branches.values():
                    temp_fake_branch.V_in_branch = V_in_branch_empty

    def get_branches(self):
        if self.in_solution_Vs:
            temp_gen = self.g.BFS_l(self.in_solution_Vs, 1, flag_v_list=True, visited=set(self.non_available_Vs))
            self.Vs_of_interest =  tuple(temp_gen)
            temp_gen = self.g.BFS_l(self.in_solution_Vs, 1, flag_v_list=True, visited=set(self.non_available_Vs))
            self.Vs_of_int_sol =  list(temp_gen)
            self.Vs_of_int_sol.extend(self.in_solution_Vs)
            self.Vs_of_int_sol = tuple(self.Vs_of_int_sol)
            #self.Vs_of_interest = self.Vs_of_int_sol
        else:
            self.Vs_of_interest =  tuple(self.available_Vs)
            self.Vs_of_int_sol = self.Vs_of_interest
        #self.Vs_of_interest =  tuple(self.available_Vs)
        self.get_first_non_ring_neighbour()
        self.branches = {}
        self.ring_branches = {}
        self.combined_FromRing_branches = {}
        self.dummy_branches = {}
        #for v0 in self.Vs_of_interest:
        for v0 in self.Vs_of_int_sol:
            self.branches[v0] = {}
            self.ring_branches[v0] = {}
#        for v0 in self.Vs_of_interest:
        for v0 in self.Vs_of_int_sol:
            if v0 in self.in_solution_dummy_match_Vs:
                self.get_dummy_branch(v0)
                continue
            if v0 in self.g.rings_map:
                for r in self.g.rings_map[v0]:
                    self.get_ring_branch_recursive(v0, r)
            for v1 in self.first_non_ring_neighbours[v0]:
                self.get_nonring_branch_recursive(v0, v1)
        self.get_ring_combined_linear_branch_tot_Nv()
        self.get_ring_break_branches()

    def get_dummy_branch(self, v0):
        if v0 in self.g.rings_map:
            protected_Vs = set()
            if v0 in self.g.rings_map:
                for r in self.g.rings_map[v0]:
                    temp_set = set(r.adj) & self._in_solution_Vs
                    if len(temp_set) == 1:
                        other_atom_in_r = temp_set.pop()
                        protected_Vs = protected_Vs | frozenset(r.adj[other_atom_in_r])
                        for ind_ring_part in r.independent_parts_map[other_atom_in_r]:
                            if v0 not in ind_ring_part:
                                protected_Vs = protected_Vs | ind_ring_part
            visited = protected_Vs | self.non_available_Vs
        else:
            visited = set(self.non_available_Vs)
        self.dummy_branches[v0] = tuple(self.g._BFS(v0, visited=visited))[1:]

    def get_nonring_branch_recursive(self, v0, v1):
        if v0 not in self.branches:
            self.branches[v0] = {}
        if v1 in self.branches[v0]:
            return self.branches[v0][v1]
        V_in_branch = [v1]
        neigh_count_list = [[1]]
        Nv_sum = 1
        if v1 in self.g.rings_map:
            for r in self.g.rings_map[v1]:
                assert v0 not in r.adj
                temp_branch = self.get_ring_branch_recursive(v1, r)
                Nv_sum = self.update_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, temp_branch)#####################
        for v_next in self.first_non_ring_neighbours[v1]:
            if v_next!=v0 and v_next in self.available_Vs:
                temp_branch = self.get_nonring_branch_recursive(v1, v_next)
                Nv_sum = self.update_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, temp_branch)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        linear_branch = Branch((v0,), (v1, ), b_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)
        self.branches[v0][v1] = linear_branch
        return linear_branch

    def get_ring_branch_recursive(self, v0, r):
        if v0 not in self.ring_branches:
            self.ring_branches[v0] = {}
        if r in self.ring_branches[v0]:
            return self.ring_branches[v0][r]
        roots = (v0,)
        visited, parents = self._r_get_visited_from_roots_parents(roots, r.adj[v0], r,
                                                self.ring_independent_parts_availability_map[r], self.non_available_Vs)
        temp_BFS = r.BFS(parents, visited = set(visited), flag_v_list = True)
        # do ring_ring_branch
        ring_ring_branch = self.get_ring_ring_branch_from_generator_VinB_recursive(roots, parents, r, temp_BFS)
        # do combined linear branch
        ring_combined_linear_branch = self.get_ring_combined_linear_branch_recursive(r, v0)
        # combine ring_ring and combined_linear into ring_branch
        ring_branch = self._get_ring_branch_combine_ring_ring_linear(ring_ring_branch, ring_combined_linear_branch)
        # do fake linear branches from each first neighbour within the ring
        ring_branch.fake_branches = {}
        for v1 in r.adj[v0]:
            if v1 in ring_branch.ring_branch.parents:
                temp_BFS = r.BFS([v1], visited = set(visited), flag_v_list = True)
                fake_ring_branch = self.get_ring_ring_branch_from_generator_recursive(roots, (v1,), r,
                                                                        ring_branch.ring_branch.V_in_branch, temp_BFS)
                ring_branch.fake_branches[v1] = fake_ring_branch
        self.ring_branches[v0][r] = ring_branch
        return ring_branch

    def get_ring_ring_branch_from_generator_VinB_recursive(self, roots, parents, r, BFS_generator):
        V_in_branch = []
        neigh_count_list = []
        Nv_sum = 0
        for v, l in BFS_generator:
            temp_comb_branch = self.get_ring_combined_FromRing_branch_recursive(r, v)
            Nv_sum = self.update_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, temp_comb_branch, l)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        return Branch(roots, parents, b_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)

    # same as above, just does not take care about V_in_branches - this is just copied
    def get_ring_ring_branch_from_generator_recursive(self, roots, parents, r, V_in_branch, BFS_generator):
        neigh_count_list = []
        Nv_sum = 0
        for v, l in BFS_generator:
            temp_comb_branch = self.get_ring_combined_FromRing_branch_recursive(r, v)
            Nv_sum = self.update_neigh_count_Nv_general(neigh_count_list, Nv_sum, temp_comb_branch, l)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        return Branch(roots, parents, b_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)

    # this would be wrong as ring_combined_linear_branch assumes v0 -> v1 + FromRing branches
    # so only atoms that have no neighbours in solution could be available for such branches
    # as it should be v0 and v1 bonded pair within the ring - which already gives 2 atoms
    """
    def __get_flag_ring_combined_linear_branch(self, r, v):
        temp__in_sol__1st_neigh = self.in_solution_Vs & set(r.adj[v])
        if temp__in_sol__1st_neigh:
            if len(temp__in_sol__1st_neigh) == 1:
                for temp_v in temp__in_sol__1st_neigh:
                    if not self.in_solution_Vs & set(r.adj[temp_v]):
                        return True
        else:
            return True
        return False
    """

    def get_ring_combined_linear_branch_recursive(self, r, v):
        neigh_count_list = []
        V_in_branch = []
        Nv_max = 0
        parents = []
        if not self._in_solution_Vs & set(r.adj[v]): # ring_combined_linear_branch assumes v0 -> v1 + FromRing branches
            for v1 in r.adj[v]:
                if v1 in self.available_Vs:
                    parents.append(v1)
                    temp_FromRing_branch = self.get_ring_combined_FromRing_branch_recursive(r, v1)
                    neigh_count_list.append(temp_FromRing_branch.neigh_count)
                    V_in_branch.extend(temp_FromRing_branch.V_in_branch)
                    if Nv_max < temp_FromRing_branch.Nv:
                        Nv_max = temp_FromRing_branch.Nv
        parents = tuple(parents)
        if parents:
            combined_neigh_count = Branch.combine_mutually_exclusive_neigh_count_list(neigh_count_list)
        else:
            combined_neigh_count = np.array([])
        return Branch((v,), parents, combined_neigh_count, V_in_branch = V_in_branch, Nv = Nv_max)

    def get_ring_combined_FromRing_branch_recursive(self, r, v):
        if r not in self.combined_FromRing_branches:
            self.combined_FromRing_branches[r] = {}
        if v in self.combined_FromRing_branches[r]:
            return self.combined_FromRing_branches[r][v]
        neigh_count_list = []
        V_in_branch = [v]
        Nv_sum = 1
        for temp_branch in self.get_FromRing_branches_recursive(r, v):
            Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, temp_branch)
        b_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        b_neigh_count = np.append([1], b_neigh_count)
        V_in_branch = tuple(V_in_branch)
        combined_FromRing_b = Branch(tuple(r.adj[v]), (v,), b_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)
        self.combined_FromRing_branches[r][v] = combined_FromRing_b
        return self.combined_FromRing_branches[r][v]

    @staticmethod
    def _get_ring_branch_combine_ring_ring_linear(ring_ring_branch, ring_combined_linear_branch):
        combined_neigh_count = Branch.combine_mutually_exclusive_neigh_count_list((ring_ring_branch.neigh_count,
                                                                            ring_combined_linear_branch.neigh_count))
        parents = frozenset(ring_ring_branch.parents) | frozenset(ring_combined_linear_branch.parents)
        V_in_branch = frozenset(ring_ring_branch.V_in_branch) | frozenset(ring_combined_linear_branch.V_in_branch)
        Nv_max = max(ring_ring_branch.Nv, ring_combined_linear_branch.Nv)
        ring_branch = Branch(ring_ring_branch.roots, parents, combined_neigh_count, V_in_branch = V_in_branch, Nv = Nv_max)
        ring_branch.ring_branch, ring_branch.combined_linear_branch = ring_ring_branch, ring_combined_linear_branch
        return ring_branch

    def get_ring_combined_linear_branch_tot_Nv(self):
        for v in self.ring_branches:
            #if v not in self.Vs_of_interest:continue
            if v not in self.Vs_of_int_sol:continue
            for r in self.ring_branches[v]:
                ring_branch = self.ring_branches[v][r]
                combined_linear_branch = ring_branch.combined_linear_branch
                temp_FromRing_branch = self.get_ring_combined_FromRing_branch_recursive(r, v)
                combined_linear_branch.tot_Nv = combined_linear_branch.Nv + temp_FromRing_branch.Nv# - 1 (this minus one is actually not wrong, but we use this to get to the total number of atoms, including the roots...)

    # ring break branches
    def get_ring_break_branches(self):
        self.break_branches = {}
        for r in self.g.rings:
            self.break_branches[r] = {}
            for temp_break in r.break_parents:
                temp_break_branches = []
                for parents in r.break_parents[temp_break]:
                    temp_gen, parents = self._r_get_BFS_generator_roots_parents_non__available_Vs(temp_break, parents,
                                            r, self.ring_independent_parts_availability_map[r], self.non_available_Vs)
                    temp_branch = self.get_ring_ring_branch_from_generator_VinB_recursive(temp_break, parents, r, temp_gen)
                    temp_branch.tot_Nv = temp_branch.Nv
                    for v in temp_branch.roots:
                        FromRing_b = self.get_ring_combined_FromRing_branch_recursive(r, v)
                        temp_branch.tot_Nv += FromRing_b.Nv
                    temp_break_branches.append(temp_branch)
                self.break_branches[r][temp_break] = temp_break_branches

    ####################################################################################################################
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    ####################################################################################################################
    @staticmethod
    def _r_get_visited_from_roots_parents(roots, parents, r, independent_parts_availability_map, non_available_Vs):
        initial_ind_parts = []
        for v in roots:
            for ind_part in r.independent_parts_map[v]:
                if independent_parts_availability_map[ind_part] and ind_part not in initial_ind_parts:
                    initial_ind_parts.append(ind_part)
        visited_ind_parts = set(ind_part for ind_part in independent_parts_availability_map
                             if not independent_parts_availability_map[ind_part]) # all ind parts that are not available
        available_Vs = set()
        for ind_part in r.independent_parts_graph._BFS(initial_ind_parts, visited = visited_ind_parts, flag_v_list = 1):
            available_Vs = available_Vs | ind_part
        visited = set(r.adj) - available_Vs
        visited = visited | set(roots) | non_available_Vs
        return visited, tuple(set(parents) - visited)

    @classmethod
    def _r_get_BFS_generator_roots_parents_non__available_Vs(cls, roots, parents, r, independent_parts_availability_map,
                                                             non_available_Vs):
        visited, parents = cls._r_get_visited_from_roots_parents(roots, parents, r, independent_parts_availability_map,
                                                                 non_available_Vs)
        return r.BFS(parents, visited = visited, flag_v_list = True), parents

    def get_FromRing_branches_recursive(self, r, v):
        for ring_outside_r in self.g.rings_map[v]:
            if ring_outside_r != r:
                yield self.get_ring_branch_recursive(v, ring_outside_r)
        for v_outside_r in self.first_non_ring_neighbours[v]:
            yield self.get_nonring_branch_recursive(v, v_outside_r)

    # same as above, but split...
    def get_FromRing_ring_branches_recursive(self, r, v):
        for ring_outside_r in self.g.rings_map[v]:
            if ring_outside_r != r:
                yield self.get_ring_branch_recursive(v, ring_outside_r)

    def get_FromRing_nonring_branches_recursive(self, v):
        for v_outside_r in self.first_non_ring_neighbours[v]:
            yield self.get_nonring_branch_recursive(v, v_outside_r)

    ####################################################################################################################
    ####################################################################################################################

    # combined branches
    def _get_combined_neigh_count_Nv_allowed_parents(self, v, allowed_parents):
        Nv_sum = 0
        neigh_count_list = []
        V_in_branch = []
        parents = []
        for b in self.branches[v].values():
            if b.parents[0] in allowed_parents:
                Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, b)
                parents.append(b.parents[0])
        if v in self.g.rings_map:
            Nv_sum_lin = Nv_sum
            neigh_count_list_lin = list(neigh_count_list)
            V_in_branch_lin = list(V_in_branch)
            parents_lin = list(parents)
            for r, b in self.ring_branches[v].items():
                flag_parents_in_allowed = True############## check DEBUGGING
                flag_1_parent_in_allowed = False
                for temp_parent in b.parents:
                    if temp_parent not in allowed_parents:
                        flag_parents_in_allowed = False
                        break
                if flag_parents_in_allowed:
                    Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, b)
                    parents.extend(b.parents)
                    Nv_sum_lin = self.combine_neigh_count_Nv_VinB_general(V_in_branch_lin, neigh_count_list_lin,
                                                                          Nv_sum_lin, b.combined_linear_branch)
                    parents_lin.extend(b.combined_linear_branch.parents)
                else:
                    for temp_parent in b.parents:
                        if temp_parent in allowed_parents:
                            fake_b = b.fake_branches[temp_parent]
                            Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, fake_b)
                            parents.extend(fake_b.parents)
                            Nv_sum_lin = self.combine_neigh_count_Nv_VinB_general(V_in_branch_lin, neigh_count_list_lin,
                                                                                  Nv_sum_lin, fake_b)
                            parents_lin.extend(b.combined_linear_branch.parents)
            temp_neigh_count_lin = Branch.sum_neigh_count_list(neigh_count_list_lin)
            ring_linear_combined_branch = Branch((v,), tuple(parents_lin), temp_neigh_count_lin,
                                                           V_in_branch = V_in_branch_lin, Nv = Nv_sum_lin)
        else:
            ring_linear_combined_branch = None
        temp_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
        combined_branch = Branch((v,), tuple(parents), temp_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)
        combined_branch.ring_linear_combined_branch = ring_linear_combined_branch
        return combined_branch

    def get_combined_neigh_count_Nv(self):
        self.combined_branches = {}
        #for v in self.Vs_of_interest:
        for v in self.Vs_of_int_sol: # combined list of Vs of int and in sol Vs
            Nv_sum = 0
            neigh_count_list = []
            V_in_branch = []
            parents = []
            for b in self.branches[v].values():
                Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, b)
                parents.extend(b.parents)
            # splitting the combined branch in ring and linear
            if v in self.g.rings_map:
                Nv_sum_lin = Nv_sum
                neigh_count_list_lin = list(neigh_count_list)
                V_in_branch_lin = list(V_in_branch)
                parents_lin = list(parents)
                for r, b in self.ring_branches[v].items():
                    Nv_sum = self.combine_neigh_count_Nv_VinB_general(V_in_branch, neigh_count_list, Nv_sum, b)
                    parents.extend(b.parents)
                    Nv_sum_lin = self.combine_neigh_count_Nv_VinB_general(V_in_branch_lin, neigh_count_list_lin,
                                                                          Nv_sum_lin, b.combined_linear_branch)
                    parents_lin.extend(b.combined_linear_branch.parents)
                temp_neigh_count_lin = Branch.sum_neigh_count_list(neigh_count_list_lin)
                ring_linear_combined_branch = Branch((v,), tuple(parents_lin), temp_neigh_count_lin,
                                                               V_in_branch = V_in_branch_lin, Nv = Nv_sum_lin)
            else:
                ring_linear_combined_branch = None
            temp_neigh_count = Branch.sum_neigh_count_list(neigh_count_list)
            self.combined_branches[v] = Branch((v,), tuple(parents), temp_neigh_count, V_in_branch = V_in_branch, Nv = Nv_sum)
            self.combined_branches[v].ring_linear_combined_branch = ring_linear_combined_branch

    ################ separate groups
    def get_separate_groups_brute_force(self):
        groups = []
        for v in self.in_solution_Vs:
            if not self.combined_branches[v].V_in_branch:
                continue
            grp = frozenset(self.combined_branches[v].V_in_branch)
            for i in range(len(groups)):
                temp_grp = groups[i]
                if grp == temp_grp:
                    grp = None
                    break
                common = grp & temp_grp
                if common:
                    if common != temp_grp:
                        groups[i] = temp_grp - grp
                    if common not in groups:
                        groups.append(common)
                    grp = grp - temp_grp
                    if not grp:
                        break
            if grp:
                groups.append(grp)
        self.check_separate_groups(groups)
        return groups

    def __get_tgp_no_dummies(self):
        return TopGraphProperties(self.g, available_Vs=self.available_Vs, in_solution_Vs=self.in_solution_Vs)

    def get_separate_groups(self, **kwargs):
        tgp_no_dummies = kwargs.get('tgp_no_dummies')
        if not tgp_no_dummies:
            if not self.in_solution_dummy_match_Vs:
                tgp_no_dummies = self
            else:
                tgp_no_dummies = self.__get_tgp_no_dummies()
        #in_solution_Vs_groups = {}
        parents__common_Vs_of_interest_map = {}
        #groups = {}
        Vs_of_interest__in_solution_Vs_map = {}
        Vs_of_interest = frozenset(self.Vs_of_interest)
        in_solution_Vs = frozenset(self._in_solution_Vs)
        for temp_parent in Vs_of_interest:
            common_Vs_in_sol = frozenset(self.g.BFS_l(temp_parent, 1)) & in_solution_Vs
            Vs_of_interest__in_solution_Vs_map[temp_parent] = common_Vs_in_sol
            """
            if common_Vs_in_sol not in in_solution_Vs_groups:
                in_solution_Vs_groups[common_Vs_in_sol] = [temp_parent]
            else:
                in_solution_Vs_groups[common_Vs_in_sol].append(temp_parent)
            """
            if temp_parent not in parents__common_Vs_of_interest_map:
                temp_group = frozenset(self.combined_branches[temp_parent].V_in_branch) | {temp_parent}
                common_Vs_of_interest = temp_group & Vs_of_interest
                for temp_temp_parent in common_Vs_of_interest:
                    parents__common_Vs_of_interest_map[temp_temp_parent] = common_Vs_of_interest
                """
                temp_group_ND = frozenset(tgp_no_dummies.combined_branches[temp_parent].V_in_branch) | {temp_parent}
                common_Vs_of_interest_ND = temp_group_ND & Vs_of_interest
                groups[common_Vs_of_interest_ND] = temp_group_ND
                """
        groups4estimates = {}
        shared_Vs_in_sol = []
        for v in self._in_solution_Vs:
            for temp_parent in self.combined_branches[v].parents:
                common_Vs_in_sol = frozenset()
                common_Vs_of_interest = parents__common_Vs_of_interest_map[temp_parent]
                for temp_temp_parent in common_Vs_of_interest:
                    common_Vs_in_sol |= Vs_of_interest__in_solution_Vs_map[temp_temp_parent]
                if common_Vs_in_sol not in groups4estimates:
                    groups4estimates[common_Vs_in_sol] = common_Vs_of_interest
                else:
                    groups4estimates[common_Vs_in_sol] |= common_Vs_of_interest
                if len(common_Vs_in_sol) != 1:
                    shared_Vs_in_sol.append(common_Vs_in_sol)
        self.parents__common_Vs_of_interest_map = parents__common_Vs_of_interest_map
        self.Vs_of_interest__in_solution_Vs_map = Vs_of_interest__in_solution_Vs_map
        #self.groups = groups
        self.shared_Vs_in_sol = shared_Vs_in_sol
        self.groups4estimates = groups4estimates
        #assert self.check_separate_groups(groups.values())################## this should be removed, debugging
        return self.groups4estimates

    def get_separate_groups_old(self, **kwargs):
        #return
        groups = {}
        groups_Vs_of_interest = {}
        Vs_of_interest = frozenset(self.Vs_of_interest)
        parents_done = set()
        for v in self.in_solution_Vs:
            if not self.combined_branches[v].V_in_branch:
                continue
            grp = frozenset(self.combined_branches[v].V_in_branch)
            common_Vs_of_interest = grp & Vs_of_interest
            N_common_Vs_of_interest = len(common_Vs_of_interest)
            if N_common_Vs_of_interest == 1:
                if grp == common_Vs_of_interest:
                    if grp not in groups:
                        groups[grp] = [v]
                        groups_Vs_of_interest[grp] = common_Vs_of_interest
                    else:
                        groups[grp].append(v)
                        groups_Vs_of_interest[grp] = groups_Vs_of_interest[grp] | common_Vs_of_interest
                else:
                    groups[grp] = [v]
                    groups_Vs_of_interest[grp] = common_Vs_of_interest
            else:
                for temp_parent in self.combined_branches[v].parents:
                    if temp_parent in parents_done:
                        continue
                    temp_group = frozenset(self.combined_branches[temp_parent].V_in_branch)
                    other_parents = temp_group & common_Vs_of_interest
                    temp_group = temp_group | {temp_parent}
                    if not other_parents:
                        if temp_group not in groups:
                            groups[temp_group] = [v]
                            groups_Vs_of_interest[temp_group] = common_Vs_of_interest
                        else:
                            groups[temp_group].append(v)
                            groups_Vs_of_interest[temp_group] = groups_Vs_of_interest[temp_group] | common_Vs_of_interest
                    else:
                        groups[temp_group] = [v]
                        #groups[temp_group].extend(other_parents) # missing the roots actually
                        groups_Vs_of_interest[temp_group] = common_Vs_of_interest
                    for other_parent in other_parents:
                        other_temp_group = frozenset(self.combined_branches[other_parent].V_in_branch) | {other_parent}
                        assert temp_group == other_temp_group
                    parents_done.add(temp_parent)
                    assert not parents_done & other_parents
                    parents_done = parents_done | other_parents
        assert self.check_separate_groups(groups)################## this should be removes, debugging
        self.groups = groups
        self.groups_Vs_of_interest = groups_Vs_of_interest
        return groups

    def check_separate_groups(self, groups, **kwargs):
        tgp = self.__get_tgp_no_dummies()
        Vs_all1 = set()
        for v in self.in_solution_Vs:
            Vs_all1 = Vs_all1 | set(tgp.combined_branches[v].V_in_branch)
        if not groups:
            assert not Vs_all1
            return True
        Vs_all2 = set()
        N = len(groups)
        groups = list(groups)
        for i in range(N - 1):
            grp1= groups[i]
            Vs_all2 = Vs_all2 | grp1
            for j in range(i+1, N):
                grp2 = groups[j]
                assert not grp1 & grp2
        Vs_all2 = Vs_all2 | groups[-1]
        assert Vs_all1 == Vs_all2
        assert Vs_all1 | set(self.in_solution_Vs) == set(self.g.adj)
        return True

    ####################################################################################################################
    ####################################################################################################################

    # dummy estimates
    def get_dummy_match_estimates_in_solution(self, **kwargs):
        self.dummy_estimates = {}
        for v in self.Vs_of_interest:
            anchor_atoms = set(self.g.adj[v]) & self._in_solution_Vs
            if len(anchor_atoms) ==  1:
                anchor_atom = anchor_atoms.pop()
                max_Nv_Vother_branches = 0
                for v_other in self.g.adj[anchor_atom]:
                    if v_other==v or v_other not in self.Vs_of_interest:continue
                    if v not in self.g.rings_map:
                        Nv4comparison = self.combined_branches[v_other].Nv + 1
                    elif v_other not in self.g.rings_map:
                        Nv4comparison = self.combined_branches[v_other].Nv + 1
                    else:
                        common_ring = set(self.g.rings_map[v]) & set(self.g.rings_map[v_other])
                        if common_ring:
                            common_ring = common_ring.pop()
                            if len(set(common_ring.adj) & self._in_solution_Vs) == 1:
                                assert (set(common_ring.adj) & self._in_solution_Vs).pop() == anchor_atom
                                if anchor_atom in common_ring.breaks_map:
                                    if v_other in common_ring.breaks_map[anchor_atom]:
                                        for break_branch in self.break_branches[common_ring][common_ring.breaks_map[anchor_atom]]:
                                            if v not in break_branch.parents:
                                                Nv4comparison = break_branch.Nv + \
                                                            self.combined_FromRing_branches[common_ring][v_other].Nv
                                    else:
                                        Nv4comparison = self.combined_FromRing_branches[common_ring][v_other].Nv
                                else:
                                    Nv4comparison = self.combined_FromRing_branches[common_ring][v_other].Nv
                            else:
                                Nv4comparison = 0
                        else:
                            Nv4comparison = self.combined_branches[v_other].Nv + 1
                    if max_Nv_Vother_branches < Nv4comparison:
                        max_Nv_Vother_branches = Nv4comparison
                self.dummy_estimates[v] = max_Nv_Vother_branches
            else:
                self.dummy_estimates[v] = 0

    def get_dummy_match_estimates(self):
        self.sorted_combined_branches = self.get_sorted_combined_branches()
        self.dummy_estimates = {}
        self.sorted_RLCBs = {}
        # sorted_RLCBs is sorted_ring_linear_combined_branches, it's called when needed
        for v in self.Vs_of_interest:
#            if v in self.non_available_Vs:continue # is this needed at all?
            self.dummy_estimates[v] = self.get_v_dummy_match_estimate(v)

    def get_v_dummy_match_estimate(self, v_dummy):
        max_Nv_neigh_branches = 0
        max_Nv_combB_vDummy = None
        for temp_combined_branch in self.sorted_combined_branches:
            if temp_combined_branch.Nv < max_Nv_neigh_branches:
                return max_Nv_neigh_branches
            if v_dummy not in temp_combined_branch.V_in_branch and v_dummy not in temp_combined_branch.roots:
                return temp_combined_branch.Nv + 1 # this is already max(temp_branch.Nv, temp_max_Nv), if above
            elif max_Nv_combB_vDummy is None:
                max_Nv_combB_vDummy = self.get_max_Nv_combB_vDummy(v_dummy)
                max_Nv_neigh_branches = max(max_Nv_neigh_branches, max_Nv_combB_vDummy)
        return max_Nv_neigh_branches

    def get_max_Nv_combB_vDummy(self, v):
        max_Nv_neigh_branches = 0
        for r in self.ring_branches[v]:
            sorted_RLCB = self.get_ring_sorted_RLCB(r)
            Vs_independent = []
            for temp_ind_ring_part in r.independent_parts:
                if v in temp_ind_ring_part: # finding independent ring parts with this atom
                    for v1 in temp_ind_ring_part: # looping over each atom of the independent part
                        if v1 in self.available_Vs:
                            if v1 in r.breaks_map:
                                if v not in r.breaks_map[v1]:
                                    for break_branch in self.break_branches[r][r.breaks_map[v1]]:
                                        if break_branch.parents and break_branch.parents[0] not in temp_ind_ring_part:
                                            if max_Nv_neigh_branches < break_branch.tot_Nv:
                                                max_Nv_neigh_branches = break_branch.tot_Nv
                            if v1 != v:
                                if self.ring_independent_parts_availability_map[r][temp_ind_ring_part]:
                                    Vs_independent.append(v1)
                                elif v1 in r.adj[v]:
                                    Vs_independent.append(v1)
            temp_max_Nv = self._r_get_max_Nv_ring_linear_combined_branch(r, v, Vs_independent, sorted_RLCB,
                                                                       self.combined_FromRing_branches)
            if max_Nv_neigh_branches < temp_max_Nv:
                max_Nv_neigh_branches = temp_max_Nv
        for b in self.branches[v].values():
            if b.Nv > max_Nv_neigh_branches:
                max_Nv_neigh_branches = b.Nv
        return max_Nv_neigh_branches

    ####################################################################################################################
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    # helper functions # helper functions # helper functions # helper functions # helper functions # helper functions ##
    ####################################################################################################################

    def get_sorted_combined_branches(self):
        #temp_list = [self.combined_branches[v] for v in self.available_Vs]
        temp_list = [self.combined_branches[v] for v in self.Vs_of_interest]
        temp_list.sort(key = lambda x:x.Nv, reverse=True)
        return temp_list

    def get_ring_sorted_RLCB(self, r):
        if r not in self.sorted_RLCBs:
            self.sorted_RLCBs[r] = self._r_get_sorted_ring_linear_combined_branches(r, self.ring_branches, self.available_Vs)
        return self.sorted_RLCBs[r]

    @staticmethod
    def _r_get_sorted_ring_linear_combined_branches(r, ring_branches, available_Vs, Vs_subset = None):
        if Vs_subset is None:
            Vs_subset = r.adj
        temp_list = [ring_branches[v][r].combined_linear_branch for v in Vs_subset if v in available_Vs]
        temp_list.sort(key = lambda x:x.tot_Nv, reverse=True)
        return temp_list

    @staticmethod
    def _r_get_max_Nv_ring_linear_combined_branch(r, v_dummy, Vs_independent, sorted_RLCB, combined_FromRing_branches):
        # sorted_RLCB is sorted_ring_linear_combined_branches
        temp_Nv_max = 0
        for b in sorted_RLCB:
            # b.roots[0] is the root vertex
            # and we are looking for its linear_combined + FromRing (tot_Nv)
            # if v_dummy is its parents, we have to recalculate it
            if b.roots[0] in Vs_independent:
                if temp_Nv_max >= b.tot_Nv:
                    return temp_Nv_max
                if v_dummy not in b.parents: # this has to be done on the level of rings and with this if, since a vertex can be involved in more ind parts
                    return b.tot_Nv # this is already max(b.Nv, temp_Nv_max)
                else:
                    Nv_max = 0
                    for v1 in b.parents:
                        if v1!=v_dummy and Nv_max < combined_FromRing_branches[r][v1].Nv:
                            Nv_max = combined_FromRing_branches[r][v1].Nv
                    Nv_max += combined_FromRing_branches[r][b.roots[0]].Nv
                    temp_Nv_max = max(temp_Nv_max, Nv_max)
        return temp_Nv_max

    ####################################################################################################################
    ####################################################################################################################
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #
    # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs # allowing for non__available_Vs #

    """
    @staticmethod
    def _g_get_rings_available_Vs(g, available_Vs = None, non_available_Vs = None, flag_root_at=False):
        if available_Vs:
            available_Vs = list(available_Vs)
        elif non_available_Vs:
            available_Vs = list(set(g.adj) - set(non_available_Vs))
        else:
            available_Vs = list(g.adj)
        sub_graph_rings = []
        sub_graphs = []
        while available_Vs:
            v = available_Vs.pop()
            rings = []
            sub_graph_rings.append(rings)
            sub_graphs.append([v])
            path = [v]
            stack = []
            for i in g.adj[v]:
                stack.append((i, 1))
            while stack:
                cur, l = stack.pop()
                path = path[:l]
                if cur in available_Vs:
                    sub_graphs[-1].append(cur)
                    path.append(cur)
                    available_Vs.remove(cur)
                    for i in g.adj[cur]:
                        if i != path[-2]:
                            stack.append((i, l + 1))
                elif cur in path:
                    new_ring = []
                    for j in reversed(path):
                        for r_i in range(len(rings)):
                            if j in rings[r_i]:
                                if flag_root_at and len(set(rings[r_i]) & set(path)) < 2:
                                    continue
                                rings[r_i].extend(new_ring)
                                new_ring = rings.pop(r_i)
                                break
                        if j not in new_ring:
                            new_ring.append(j)
                        if j == cur:
                            break
                    rings.append(new_ring)
        return sub_graphs, sub_graph_rings
    """

class AlchemicalSolution_base:
    def __init__(self, tops, sol, available_atoms, flag_top_not_in_sol, sol_atoms_map, tried_pairs = None):
        self.tops = tops
        self._sol = sol
        self.available_atoms = available_atoms
        self.flag_top_not_in_sol = flag_top_not_in_sol
        self.sol_atoms_map = sol_atoms_map
        self.tried_pairs = tried_pairs

    ####################################################################################################################
    #                                    updates and copying
    ####################################################################################################################
    def copy(self):
        new_sol = self.copy4update()
        new_sol.tried_pairs = copy.copy(self.tried_pairs)
        self.update_copy_sol_toptp(new_sol)
        return new_sol

    def copy_add_state(self):
        new_sol = self.copy4update()
        new_sol.tried_pairs = copy.copy(self.tried_pairs)
        self.update_copy_sol_toptp_add_state(new_sol)
        new_sol._sol[new_sol._sol.columns.max() + 1] = None
        new_sol.flag_top_not_in_sol.append(True)
        return new_sol

    def copy4update(self):
        return self.copy4update_new_sol_df(self._sol.copy())

    def copy4update_new_sol_df(self, new_sol_df):
        return AlchemicalSolution_base(self.tops, new_sol_df, list(self.available_atoms),
                                          list(self.flag_top_not_in_sol), self.sol_atoms_map.copy(), list())

    def __update_copy_sol_toptp(self, new_sol):
        new_sol.toptp = self.toptp.reduce()
        new_sol.toptp.ptp_atom = self.toptp.ptp_atom.copy()
        for at in new_sol.toptp.get_atoms():
            for k in at.ND_m_states:
                at.ND_m_states[k] = set(at.ND_m_states[k])
            for k in at.ND_pch_states:
                at.ND_pch_states[k] = set(at.ND_pch_states[k])
            for k in at.ND_a_type_states:
                at.ND_a_type_states[k] = set(at.ND_a_type_states[k])
        new_sol.toptp.ptp_int = self.toptp.ptp_int.copy()
        new_sol.toptp.ptp_make_break_int = self.toptp.ptp_make_break_int.copy()
        new_sol.toptp.ptp_excl_pair = self.toptp.ptp_excl_pair
        new_sol.toptp.excl_pair = dict()
        #print(len(self.toptp.excl_pair))
        for atom_pair in self.toptp.excl_pair:
            new_atom_pair = frozenset(new_sol.toptp.atoms[at.id] for at in atom_pair)
            new_sol.toptp.excl_pair[new_atom_pair] = self.toptp.excl_pair[atom_pair]
            at1, at2 = new_atom_pair
            new_sol.toptp.add_atom_pair2EP_l(at1, at2)

    def __copy_interaction_states_map(self, flag_add_new_state=False):
        N_states = len(self.toptp.interaction_states_map)
        if flag_add_new_state:
            N_states+=1
        new_interaction_states_map = tuple(dict() for _ in range(N_states))
        done_interactions = set()
        for top_i, int_st_map in enumerate(self.toptp.interaction_states_map):
            for top_interaction in int_st_map: # interaction (e.g. bond between 2 atoms in one of the states/tops)
                sol_int = int_st_map[top_interaction] # solution interaction that stores ND_states (non-redundant) and int_states (interactions for each top)
                if sol_int not in done_interactions:
                    new_sol_int = top_interaction.__class__(top_interaction.int_type)
                    done_interactions.add(sol_int)
                    new_sol_int.ND_states={}
                    new_sol_int.int_states = [None] * N_states
                    for ND_state, top_set in sol_int.ND_states.items():
                        new_sol_int.ND_states[ND_state] = set(top_set)
                        for temp_top_i in top_set:
                            temp_int_top_i = sol_int.int_states[temp_top_i]
                            new_interaction_states_map[temp_top_i][temp_int_top_i] = new_sol_int
                            new_sol_int.int_states[temp_top_i] = temp_int_top_i
        return new_interaction_states_map

    def update_copy_sol_toptp(self, new_sol):
        if hasattr(self, 'toptp'):
            self.__update_copy_sol_toptp(new_sol)
            #new_sol.toptp.interaction_states_map = tuple(int_st_map.copy() for int_st_map in self.toptp.interaction_states_map)
            new_sol.toptp.interaction_states_map = self.__copy_interaction_states_map()

    def update_copy_sol_toptp_add_state(self, new_sol):
        if hasattr(self, 'toptp'):
            self.__update_copy_sol_toptp(new_sol)
            new_sol.toptp.interaction_states_map = self.__copy_interaction_states_map(flag_add_new_state=True)

    """
    def __add_update_match2sol_atom_None_pair(self, matched_pair):
        new_row = pd.DataFrame([[None] * self._sol.shape[1]], dtype = object)
        if self._sol.shape[0]:
            new_row.index = [self._sol.index[-1] + 1]
        new_row[matched_pair[0][0]] = matched_pair[0][1]
        new_sol_df = self._sol.append(new_row)
        new_sol = AlchemicalSolution_base(self.tops, new_sol_df, list(self.available_atoms),
                                          list(self.flag_top_not_in_sol), self.sol_atoms_map.copy(), list())
    """

    def __add_update_match2sol_atom_atom_pair_with_Dummy_None(self, matched_pair):
        new_row = pd.DataFrame([[None] * self._sol.shape[1]], dtype = object)
        if self._sol.shape[0]:
            new_row.index = [self._sol.index[-1] + 1]
        for i in range(2):
            if matched_pair[i][0] is not None:
                new_row[matched_pair[i][0]] = matched_pair[i][1]
        #new_sol_df = self._sol.append(new_row)
        new_sol_df = pd.concat([self._sol, new_row])
        new_sol = self.copy4update_new_sol_df(new_sol_df)
        for i in range(2):
            if matched_pair[i][1] not in (Dummy, None):
                top_available_atoms = new_sol.available_atoms[matched_pair[i][0]]
                pos = top_available_atoms.index(matched_pair[i][1])
                new_sol.available_atoms[matched_pair[i][0]] = top_available_atoms[:pos] + top_available_atoms[pos + 1:]
                new_sol.flag_top_not_in_sol[matched_pair[i][0]] = False
                new_sol.sol_atoms_map[matched_pair[i]] = new_row.index[0]
        new_row_pair_identifier = (None, new_row.index[0])
        if matched_pair[1][1] is None:
            new_sol.tried_pairs = list(self.tried_pairs)
        else:
            for temp_tried_pair in self.tried_pairs:
                flag = True
                for i in range(2):
                    if temp_tried_pair[i] in matched_pair and temp_tried_pair[i][1]!=Dummy:
                        updated_tried_pair = (new_row_pair_identifier, temp_tried_pair[(i+1) % 2]) # the other in pair
                        if updated_tried_pair[1][0] is None:
                            updated_tried_pair = tuple(sorted(updated_tried_pair))
                        if updated_tried_pair not in new_sol.tried_pairs:
                            new_sol.tried_pairs.append(updated_tried_pair)
                        flag = False
                        break
                if flag:
                    new_sol.tried_pairs.append(temp_tried_pair)
        return new_sol

    def __add_update_match2sol_row_atom_pair(self, matched_pair):
        new_sol = self.copy4update()
        new_sol._sol[matched_pair[1][0]][matched_pair[0][1]] = matched_pair[1][1]
        if matched_pair[1][1]!=Dummy:
            top_available_atoms = new_sol.available_atoms[matched_pair[1][0]]
            pos = top_available_atoms.index(matched_pair[1][1])
            new_sol.available_atoms[matched_pair[1][0]] = top_available_atoms[:pos] + top_available_atoms[pos + 1:]
            new_sol.flag_top_not_in_sol[matched_pair[1][0]] = False
            new_sol.sol_atoms_map[matched_pair[1]] = matched_pair[0][1]
        if matched_pair[1][1] != Dummy:
            row_identifier = (None, matched_pair[0][1])
            for temp_tried_pair in self.tried_pairs:
                flag = True
                for i in range(2):
                    if temp_tried_pair[i] == matched_pair[1]:
                        updated_tried_pair = (row_identifier, temp_tried_pair[(i+1) % 2]) # the other in pair
                        if updated_tried_pair[1][0] is None:
                            updated_tried_pair = tuple(sorted(updated_tried_pair))
                        if updated_tried_pair not in new_sol.tried_pairs:
                            new_sol.tried_pairs.append(updated_tried_pair)
                        flag = False
                        break
                if flag:
                    if matched_pair[0] in temp_tried_pair:
                        if temp_tried_pair not in new_sol.tried_pairs:
                            new_sol.tried_pairs.append(temp_tried_pair)
                    else:
                        new_sol.tried_pairs.append(temp_tried_pair)
        else:
            new_sol.tried_pairs = list(self.tried_pairs)
        return new_sol

    @staticmethod
    def __update_row(row, rows):
        """
        if row < rows[1]:
            return row
        elif row == rows[1]:
            return rows[0]
        else:
            return row - 1
        """
        if row == rows[1]:
            return rows[0]
        return row

    def __updated_row_row_tried_pair(self, tried_pair, rows):
        if tried_pair[0][0] != None:
            return tried_pair
        else:
            row_0 = self.__update_row(tried_pair[0][1], rows)
            if tried_pair[1][0] != None:
                return (None, row_0), tried_pair[1]
            else:
                row_1 = self.__update_row(tried_pair[1][1], rows)
                updated_pair = tuple(sorted(((None, row_0), (None, row_1))))
                return updated_pair

    def __add_update_match2sol_row_row_pair(self, matched_pair):
        new_sol = self.copy4update()
        rows = matched_pair[0][1], matched_pair[1][1] # has to be sorted!!!!!!!!!!!!
        row4merging = new_sol._sol.loc[rows[1]]
        row4merging = row4merging[(row4merging.values != None)]
        new_sol._sol.loc[rows[0], row4merging.index] = row4merging.values
        new_sol._sol.drop(rows[1], inplace = True)
        for temp_top_i_atom in zip(row4merging.index, row4merging):
            if temp_top_i_atom[1] != Dummy:
                assert new_sol.sol_atoms_map[temp_top_i_atom] == matched_pair[1][1]
                # update the sol_atoms_map
                new_sol.sol_atoms_map[temp_top_i_atom] = matched_pair[0][1]
        for temp_tried_pair in self.tried_pairs:
            updated_tried_pair = self.__updated_row_row_tried_pair(temp_tried_pair, rows)
            flag = True
            for i in range(2):
                if matched_pair[i] in temp_tried_pair:
                    if updated_tried_pair not in new_sol.tried_pairs:
                        new_sol.tried_pairs.append(updated_tried_pair)
                    flag = False
                    break
            if flag:
                new_sol.tried_pairs.append(updated_tried_pair)
        return new_sol

    """
    def __add_update_match2sol_row_row_pair_with_reset(self, matched_pair):
        new_sol = AlchemicalSolution_base(self.tops, self._sol.copy(), list(self.available_atoms),
                                          list(self.flag_top_not_in_sol), self.sol_atoms_map.copy(), list())
        rows = matched_pair[0][1], matched_pair[1][1] # has to be sorted!!!!!!!!!!!!
        row4merging = new_sol._sol.iloc[rows[1]]
        row4merging = row4merging[(row4merging.values != None)]
        new_sol._sol.iloc[rows[0], row4merging.index] = row4merging.values
        new_sol._sol.drop(rows[1], inplace = True)
        for temp_top_i_atom in zip(row4merging.index, row4merging):
            if temp_top_i_atom[1] != Dummy:
                assert new_sol.sol_atoms_map[temp_top_i_atom] == matched_pair[1][1]
            #new_sol.sol_atoms_map[temp_top_i_atom] = matched_pair[0][1]
            # this is in principle the update if no reset is done...
        new_sol._sol.reset_index(drop = True, inplace = True)
        new_sol._get_sol_atoms_map()
        for temp_tried_pair in self.tried_pairs:
            updated_tried_pair = self.__updated_row_row_tried_pair(temp_tried_pair, rows)
            flag = True
            for i in range(2):
                if temp_tried_pair[i] == matched_pair[1]:
                    # if things get done without iloc and reseting the index, update has to be here...
                    if updated_tried_pair not in new_sol.tried_pairs:
                        new_sol.tried_pairs.append(updated_tried_pair)
                    flag = False
                    break
            if flag:
                new_sol.tried_pairs.append(updated_tried_pair)
        return new_sol
    """

    def add_match2sol(self, matched_pair):
        """
        adds a matched pair (of atoms) to the current solution and creates a new one
        updates the other important variables (all from self.__init__ and tried_pairs)
        :param matched_pair:
        :return:
        """
        # matched_pair is of this format (top_index, atom_key), (top_index, atom_key) -> sorted by top_index
        # or (None, row_index), (top_index, atom_key) -> always None, row_index first
        # or (None, row_index), (None, row_index) -> row index not sorted
        if matched_pair[0][0] is not None:
            return self.__add_update_match2sol_atom_atom_pair_with_Dummy_None(matched_pair)
        elif matched_pair[1][0] is not None:
            return self.__add_update_match2sol_row_atom_pair(matched_pair)
        else:
            return self.__add_update_match2sol_row_row_pair(matched_pair)
#        return new_solution

    ####################################################################################################################
    #                                           END
    #                                    updates and copying
    ####################################################################################################################

    def _get_sol_atoms_map(self):
        self.sol_atoms_map = {}
        for row_i, row in self._sol.iterrows():
            for top_i in self._sol.columns:
                if row[top_i] not in (None, Dummy):
                    self.sol_atoms_map[(top_i, row[top_i])] = row_i

    def transform_top_index_atom(self, top_index_atom):
        if top_index_atom in self.sol_atoms_map:
            return None, self.sol_atoms_map[top_index_atom]
        else:
            return top_index_atom

    def transform_top_index_atom_pair(self, matched_pair):
        transformed_matched_pair = tuple(self.transform_top_index_atom(temp_top_i_at) for temp_top_i_at in matched_pair)
        if transformed_matched_pair[1][0] is None:
            if transformed_matched_pair[0][0] is None:
                transformed_matched_pair = tuple(sorted(transformed_matched_pair))
            else:
                transformed_matched_pair = transformed_matched_pair[::-1]
        return transformed_matched_pair

    def update_tried_pairs(self):
        updated_tried_pairs = set()
        for temp_tried_pair in self.tried_pairs:
            updated_tried_pairs.add(self.transform_top_index_atom_pair(temp_tried_pair))
        self.tried_pairs = list(updated_tried_pairs)

    def states2df(self, **kwargs):
        temp_df = self._sol.copy()
        temp_df.columns = pd.MultiIndex.from_product([['top'], temp_df.columns])
        temp_df.index = pd.MultiIndex.from_product([['atom'], temp_df.index])
        if kwargs.get('mask_dummies', True):
            temp_df = temp_df.mask(temp_df.values == Dummy, 'DUM')
        return temp_df

    def __repr__(self):
        df = self.states2df()
        return str(df)

    def find_common_atom_pair(self, top_atom0_index, top1_index):
        row = self.sol_atoms_map[top_atom0_index]
        return self._sol.loc[row, top1_index]

    def find_sol_atom(self, top_atom0_index):
        row = self.sol_atoms_map[top_atom0_index]
        return self.toptp.atoms[row]

class AlchemicalSolution(AlchemicalSolution_base):
    def __init__(self, tops, **kwargs):
        """
        :param tops:
        :param kwargs:
            available_atoms
                ((n-atoms of 1st top), (m-atoms of 2nd top), (l-atoms of 3rd top), ...)
                e.g. (None, (1,2,3,4,6,8), (), None)
                    None - all atoms
                    4 topologies, all atoms from the 1st and 4th topology, (1,2,3,4,6,8) from 2nd, none from 3rd
                    () - no atoms (needed for divide-and-conquer implementation
            common_atoms
                ((n-atoms of 1st top), (n-atoms of 2nd top), (n-atoms of 3rd top), ...)
                e.g. [[1,3,8], [1,4,8], [1, Dummy, 7], [2, 4, None]]
                    4 topologies, atoms (1, 1, 1, 2); (3, 4, dummy, 4) and (8, 8, 7, ?) are common
                    Dummy means dummy, None means not decided yet
            tried_pairs
                [((top_index, atom), (top_index, atom)), ((top_index, atom), (top_index, atom)), ...]
                e.g. [((0,1), (1,1)), ((0,1), (2,2)), ((1,2), (2,3))]
                three pairs - tops[0].atoms[1] and tops[1].atoms[1]; tops[0].atoms[1] and tops[2].atoms[2], ...
        """
        self.tops = tops
        common_atoms = kwargs.get('common_atoms')
        if common_atoms:
            assert len(common_atoms) == len(tops)
            self._sol = pd.DataFrame(common_atoms, dtype = object).T
        else:
            self._sol = pd.DataFrame(np.empty((0, len(tops)), dtype = object))
        self._get_sol_atoms_map()
        available_atoms = kwargs.get('available_atoms')
        if available_atoms:
            assert len(available_atoms) == len(tops)
        else:
            available_atoms = [None] * len(tops)
        self.available_atoms = []
        for top_i, top_av_atoms in enumerate(available_atoms):
            if available_atoms[top_i] is None:
                temp_set = set(self.tops[top_i].adj)
            else:
                temp_set = set(available_atoms[top_i])
            temp_atoms_in_solution = set(self._sol[top_i])
            self.available_atoms.append(tuple(temp_set - temp_atoms_in_solution))
        self.__get_flag_top_not_in_sol()
        self.tried_pairs = kwargs.get('tried_pairs', [])
        self.update_tried_pairs()
        self.matches_layers = [[] for _ in range(kwargs.get('N_layers', 3))]

    def __get_flag_top_not_in_sol(self): # shape = len(tops); True if not a single atom in solution yet
        if self._sol.shape[0]:
            self.flag_top_not_in_sol = [self._sol[top_i].isin([None, Dummy]).all() for top_i in self._sol.columns]
        else:
            self.flag_top_not_in_sol = [True] * self._sol.shape[1] #

"""
G - graph defined as a dictionary (G.adj)
v - vertex
e - edges
"""

from SMArt.incl import defaultdict, __VersionCompatibility, np, pd, minimize, plt, deque
from SMArt.geometry import generate_new_coordinates


class Graph(__VersionCompatibility):
    """graph defined as a dictionary (adjacency list)"""
    def __init__(self, adj_dict_edges=None, flag_edges=False, adj_type=defaultdict, adj_format_new_vertex=set, **kwargs):
        """
        :param adj_dict_edges:
        :param flag_edges:
        :param adj_type:
        :param adj_format_new_vertex:
        :param kwargs:
             set_adj_type - change the type of self.adj - dict is default
             set_adj_format - change the type of each element in the self.adj - tuple is default
             fnc2add - used for self._set_adj_type_format
        """
        self._set_adj_type_format(adj_type, adj_format_new_vertex, **kwargs)
        if adj_dict_edges:
            if flag_edges:
                self.edges2adj(adj_dict_edges)
            else:
                self.adj.update(adj_dict_edges)
        kwargs['set_adj_type'] = kwargs.get('set_adj_type', dict)
        kwargs['set_adj_format'] = kwargs.get('set_adj_format', tuple)
        self.set_adj_type_format(**kwargs)

    def _set_adj_type_format(self, adj_type = None, adj_format_new_vertex = None, **kwargs):
        """
        used in principle only at the initialization - sets the type of the self.adj and type for new elements
        :param adj_type:
        :param adj_format_new_vertex:
        :param kwargs:
            fnc2add - if adj_format_new_vertex is custom-made class, fnc2add would be a function that adds new elements
        :return:
        """
        if adj_format_new_vertex:
            self.__adj_format_new_vertex = adj_format_new_vertex
            if adj_format_new_vertex not in (list, set):
                self.__fnc2add = kwargs.get('fnc2add')
                assert self.__fnc2add
        elif not hasattr(self, '_Graph__adj_format_new_vertex'):
            self.__adj_format_new_vertex = set
        if adj_type:
            self.__adj_type = adj_type
            if not hasattr(self, 'adj'):
                self.adj = {}
            if adj_type is defaultdict:
                self.adj = defaultdict(self.__adj_format_new_vertex, self.adj)
            else:
                self.adj = adj_type(self.adj)
        if self.__adj_type == defaultdict and self.__adj_format_new_vertex == set:
            self.add_edge = self._add_edge_defaultdict_set
        elif self.__adj_format_new_vertex == set:
            self.add_edge = self._add_edge_set
        elif self.__adj_format_new_vertex == list:
            self.add_edge = self._add_edge_list
        else:
            self.add_edge = self._add_edge_general

    def _add_edge_defaultdict_set(self, edge):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            self.adj[edge[i]].add(v_another)

    def _add_edge_set(self, edge):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            if edge[i] not in self.adj:
                self.adj[edge[i]] = {v_another}
            else:
                self.adj[edge[i]].add(v_another)

    def _add_edge_list(self, edge):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            if edge[i] not in self.adj:
                self.adj[edge[i]] = [v_another]
            elif v_another not in self.adj[edge[i]]:
                self.adj[edge[i]].append(v_another)

    def _add_edge_general(self, edge, **kwargs):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            if edge[i] not in self.adj:
                self.adj[edge[i]] = self.__adj_format_new_vertex()
            if v_another not in self.adj[edge[i]]:
                getattr(self.adj[edge[i]], self.__fnc2add)(v_another)

    def set_adj_type_format(self, set_adj_type = None, set_adj_format = None, **kwargs):
        """
        :param set_adj_type: self.adj type
        :param set_adj_format:  type of each element in self.adj
        :param kwargs:
        :return:
        """
        if set_adj_type:
            self._set_adj_type_format(adj_type=set_adj_type)
        if set_adj_format:
            for v in self.adj:
                self.adj[v] = set_adj_format(self.adj[v])

    def check_G(self):
        """check if adjacency list is correct"""
        for v in self.adj:
            for v2 in self.adj[v]:
                if v not in self.adj[v2]: return False
        return True

    def get_edges(self):
        """from graph to edges - generator"""
        checked_v = []
        for v in self.adj:
            checked_v.append(v)
            for v_n in self.adj[v]:
                if v_n not in checked_v:
                    yield (v, v_n)

    def add_edge(self, edge):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            self.adj[edge[i]].add(v_another)

    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge)

    def remove_edge(self, edge):
        for i in range(2):
            v_another = edge[(i + 1) % 2]
            self.adj[edge[i]].remove(v_another)

    def edges2adj(self, edges):
        """from edges to adjacency dictionary"""
        self.add_edges(edges)

    def update_adj(self, additional_adj):
        tempG = Graph(additional_adj)
        for edge in tempG.get_edges():
            self.add_edge(edge)

    def sub_graph(self, vertices, **kwargs):
        """
        creates a subgraph
        :param vertices:
        :param kwargs:
            adj_fromat (set, tuple)
        :return:
        """
        new_adj = defaultdict(set)
        for v in vertices:
            new_adj[v]
            for at_neigh in self.adj[v]:
                if at_neigh in vertices:
                    new_adj[v].add(at_neigh)
        if kwargs.get('flag_directed'):
            return GraphDirected(new_adj, **kwargs)
        return Graph(new_adj, **kwargs)

    def __get_stack_visited(self, v, **kwargs):
        if 'stack_visited' in kwargs:
            return kwargs['stack_visited']
        if 'visited' in kwargs:
            visited = kwargs['visited']
        else:
            visited = set()
        if kwargs.get('flag_v_list'):
            if kwargs.get('flag_stack_level', True):
                stack = deque(zip(v, [0]*len(v)))
            else:
                stack = deque(v)
            if kwargs.get('visited_from_stack', True):
                visited |= set(v)
        else:
            if kwargs.get('flag_stack_level', True):
                stack = deque([(v, 0)])
            else:
                stack = deque([v])
            if kwargs.get('visited_from_stack', True):
                visited.add(v)
        return stack, visited

    def _BFS(self, v, **kwargs):
        """
        Breath-first search - generator (no levels)
        :param v:
        :param kwargs:
            flag_stack_level [True]
            flag_v_list
            stack_visited
            visited
            visited_from_stack [True]
        :return:
        """
        kwargs['flag_stack_level'] = False
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur = stack.popleft()
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append(i)
                    visited.add(i)
            yield cur

    def BFS(self, v, **kwargs):
        """
        Breath-first search - generator
        :param v:
        :param kwargs:
            flag_stack_level [True]
            flag_v_list
            stack_visited
            visited
            visited_from_stack [True]
        :return:
        """
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur, l = stack.popleft()
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append((i, l + 1))
                    visited.add(i)
            yield cur, l

    def BFS_d(self, v, d, **kwargs):
        """Breath-first search with depth - finds first d neighbours (all from 1 to d)
        :param v:
        :param kwargs:
            flag_stack_level [True]
            flag_v_list
            stack_visited
            visited
            visited_from_stack [True]
        """
        for cur, l in self.BFS(v, **kwargs):
            if l > d:
                break
            yield cur

    def BFS_l(self, v, level, **kwargs):
        """Breath-first search with level - finds all l-th neighbours (not including <l)
        :param v:
        :param kwargs:
            flag_stack_level [True]
            flag_v_list
            stack_visited
            visited
            visited_from_stack [True]
        """
        for cur, l in self.BFS(v, **kwargs):
            if l < level:
                continue
            if l == level:
                yield cur
            else:
                break

    def BFS_f(self, v, f, *args, **kwargs):
        """Breath-first search - with an on-the-fly-function

example:
def f(*args, **kwargs):
        cur = kwargs['BFS_local'][0]
        kwargs["all_v"].append(cur)
G.BFS_f(v,f,all_v=[])

returns ((), {'all_v':[v1,v2,v3,...]})
        """
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur, l = stack.popleft()
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append((i, l + 1))
                    visited.add(i)
            kwargs['BFS_local'] = (cur, l, stack, visited)
            f(*args, **kwargs)
        return args, kwargs

    def BFS_d_f(self, v, d, f, *args, **kwargs):
        """Breath-first search with depth (see BFS_d) - together with an on-the-fly-function

example:
def f(*args, **kwargs):
        cur = kwargs['BFS_local'][0]
        kwargs["all_v"].append(cur)
self.BFS_d_f(v,d,f,all_v=[])

returns ((), {'all_v':[v1,v2,v3,...]})"""
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur, l = stack.popleft()
            if l > d:
                break
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append((i, l + 1))
                    visited.add(i)
            kwargs['BFS_local'] = (cur, l, stack, visited)
            f(*args, **kwargs)
        return args, kwargs

    def dist(self, v1, v2):
        """distance between v1 and v2"""
        for cur, l in self.BFS(v1):
            if cur == v2:
                return l

    def DFS(self, v, **kwargs):
        """Depth-first search - generator"""
        kwargs['flag_stack_level'] = False
        kwargs['visited_from_stack'] = False
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur = stack.pop()
            if cur not in visited:
                visited.add(cur)
                for i in self.adj[cur]:
                    stack.append(i)
                yield cur

    def DFS_f(self, v, f, *args, **kwargs):
        """Depth-first search - with an on-the-fly-function

example:
def f(*args, **kwargs):
    cur = kwargs['DFS_local'][0]
    kwargs["all_v"].append(cur)
self.DFS_f(v,f,all_v=[])
returns ((), {'all_v':[v1,v2,v3,...]})
        """
        path = []
        kwargs['visited_from_stack'] = False
        stack, visited = self.__get_stack_visited(v, **kwargs)
        while stack:
            cur, l = stack.pop()
            if cur not in visited:
                path = path[:l]
                path.append(cur)
                visited.add(cur)
                for i in self.adj[cur]:
                    if i != path[-2]:
                        stack.append((i, l + 1))
                kwargs['DFS_local'] = (cur, l, stack, visited, path)
                f(*args, **kwargs)
        return args, kwargs

    def get_Gs(self):
        """if G is not connected, finds separate Gs"""
        Vs = set(self.adj.keys())
        while Vs:
            v = next(iter(Vs))
            temp_vs = []
            for v in self.DFS(v):
                temp_vs.append(v)
                Vs.remove(v)
            yield self.sub_graph(temp_vs)

    def _find_rings(self, flag_root_at=False):
        """finds all rings in a graph"""
        v = next(iter(self.adj.keys()))
        rings = []
        path = [v]
        stack, visited = self.__get_stack_visited(v)
        for i in self.adj[v]:
            stack.append((i, 1))
        while stack:
            cur, l = stack.pop()
            path = path[:l]
            if cur not in visited:
                path.append(cur)
                visited.add(cur)
                for i in self.adj[cur]:
                    if i != path[-2]:
                        stack.append((i, l + 1))
            elif cur in path:
                new_ring = []
                for j in reversed(path):
                    for r_i in range(len(rings)):
                        if j in rings[r_i]:
                            if flag_root_at:
                                if len(set(rings[r_i]) & set(path)) < 2:
                                    continue
                                if j==cur and len(set(rings[r_i]) & set(new_ring)) < 2:
                                    continue
                            rings[r_i].extend(new_ring)
                            new_ring = rings.pop(r_i)
                            break
                    if j not in new_ring:
                        new_ring.append(j)
                    if j == cur:
                        break
                rings.append(new_ring)
        if len(visited) == len(self.adj):
            return tuple(rings)
        return False

    def find_rings(self, flag_root_at=False, **kwargs):
        rings = []
        for temp_G in self.get_Gs():
            rings.extend(temp_G._find_rings(flag_root_at))
        return tuple(rings)

    def find_rings_graph(self, flag_root_at=False, **kwargs):
        rings = self.find_rings(flag_root_at, **kwargs)
        for r in rings:
            yield self.sub_graph(r)

    def get_2D_repr(self, vertices=None, d_12=1., d_13=1.7,  Fk=(100, 50, 5), rest_pow = 1, Nmaxiter=25, **kwargs):
        """
        :param vertices:
        :param kwargs: fig_name, step, init_struc
        :return:
        """
        fig_name = kwargs.get('fig_name')
        step = kwargs.get('step', False)
        init_struc = kwargs.get('init_struc', False)
        if vertices is None:
            vertices = list(self.adj)
        temp_G = self.sub_graph(vertices, flag_directed = True, parents = vertices[0])
        coord_df = pd.DataFrame(np.zeros((len(vertices), 2)), pd.Index(temp_G.adj.keys(), name = 'idx'))
        coord_df.columns = ['x', 'y']
        vertices = [vertices[0]]
        grp_add, added_paths, anchor_points = temp_G.get_add_missing_groups(**kwargs)
        for i in range(len(added_paths)):
            path2add, anchor_p = added_paths[i], anchor_points[i]
            generate_new_coordinates(path2add, anchor_p, coord_df, v_fact=1., fnc_ring_path=1)
            vertices += grp_add[i]
            sub_temp_G = temp_G.sub_graph(vertices)
            temp_coorf_df = coord_df.loc[vertices]
            if i + 1 != len(added_paths):
                temp_min = self.__do_min(sub_temp_G, temp_coorf_df, 16, Nmaxiter, d_12, d_13,  Fk, rest_pow)
            else:
                temp_min = self.__do_min(sub_temp_G, temp_coorf_df, 16, Nmaxiter, d_12, d_13,  Fk, rest_pow)
            if step and fig_name and init_struc:
                self.__plot_graph(temp_G, coord_df, fig_name + '_step_' + str(i) + '_init.png', **kwargs)
            coord_df.loc[vertices] = temp_min.x.reshape((len(vertices),2))
            if step and fig_name:
                self.__plot_graph(temp_G, coord_df, fig_name + '_step_' + str(i) + '_min.png', **kwargs)
        sub_temp_G = temp_G.sub_graph(vertices)
        temp_coorf_df = coord_df.loc[vertices]
        temp_min = self.__do_min(sub_temp_G, temp_coorf_df, 25, None, d_12, d_13, Fk, rest_pow)
        coord_df.loc[vertices] = temp_min.x.reshape((len(vertices),2))
        if fig_name:
            self.__plot_graph(temp_G, coord_df, fig_name + '_min.png', **kwargs)
        else:
            self.__plot_graph(temp_G, coord_df, **kwargs)
        return coord_df

    @staticmethod
    def __plot_graph(temp_G, coord_df, ffig = None, flag_name_attribute=None, **kwargs):
        x = coord_df.x.values
        y = coord_df.y.values
        plt.scatter(x,y,color = "k")
        v_list = list(coord_df.index)
        for i in range(len(v_list)):
            if flag_name_attribute:
                #plt.text(x[i], y[i], getattr(v_list[i], flag_name_attribute), ha='center', va='center')
                plt.text(x[i], y[i], getattr(v_list[i], flag_name_attribute))
            else:
                #plt.text(x[i], y[i], str(v_list[i]), ha='center', va='center')
                plt.text(x[i], y[i], str(v_list[i]))
        for b in temp_G.get_edges():
            temp_coord = coord_df.loc[list(b)]
            temp_x = temp_coord.x
            temp_y = temp_coord.y
            plt.plot(temp_x, temp_y, color='k')
        plt.xlim(min(x) - 1, max(x) + 1)
        plt.ylim(min(y) - 1, max(y) + 1)
        plt.gca().set_aspect('equal')
        if ffig is None:
            if kwargs.get('flag_show', True):
                plt.show()
        else:
            plt.savefig(ffig)
            plt.close()

    def __do_min(self, temp_G, coord_df, cutoff2, Nmaxiter, d_12, d_13,  Fk, rest_pow):
        l12, l13, lrest = self.__get_123(temp_G, coord_df, cutoff2)
        coord_flat = coord_df.values
        coord_flat = coord_flat.reshape(coord_flat.shape[0] * 2)
        if Nmaxiter is None:
            res = minimize(self.__E, coord_flat, args = (l12, l13, lrest, d_12, d_13,  Fk, rest_pow))
        else:
            res = minimize(self.__E, coord_flat, args = (l12, l13, lrest, d_12, d_13,  Fk, rest_pow),
                           options={'maxiter': Nmaxiter})
        return res

    @staticmethod
    def __get_123(temp_G, coord_df, cutoff2 = None):
        checked_pairs = []
        l12 = []
        l13 = []
        lrest = []
        v_list = list(coord_df.index)
        v_index_map = dict(zip(v_list, range(len(v_list))))
        for v1 in temp_G.adj:
            for (v2, l) in temp_G.BFS(v1):
                if v1 != v2:
                    temp_p = tuple(sorted((v_index_map[v1], v_index_map[v2])))
                    if temp_p not in checked_pairs:
                        checked_pairs.append(temp_p)
                        if l == 1:
                            l12.append(temp_p)
                        elif l == 2:
                            l13.append(temp_p)
                        else:
                            if cutoff2:
                                if np.sum(np.diff(coord_df.loc[[v1,v2]].values)**2) < cutoff2:
                                    lrest.append(temp_p)
                            else:
                                lrest.append(temp_p)
        return l12, l13, lrest

    @staticmethod
    def __E(coord_flat, l12, l13, lrest, d_12, d_13,  Fk, rest_pow):
        E = 0
        coord = coord_flat.reshape(coord_flat.shape[0] // 2, 2)
        temp_ind_lists = list(zip(*l12))
        temp_distances = coord[temp_ind_lists[0],] - coord[temp_ind_lists[1],]
        temp_distances = np.linalg.norm(temp_distances, axis = 1)
        E += 0.5 * Fk[0] * np.linalg.norm(temp_distances - np.array([d_12] * temp_distances.shape[0]))
        if l13:
            temp_ind_lists = list(zip(*l13))
            temp_distances = coord[temp_ind_lists[0],] - coord[temp_ind_lists[1],]
            temp_distances = np.linalg.norm(temp_distances, axis = 1)
            E += 0.5 * Fk[1] * np.linalg.norm(temp_distances - np.array([d_13] * temp_distances.shape[0]))
        if lrest:
            temp_ind_lists = list(zip(*lrest))
            temp_distances = coord[temp_ind_lists[0],] - coord[temp_ind_lists[1],]
            temp_distances = Fk[2] / np.sum(temp_distances**(2*rest_pow), axis = 1)
            E += np.sum(temp_distances)
        return E


class Vertex(__VersionCompatibility):
    __slots__ = ('p', 'c', 'l')

    def __init__(self, l, p=None):
        self.l = l # level
        if p is None:
            self.p = set()
        else:
            self.p = {p}
        self.c = set()

    def add_p(self, p):
        self.p.add(p)

    def add_c(self, c):
        self.c.add(c)


class GraphDirected(Graph):
    def __init__(self, adj_dict_edges = None, flag_edges = False, parents = None, **kwargs):
        """
        directed graph - it could have more than 1 vertex as roots (parents)
        :param adj_dict_edges:
        :param flag_edges:
        :param parents:
        :param kwargs:
        """
        super(GraphDirected, self).__init__(adj_dict_edges, flag_edges, **kwargs)
        if parents is not None:
            try:
                self.set_parents(*parents)
            except TypeError:
                self.set_parents(parents)
            self.get_p_c(**kwargs)

    def directed_BFS(self):
        return self.BFS(self.parent_v, flag_v_list = True)

    def set_parents(self, *v, **kwargs):
        self.parent_v = v # in principle this should be 1 vertex - but here it's generalized to many
        self.get_p_c(**kwargs)

    def get_p_c(self, **kwargs):
        """
        iterate over vertices and get parents and children
        :param kwargs:
            level - max level
            flag_children_level_map - if and how many levels to add to the map (default 0; -1 adds all levels)
        :return:
        """
        self.Gv = {}
        level_max = kwargs.get('level')
        visited = set()
        stack = deque()
        for vv in self.parent_v:
            self.Gv[vv] = Vertex(0)
            visited.add(vv)
            stack.append((vv, 0))
        while stack:
            cur, l = stack.popleft()
            if level_max is not None and l==level_max:
                break
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append((i, l + 1))
                    visited.add(i)
                    self.Gv[i] = Vertex(l + 1, cur)
                    self.Gv[cur].add_c(i)
                elif l+1 == self.Gv[i].l and cur not in self.Gv[i].p:
                    self.Gv[i].add_p(cur)
                    self.Gv[cur].add_c(i)

    def get_parents(self, v, path=None):
        if path is None:
            path = []
        path.append(v)
        if self.Gv[v].p:
            for p in self.Gv[v].p:
                self.get_parents(p, path)
        return path

    def find_shortest_rpath(self):
        visited = set()
        stack = deque()
        for vv in self.parent_v:
            visited.add(vv)
            stack.append((vv, 0))
        while stack:
            cur, l = stack.popleft()
            for i in self.adj[cur]:
                if i not in visited:
                    stack.append((i, l + 1))
                    visited.add(i)
#                    self.Gv[i] = Vertex(l + 1, cur)
                elif self.Gv[i].l >= l and self.Gv[i].l + l != 0:
                    if len(self.Gv[i].p) == 2:
                        p1, p2 = list(self.Gv[i].p)
                        path = list(reversed(self.get_parents(p1)))
                        path.append(i)
                        path.extend(self.get_parents(p2))
                    else:
                        path = list(reversed(self.get_parents(cur)))
                        path.extend(self.get_parents(i))
                    if path[0] == path[-1]:
                        path = path[:-1]
                    return path

    def __get_ap_atoms(self, anchor_point, set_parents_v = None):
        if set_parents_v is None:
            set_parents_v = set(self.parent_v)
        return set(self.adj[anchor_point]) & set_parents_v

    def __check_partial_ring(self, rings, rings_map, **kwargs):
        temp_rings = set()
        ring_atoms_in_parents = []
        for v in self.parent_v:
            if v in rings_map:
                temp_rings.add(rings_map[v])
                ring_atoms_in_parents.append(v)
        atoms_in_parents_complete_rings = []
        for r_index in temp_rings:
            atoms_in_parents_complete_rings += rings[r_index]
        return not(len(ring_atoms_in_parents) == len(atoms_in_parents_complete_rings))

    def __get_add_missing_groups(self, **kwargs):
        rings = self.find_rings()
        rings_map = {}
        for rc, r in enumerate(rings):
            for ring_vertex in r:
                rings_map[ring_vertex] = rc
        grp_add = []
        added_paths = []
        anchor_points = []
        while len(self.parent_v)!=len(self.adj):
            grp_add.append([])
            anchor_points.append({})
            added_paths.append({})
            # add rind shortes path
            flag_ring_path = None
            if self.__check_partial_ring(rings, rings_map):
                temp_path = self.find_shortest_rpath()
                temp_path2add = []
                temp_anch_points = []
                if temp_path:
                    flag_ring_path = rings_map[temp_path[0]]
                    for temp_v in temp_path:
                        if temp_v not in self.parent_v:
                            temp_path2add.append(temp_v)
                        else:
                            temp_anch_points.append(temp_v)
                    temp_anch_points = tuple(temp_anch_points)
                    if temp_anch_points not in added_paths[-1]:
                        added_paths[-1][temp_anch_points] = []
                    added_paths[-1][temp_anch_points].extend(temp_path2add)
                    grp_add[-1].extend(temp_path2add)
            # add other first neighbours
            for temp_child, level in self.directed_BFS():
                if level==0:continue
                if level == 2:break
                if temp_child not in grp_add[-1]:
                    if temp_child in rings_map and rings_map[temp_child]==flag_ring_path:continue
                    assert len(self.Gv[temp_child].p) == 1
                    temp_parent = list(self.Gv[temp_child].p)[0]
                    grp_add[-1].append(temp_child)
                    if (temp_parent,) not in added_paths[-1]:
                        added_paths[-1][(temp_parent,)] = []
                    added_paths[-1][(temp_parent,)].append(temp_child)
            temp_set_parents = set(self.parent_v)
            for anchor_point in added_paths[-1]:
                anchor_points[-1][anchor_point] = []
                for an_p in anchor_point:
                    anchor_points[-1][anchor_point].append(self.__get_ap_atoms(an_p, temp_set_parents))
            new_parents = temp_set_parents | set(grp_add[-1])
            self.set_parents(*new_parents)
        return grp_add, added_paths, anchor_points


    def get_add_missing_groups(self, **kwargs):
        return self.__get_add_missing_groups(**kwargs)

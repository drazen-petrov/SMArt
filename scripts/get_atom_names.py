import SMArt
from SMArt.md import Topology, parse_mtb, parse_top
from SMArt import alchemy
from SMArt.md.gro2gro import g2g
import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-qc'])

SMArt.incl.set_print_warnings(False)

def find_matching(in_file, ifp_file, bb_id, pdb_file, selection='all', flag_topology=False):
    """
    uses SMArt.alchemy to find possible matches between a pdb and a building block
    
    selection is a string defining a pymol selection (e.g. resn FMN)
    in_file is expecting a mtb file, however if flag_topology==True; it would assume that the in_file is a topo)
    also, if flag_topology==True, ifp_file is disregarded

    bb_id is the name of the building block, or residue number if a topology is given
    """

    # get the info from the 
    cmd.load(pdb_file)
    t = Topology()
    cmd.iterate(selection, 'a = t.add_atom(index, name, a_type=elem.upper())', space=dict(t=t, a=None))
    for a in t.atoms.values():
        a.element = a.a_type
    # get the bonds
    edges = []
    for b in cmd.get_bond('all', selection)[0][1]:
        edges.append((t.atoms[b[0]], t.atoms[b[1]]))
    G = SMArt.graph.Graph(edges, True)
    t.adj = G.adj

    # get the building block / residue
    if flag_topology:
        top = parse_top(in_file)
        atms = [at.id for at in top.residues[bb_id].atoms]
        bb_res=top.reduce_top(atms)
    else:
        mtb = parse_mtb(in_file, ifp_file = ifp_file)
        bb_res = mtb.bb[bb_id]
    bb_res.gen_graph()

    for a in bb_res.get_atoms():
        element = int(a.gr_get_element())
        element = g2g._AtomicNum__element_map[element]
        a.element = element

    # let's do the matching
    mcs = alchemy.MCS(t, bb_res)
    mcs.enumerate_stepwise_sorted()
    return mcs

def check_elements(mcs, sol_i = 0, sort_by_top = 0):
    """
    get only atom pairs that are have matching elements
    sol_i - solution index
    sort_by_top [0,1] - sort atoms according to the first/second topology (pdb or mtb/top)
    """
    s = mcs.solutions[sol_i]
    diff_same = [], []
    for i, row in s._sol.iterrows():
        if alchemy.incl.Dummy not in row.values:
            diff_same[row[0].element == row[1].element].append(i)
    match = s._sol.iloc[diff_same[1]]
    match['pos'] = [int(a.id) for a in match.iloc[:,sort_by_top]]
    match = match.sort_values('pos').drop('pos', axis=1).reset_index(drop=True)
    return [len(diff_same[i]) for i in range(2)], diff_same, match

if __name__ == '__main__':
    import argparse
    import os

    desc = 'run in an interactive mode: python -i ...'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-in_file', help='mtb or top file', type=str, required=True)
    parser.add_argument('-ifp', help='ifp file', type=str, default=None)
    parser.add_argument('-flag_top', help='if defined, in_file will be parsed as a top', default=False, action='store_true')
    parser.add_argument('-bb', help='building block ID', type=str, required=True)
    parser.add_argument('-pdb', help='pdb file', type=str, required=True)
    parser.add_argument('-sel', help='selection (using pymol, e.g. resn FMN)', type=str, default='all')
    
    args = parser.parse_args()

    assert os.path.isfile(args.in_file)
    assert os.path.isfile(args.pdb)
    if args.flag_top is False:
        assert os.path.isfile(args.ifp)

    in_file = args.in_file
    ifp_file = args.ifp
    bb_id = args.bb
    pdb_file = args.pdb
    selection = args.sel
    flag_topology = args.flag_top

    # get the matches
    print('\n\n\tRUNNING...\n')
    mcs = find_matching(in_file, ifp_file, bb_id, pdb_file, selection=selection, flag_topology=flag_topology)

    print('\n\nsolution number\t\tnumber of matched atoms')
    results = []
    for i in range(len(mcs.solutions)):
        res = check_elements(mcs, i)
        results.append((i, res[0][1]))
    results.sort(key=lambda x:x[1], reverse=True)

    for i, res in results:
        print(i, '\t\t\t\t', res)

    print('\n\nsolution with the most matches is: {}; and has {:} matched atom'.format(*results[0]))
    print('run the following commands to get the matches:\nsol_index = {:}'.format(results[0][0]))
    print('res = check_elements(mcs, sol_index)\nres[2]')
    print('the actual match is in res[2]._sol (pandas DataFrame)')
    print('\nto check other solutions, change sol_index')
    print('variable results is a sorted list of solutions by the number of matched atoms')

    sol_index = results[0][0]
    res = check_elements(mcs, sol_index)
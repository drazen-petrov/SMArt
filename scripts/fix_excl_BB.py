from SMArt.md import parse_mtb

if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--mtb', help='mtb file', type=str, required=True)
    parser.add_argument('-b', '--bb', help='building block to add exclusions', type=str, required=True)
    parser.add_argument('-n', '--nexcl', help='number of excluded atoms (based on bonds)', type=int)
    parser.add_argument('-a', '--atom_subset', help='subset of atoms to add exclusions to (based on atom names)', type=str, nargs='*')
    parser.add_argument('--flag_1_neigh_subset', default=False, action='store_true', help='add first neighbours of the atom_subset to the subset')
    parser.add_argument('--flag_rings', default=False, action='store_true', help='use rings to generate atom_subset')
    parser.add_argument('-o', '--out_mtb', type=str, help='output file to write out fixed mtb file')
    parser.add_argument('--print_rings', default=False, action='store_true', help='finds the rings in the structure and prints it out (and exit)')
    
    args = parser.parse_args()
    fpath = os.path.abspath(args.mtb)

    mtb = parse_mtb(fpath)
    bb = mtb.bb[args.bb]
    bb.gen_graph()

    flag_run = True
    if args.print_rings:
        rings = bb.find_rings()
        print('rings:')
        for r in rings:
            print('\t',r)
        flag_run = False
    
    if flag_run:
        assert args.nexcl # provide the number of exclusions if you want to add them
        # prepare atom_subset
        if args.atom_subset and args.flag_rings:assert False # provide only one of these
        atom_subset = set()
        if args.atom_subset:
            atom_subset = set(bb.find_atoms(*args.atom_subset))
        if args.flag_rings:
            atom_subset = set()
            for r in bb.find_rings():
                atom_subset |= set(r)
        if args.flag_1_neigh_subset:
            assert atom_subset # using `flag_1_neigh_subset` makes only sense if atom_subset is defined
            for at in list(atom_subset):
                atom_subset |= bb.adj[at]
        # add exclusions
        bb.get_EP_l()
        bb.add_exclusions_neigh(args.nexcl, atom_subset=atom_subset)
        bb.transform_EP_l()
        if args.out_mtb:
            fout_path = args.out_mtb
            if not fout_path.endswith(".mtb"):
                fout_path += ".mtb"
            print('writing out to', fout_path)
            mtb.write_mtb(fout_path)
        else:
            print("WARN: no output file defined")

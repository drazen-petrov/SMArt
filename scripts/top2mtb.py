from SMArt.md import parse_mtb, parse_top, parse_ff
from SMArt.incl import GeneralContainer
from SMArt.md.gromos.io.incl import GromosFile
from SMArt.md.data_st import BuildingBlock, BBAtom, Interaction

if __name__ == '__main__':
    import argparse
    import os

    excl_txt = "add exclusions between BB atoms and preceding/trailing atoms based on the template building block"
    excl_bond = "add bonded interactions between BB atoms and preceding/trailing atoms based on the template building block"
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--top', help='top file', type=str, required=True)
    parser.add_argument('--bb_name', help='building block name', type=str)
    parser.add_argument('-a', '--atom_subset', help='subset of atoms (reduced topology)', type=str, nargs='*')
    #parser.add_argument('--no_reorder', default=False, action='store_true', help='do not change the order of atoms')
    parser.add_argument('-o', '--out_mtb', type=str, help='output file to write out fixed mtb file')
    parser.add_argument('--print_atoms', default=False, action='store_true', help='prints all atoms and exits')
    parser.add_argument('--mtb', help='mtb file', type=str)
    parser.add_argument('-b', '--bb', help='template building block for preceding/trailing atoms', type=str)
    parser.add_argument('--ifp', help='ifp file', type=str)
    parser.add_argument('--flag_excl', default=False, action='store_true', help=excl_txt)
    parser.add_argument('--flag_bond', default=False, action='store_true', help=excl_bond)
    parser.add_argument('--m_diff', default=0.01, help="threshold for printing a warning (difference to the closest m_type)")

    args = parser.parse_args()
    top_path = os.path.abspath(args.top)
    
    top = parse_top(top_path)
    if args.print_atoms:
        s_out = ""
        for at in top.get_atoms():
            print(at.gr_id, at.name)
            s_out += at.name + " "
        print(s_out)
        exit()
    else:
        assert args.ifp # needed for mass types

    # mtb file
    bb = None
    if args.bb:
        assert args.mtb # if bb provided, mtb is also needed
    if args.mtb:
        mtb = parse_mtb(args.mtb, ifp_file=args.ifp)
        assert args.bb
        bb = mtb.bb[args.bb]
    # force field (for mass types)
    if args.ifp:
        ff = parse_ff(args.ifp)
    # prep the topology (reduce if needed)
    rt = top
    residues = list(top.residues.values())
    if args.atom_subset:
        assert len(residues)==1
        top_atom_subset = list(residues[0].find_atoms(*args.atom_subset))
        top_atom_id_subset = [at.id for at in top_atom_subset]
        rt = top.reduce_top(top_atom_id_subset)
        # reoder the atoms
        temp_atms = GeneralContainer()
        for at_id in top_atom_id_subset:
            temp_atms.add2container(rt.atoms[at_id], create=True)
        temp_atms.renumber_container('atoms')
        rt.sort_container('atoms')
    
    # make a new bb
    if args.bb_name:
        bb_name = args.bb_name
    else:
        bb_name = residues[0].name
    new_bb = BuildingBlock(bb_name)
    new_bb.n_at = len(rt.atoms)
    
    ### atoms ###
    # atoms from the template BB
    post_at = []
    pre_post_bb_atoms = []
    if bb:
        new_bb.n_pre_at = bb.n_pre_at
        for at in bb.get_atoms():
            if int(at.id)<1:
                new_bb.add2container(at, create=True)
            elif at.flag_bb is False:
                post_at.append(at)
            if at.flag_bb:
                if int(at.id) <= bb.n_pre_at or at.flag_excl is False:
                    pre_post_bb_atoms.append(at)
    else:
        new_bb.n_pre_at = 0
    # atoms from the topology
    pre_post_new_bb_atoms = []
    for i,at in enumerate(rt.get_atoms()):
        if i < new_bb.n_pre_at:
            pre_post_new_bb_atoms.append(at)
        if i < new_bb.n_at-new_bb.n_pre_at:
            at.flag_excl = True
        else:
            at.flag_excl = False
            pre_post_new_bb_atoms.append(at)
        at.flag_bb = True
        new_bb.add2container(at, create=True)
    last_id = int(at.id)
    for at in post_at:
        last_id +=1
        at.id = str(last_id)
        new_bb.add2container(at)
    for at in new_bb.get_atoms():
        at.gr_id = at.id
    # map between template bb and new bb (for preceding/trailing atoms)
    pre_post_bb_map = {}
    for i, at in enumerate(pre_post_bb_atoms):
        pre_post_bb_map[at.id] = pre_post_new_bb_atoms[i].id
    # function to get new bb atom from the template bb (for preceding/trailing atoms)
    def get_new_bb_at(at):
        if at in pre_post_bb_atoms:
            new_bb_at = new_bb.atoms[pre_post_bb_map[at.id]]
        else:
            new_bb_at = new_bb.atoms[at.id]
        return new_bb_at
    ### excl ###
    e_l, _ = rt._generate_atom_excl_pair_list()
    for at in new_bb.get_atoms():
        if at.flag_bb:
            if at.flag_excl:
                at.e_l = list(e_l[at])
        else:
            new_e_l = []
            for temp_at in at.e_l:
                if temp_at.flag_bb:
                    # `at` is either pre or post atom
                    # if temp_at is pre/post atom or one of the first/last n_pre_at -> add exclusion
                    #if not temp_at.flag_bb or int(temp_at.id) <= new_bb.n_pre_at or not temp_at.flag_excl:
                    if not temp_at.flag_bb or new_bb.atoms[temp_at.id] in pre_post_new_bb_atoms:
                        new_e_l.append(temp_at)
                    # if temp_at is one of the building block atoms that are not among first/last n_pre_at 
                    # decision based on flag_excl
                    else:
                        print("exclusion between:",at, temp_at,
                              "\t\t\tin new bb:", new_bb.atoms[temp_at.id], "\t\tflag_excl:", args.flag_excl)
                        if args.flag_excl:
                            new_e_l.append(new_bb.atoms[temp_at.id])
            at.e_l = new_e_l
    if bb:
        for at in bb.get_atoms():
            if at.flag_bb:
                new_bb_at = get_new_bb_at(at)
                for temp_at in list(at.e_l):
                    if not temp_at.flag_bb:
                        print("exclusion between:",at, temp_at, "\t\tflag_excl:", args.flag_excl)
                        if int(at.id) > new_bb.n_pre_at:
                            if args.flag_excl:
                                new_bb_at.e_l.append(temp_at)
                        else:
                            new_bb_at.e_l.append(temp_at)
    ### bonded interactions ###
    for int_cont in rt.get_interaction_containers():
        for b in int_cont:
            new_bb.add2container(b, create=True, db_type=list)
    if bb:
        for int_cont in bb.get_interaction_containers():
            for b in int_cont:
                flag_non_bb_atms = False
                flag_non_pre_post_atm = False
                for at in b.atoms:
                    if at.flag_bb:
                        new_bb_at = get_new_bb_at(at)
                    else:
                        flag_non_bb_atms = True
                if flag_non_bb_atms:
                    new_bb_atoms = []
                    for at in b.atoms:
                        if at.flag_bb:
                            new_bb_at = get_new_bb_at(at)
                        else:
                            new_bb_at = at
                        new_bb_atoms.append(new_bb_at)
                    new_b = Interaction(b.int_type, new_bb_atoms, b.states)
                    print(b.int_type, b.atoms, "\t\t\tin new bb:", new_bb_atoms,"\t\tflag_bond:", args.flag_bond)
                    if flag_non_pre_post_atm:
                        if args.flag_bond:
                            new_bb.add2container(b, create=True, db_type=list)
                    else:
                        new_bb.add2container(b, create=True, db_type=list)
    ### fix bb atoms ###
    new_bb_atoms = list(new_bb.get_atoms())
    for at in new_bb_atoms:
        if at.flag_bb:
            new_at = BBAtom(at.id)
            new_at.gr_id = at.id
            new_at.flag_bb = at.flag_bb
            new_at.flag_excl = at.flag_excl
            new_at.name = at.name
            new_at.p_ch = at.p_ch
            new_at.a_type = at.a_type
            new_at.mk_cg = at.mk_cg
            new_at.e_l = at.e_l
            # find the best fitting m_type
            min_m_diff = 10000
            min_m_type = None
            for m_type in mtb.ff.m_type.values():
                if abs(at.m - m_type.m) < min_m_diff:
                    min_m_type = m_type
                    min_m_diff = abs(at.m - m_type.m)
            if min_m_diff > args.m_diff:
                print("WARN: no identical mass type found for", at, at.m, "\t\tclosest:", min_m_type, "\t\tdiff:", min_m_diff)
            new_at.m_type = min_m_type
            new_bb.atoms[new_at.id] = new_at
    ### fix lj_exceptions ###
    new_bb.lj_exceptions = ['    0\n']

    ### write out ###
    assert args.out_mtb
    f_new = GromosFile(args.out_mtb, 'w')
    new_bb.write_bb(f_new)
    f_new.f.close()

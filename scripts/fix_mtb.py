import SMArt
from SMArt.md.data_st import BBdb
from SMArt.graph import Graph, GraphDirected
from SMArt.md.gromos.io.incl import GromosFile

def fix_mtb_G_names(bb, flag_add_neg_atoms=False):
    m = 1
    if flag_add_neg_atoms:
        for at in bb.adj:
            if int(at.id) < m:
                m = int(at.id)
        for i in range(m, 1):
            at1 = str(i)
            if at1 in bb.atoms:
                at1 = bb.atoms[at1]
            at2 = str(i + 1)
            if at2 in bb.atoms:
                at2 = bb.atoms[at2]
            bb.add_edge((at1, at2))
    else:
        at_m = []
        for at in bb.adj:
            if int(at.id) < 1:
                at_m.append(at)
            if int(at.id) == 1:
                AT_1 = at
        for at_i in range(len(at_m) - 1):
            at1 = at_m[at_i]
            at2 = at_m[at_i + 1]
            bb.add_edge((at1, at2))
        at1 = at_m[-1]
        bb.add_edge((at1, AT_1))
    for at in bb.atoms:
        if not hasattr(bb.atoms[at], 'name'):
            bb.atoms[at].name = at


if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='f', help='mtb file', type=str, required=True)
    parser.add_argument('-bb', help='building blocks to fix', type=str, nargs='*', default=None)
    parser.add_argument('-add_neg_atoms',
                        help='add atoms with negative IDs, e.g. if only -2 is there and one wants to have -1',
                        default=False, action='store_true')
    parser.add_argument('-plot', dest='plot', help='plot structure', default=False, action='store_true')
    parser.add_argument('-step', dest='step', help='plot structure at each step', default=False, action='store_true')
    parser.add_argument('-init', dest='init', help='plot initial structure at each step', default=False,
                        action='store_true')
    parser.add_argument('-mov', help='make a movie from structures at each step', default=False, action='store_true')
    parser.add_argument('-at_names', help='use atom names, not atom ids', default=False, action='store_true')
    parser.add_argument('-v', help='verbose', default=False, action='store_true')

    args = parser.parse_args()
    fpath = os.path.abspath(args.f)

    mtb = SMArt.md.parse_mtb(fpath)
    if args.bb is None:
        if len(mtb.bb) > 1:
            raise Exception('mtb file has more than 1 BB, use -bb flag')
        bb_names = [list(mtb.bb.keys())[0]]
    else:
        if args.bb[0] == 'all':
            bb_names = list(mtb.bb.keys())
        else:
            bb_names = list(args.bb)
    for bb_name in bb_names:
        if bb_name not in mtb.bb:
            SMArt.incl.do_warn('BB {:} not in the mtb file'.format(bb_name))
            continue
        print(bb_name)
        if args.bb is None:
            suf_fname = ''
        else:
            suf_fname = '_' + bb_name
        mtb.bb[bb_name].create_adj_list()
        fix_mtb_G_names(mtb.bb[bb_name], args.add_neg_atoms)
        new_G = GraphDirected(mtb.bb[bb_name].adj)
        fd, fname = os.path.split(fpath)

        fig_name = None
        if args.plot:
            fig_name = 'structure_' + fname[:-4] + suf_fname
        flag_names = 'id'
        if args.at_names:
            flag_names = 'name'

        kwargs = dict(fig_name=fig_name, step=args.step, init_struc=args.init, flag_name_attribute=flag_names)
        kwargs['verbose'] = args.v

        mc = new_G.get_2D_repr(**kwargs)

        if args.mov and args.plot:
            comm = 'convert -delay  100 ' + fig_name + '* ' + fig_name + '_movie.gif'
            os.system(comm)

        f_new = GromosFile('fixed_' + fname[:-4] + suf_fname + '.mtb', 'w')
        pos_txt = 'XYPOSITION\n' + str(len(mtb.bb[bb_name].atoms)) + '\n'
        for at in mtb.bb[bb_name].get_atoms():
            pos_txt += '{:<6} {:10.2f} {:10.2f}\n'.format(at.id, mc.loc[at,'x'] * 35, mc.loc[at,'y'] * 35)

        if args.bb is None:
            f = open(fpath)
            for i in f:
                f_new.f.write(i)
            f.close()
            f_new.f.write(pos_txt)
            f_new.f.write('END')
        else:
            f_new.write_block('MTBUILDBLSOLUTE', mtb.bb[bb_name].write_bb(get_str=True) + pos_txt)
        f_new.f.close()




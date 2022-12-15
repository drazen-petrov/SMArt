import numpy as np
from SMArt.md import parse_trc 
import mdtraj as md

def convert_trj_gr2gm(trc, gro, out):
    """
    convert gromos trajectory (trc) into gromacs (e.g. xtc format) - it only works with cubic box!

    :param trc: gromos trajectory file (*.trc)
    :param gro: gro or pdb file that corresponds to the trajectory
    :param out: output file (*xtc or *trr)

    :return: None
    """
    # read the trajectory (*trc)
    trc = parse_trc(trc)

    # create mdtraj object
    temp_md_trj = md.load(gro)
    trj_md = md.core.trajectory.Trajectory(trc.trj['coord'], temp_md_trj.top)

    # add the box (only works with cubic box)
    trj_md.unitcell_angles = trc.trj['box'][:,1]
    trj_md.unitcell_lengths = trc.trj['box'][:,0]
    unitcell_vectors = []
    for frame_box in trc.trj['box']:
	    unitcell_vectors.append(np.diag(frame_box[0]))
    trj_md.unitcell_vectors = np.array(unitcell_vectors)

    # write out
    if args.out.endswith('xtc'):
    	trj_md.save_xtc(args.out)
    elif args.out.endswith('trr'):
	    trj_md.save_trr(args.out)

if __name__ == "__main__":
    # parsing arguments
    #import argparse
    #parser = argparse.ArgumentParser()
    from SMArt.incl import ArgParser # same as argparse with a fromfile_prefix_chars fix
    parser = ArgParser(fromfile_prefix_chars='@')
    parser.add_argument('-trc', type=str, required = True, help='trc file')
    parser.add_argument('-gro', type=str, required = True, help='gro file - first frame or equilibrated structure (could be a pdb file as well)')
    parser.add_argument('-out', type=str, default='out.xtc', help='output file name (*xtc or *trr)')

    args = parser.parse_args()

    # convert
    convert_trj_gr2gm(args.trc, args.gro, args.out)


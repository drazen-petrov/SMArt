from SMArt.md import Trajectory
import mdtraj as md

def convert_trj_gr2gm(gmx_trj_path, gro, out_f, file_format='trc'):
    """
    convert gromacs trajectory (e.g. xtc) into gromos (e.g. trc format) - it was tested with cubic boxes only!

    :param gmx_trj_path: gromacs trajectory file (e.g. *.xtc)
    :param gro: gro or pdb file that corresponds to the trajectory (used to load the trajectory with mdtraj)
    :param out_f: output file (e.g. *trc)
    :param file_format: trc or npz (npz is a numpy compressed version of the trajectory)

    :return: None
    """
    # read the trajectory (*trc)
    gmx_trj = md.load(gmx_trj_path, top=gro)

    # create Trajectory object and copy the data
    N_frm, N_atm = gmx_trj.xyz.shape[:2]
    gro_trj = Trajectory(N_atoms=N_atm)
    gro_trj.add_gr_title(['from files: {:}, {:}'.format(gmx_trj_path, gro)])
    gro_trj.set_N_frames(N_frm)
    gro_trj.trj['time_step'][:,0] = list(range(N_frm))
    gro_trj.trj['time'][:,0] = gmx_trj.time.copy()
    gro_trj.trj['coord'][:] = gmx_trj.xyz.copy()
    print('WARN: tested for cubic boxes!')
    gro_trj.trj['box'][:,0,:] = gmx_trj.unitcell_lengths.copy()
    gro_trj.trj['box_type'] = 1
    gro_trj.trj['box'][:,1,:] = gmx_trj.unitcell_angles.copy()
    gro_trj.trj['box'][:,2:,:] = 0

    # write out
    assert file_format in ('trc', 'npz')
    if file_format=='trc':
        assert out_f.endswith('trc')
        gro_trj.write_trc(out_f)
    elif file_format=='npz':
	    gro_trj.write_trc_npz(out_f)

if __name__ == "__main__":
    # parsing arguments
    #import argparse
    #parser = argparse.ArgumentParser()
    from SMArt.incl import ArgParser # same as argparse with a fromfile_prefix_chars fix
    parser = ArgParser(fromfile_prefix_chars='@')
    parser.add_argument('-gmx_trj', type=str, required = True, help='gromacs trajectory file (e.g. *.xtc)')
    parser.add_argument('-gro', type=str, required = True, help='gro or pdb file that corresponds to the trajectory (used to load the trajectory with mdtraj)')
    parser.add_argument('-out', type=str, default='out.trc', help='output file (e.g. *trc)')
    parser.add_argument('-file_format', type=str, default='trc', help='trc or npz')

    args = parser.parse_args()

    # convert
    convert_trj_gr2gm(args.gmx_trj, args.gro, args.out, args.file_format)


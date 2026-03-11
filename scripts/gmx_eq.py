import os
import shutil
import numpy as np
from SMArt.md.wrappers import run_gm_prog
from SMArt.md.data_st import MD_Parameters
from SMArt.md import parse_top

def get_num_atoms(template_gro):
    with open(template_gro) as f:
        f.readline()  # skip first line
        n_atoms = int(f.readline().strip())
    return n_atoms

def generate_posre(atoms, fk, posre_out, exclude_atoms=None):
    if exclude_atoms is None:
        exclude_atoms = set()
    f=open(posre_out, "w")
    f.write("[ position_restraints ]\n")
    for at in atoms:
        if at not in exclude_atoms:
            temp_l=str(at).rjust(6)+"     1"
            for _ in range(3):
                temp_l += "{:>10.1e}".format(fk)
            f.write(temp_l+"\n")
    f.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--fd', help='root folder to generate files')
    parser.add_argument('-c', '--conf', required = True)
    parser.add_argument('-t', '--top', required = True)
    parser.add_argument('-m', '--mdp', required = True)
    parser.add_argument('-ndx', help='index file')
    parser.add_argument('-t_eq', default=0.1, type=float, help='equilibration time per simulation in ns')
    parser.add_argument('-ref_conf', help='reference configuration file for position restraints - if not given, use conf')
    parser.add_argument('-posres_first_N', type=int, help='position restrain first N atoms from the reference configuration')
    parser.add_argument('-posres_first_N_gro_file', type=str, help='position restrain first N atoms based on a gro file (just total number of atoms)')
    parser.add_argument('-posres_from_top', type=str, help='position restrain deretmined from a topology (first molecule)')
    parser.add_argument('-posres_heavy_only', default=False, action='store_true', help='position restrain heavy atoms only')
    parser.add_argument('-posres_filename', default='posre.itp', help='poseres filename that is included in the topology')
    parser.add_argument('-posres_k', default = np.power(10.,np.arange(4,-2,-1)) * 2.5, help='position restraint force constants for eq steps', type = float, nargs = '+')
    parser.add_argument('-T_eq', default = np.arange(50,301,50), type = float, nargs = '+', help='temperature for equilibration steps')
    parser.add_argument('-no_npt_step', default=False, action='store_true', help='skip NPT equilibration step at the end')
    parser.add_argument('-n', '--nsim', type = int, help = 'number of ind replicas', default = 3)
    parser.add_argument('-ncpu', dest='ncpu', type = int, help = 'number of cpu', default = 4)
    parser.add_argument('-remove_eq_fd', default=False, action='store_true', help='remove existing equilibration folders')
    args = parser.parse_args()

    # get absolute paths
    args.fd = os.path.abspath(args.fd)
    args.conf = os.path.abspath(args.conf)
    args.top =  os.path.abspath(args.top)
    args.mdp = os.path.abspath(args.mdp)
    if args.ref_conf:
        args.ref_conf = os.path.abspath(args.ref_conf)
    if args.ndx:
        args.ndx = os.path.abspath(args.ndx)

    atoms2posres = None
    n_posres_atoms = None
    if args.posres_first_N_gro_file is not None:
        n_posres_atoms = get_num_atoms(args.posres_first_N_gro_file)
        atoms2posres = list(range(1, n_posres_atoms + 1))
    if args.posres_first_N is not None:
        assert atoms2posres is None
        atoms2posres = list(range(1, args.posres_first_N + 1))
    
    if args.posres_from_top:
        top = parse_top(args.posres_from_top, format_type='gm')
        
        if args.posres_heavy_only:
            gmx_atoms = top.molecules[0].mol_type.get_HH()[1]
        else:
            gmx_atoms = top.molecules[0].mol_type.get_atoms()
        atoms2posres = [at.id for at in gmx_atoms]
    
    assert len(args.posres_k) == len(args.T_eq), "number of position restraint force constants must be equal to number of temperatures"

    mdp = MD_Parameters(args.mdp)
    top_fname = os.path.basename(args.top)
    
    for sim_i in range(args.nsim):
        target_dir = os.path.join(args.fd, 'eq_'+str(sim_i))
        if args.remove_eq_fd:
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir, exist_ok=True)
        os.system('cp {} {}/'.format(args.top, target_dir))
        run_file_path = os.path.join(target_dir, 'run_eq_{}.sh'.format(sim_i))
        with open(run_file_path, 'w') as f_sh:
            f_sh.write('cd {}\n'.format(target_dir))
            f_sh.write('\n')
            # generate and run equilibration simulations for each step
            current_conf = args.conf
            for eq_sim_i in range(len(args.posres_k)):
                posres_fk = args.posres_k[eq_sim_i]
                temp_T = args.T_eq[eq_sim_i]
                # generate posre file
                posre_out = None
                if posres_fk>1:
                    posre_out = os.path.join(target_dir, 'posre_{}.itp'.format(eq_sim_i))
                    generate_posre(atoms2posres, posres_fk, posre_out)

                # generate mdp with temp_T etc
                mdp_out = os.path.join(target_dir, 'mdp_{}.mdp'.format(eq_sim_i))
                mdp_kw = dict(nsteps = int(1000 * args.t_eq / 0.002))  # 2 fs time step
                mdp_kw['ref_t'] = f"{temp_T} {temp_T}"
                mdp_kw['gen_temp'] = temp_T
                mdp_kw['gen_vel'] = "yes"
                mdp_kw['pcoupl'] = 'no'  # turn off pressure coupling during equilibration
                if posre_out:
                    mdp_kw['define'] = '-DPOSRES'
                mdp.change_mdp(mdp_kw, mdp_out = mdp_out)

                # generate tpr
                run_name = 'eq_{}_{}'.format(sim_i, eq_sim_i)
                additional_gmx_kwargs = {}
                if posre_out:
                    f_sh.write('cp {} {}\n'.format(posre_out, args.posres_filename))
                    if args.ref_conf:
                        additional_gmx_kwargs['r'] = args.ref_conf
                    if args.ndx:
                        additional_gmx_kwargs['n'] = args.ndx
                comm = run_gm_prog('grompp', c=current_conf, p=top_fname, f=mdp_out,
                                o=run_name, maxwarn=3, **additional_gmx_kwargs)
                f_sh.write(comm + '\n')

                # run md
                comm = run_gm_prog('mdrun', deffnm=run_name, nt=args.ncpu)
                f_sh.write(comm + '\n')
                f_sh.write('\n')
                current_conf = run_name + '.gro'
            
            if not args.no_npt_step:
                # final NPT equilibration step without position restraints
                f_sh.write('\n# final NPT equilibration step without position restraints\n')
                mdp_out = os.path.join(target_dir, 'mdp_{}.mdp'.format('npt'))
                mdp_kw = dict(nsteps = int(1000 * args.t_eq / 0.002))  # 2 fs time step
                mdp_kw['ref_t'] = f"{args.T_eq[-1]} {args.T_eq[-1]}"
                mdp_kw['gen_vel'] = "no"
                mdp_kw['pcoupl'] = 'C-rescale'  # turn on pressure coupling during NPT equilibration
                mdp.change_mdp(mdp_kw, mdp_out = mdp_out)

                # generate tpr
                run_name = 'eq_{}_npt'.format(sim_i)
                additional_gmx_kwargs = {}
                if args.ref_conf:
                    additional_gmx_kwargs['r'] = args.ref_conf
                if args.ndx:
                    additional_gmx_kwargs['n'] = args.ndx
                comm = run_gm_prog('grompp', c=current_conf, p=top_fname, f=mdp_out,
                                o=run_name, maxwarn=3, **additional_gmx_kwargs)
                f_sh.write(comm + '\n')

                # run md
                comm = run_gm_prog('mdrun', deffnm=run_name, nt=args.ncpu)
                f_sh.write(comm + '\n')
                f_sh.write('\n')
        os.system("chmod u+x " + run_file_path)

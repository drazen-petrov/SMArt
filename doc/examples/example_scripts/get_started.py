import sys
import os

smart_fd = os.path.abspath(sys.argv[0])
smart_fd = os.path.split(smart_fd)[0]
out_data_fd = os.path.join(smart_fd, 'out_data')

for _ in range(3):
    smart_fd = os.path.split(smart_fd)[0]

sys.path.insert(0, smart_fd)

import SMArt
from SMArt.md import parse_top, parse_ff
from SMArt.md import FF, BBdb, Topology, Configuration, MolSystem

from SMArt.md.gromos.io.incl import GromosFile, GromosString
from SMArt.md.gromacs.io.incl import GromacsFile, GromacsString

SMArt.incl.set_print_warnings(False)

if not os.path.isdir(out_data_fd):
    os.mkdir(out_data_fd)
os.chdir(out_data_fd)

### gromos
# parsing from string

ff = FF()
ff.parse_gro('TITLE\ntest title\nEND', parse_from_file = False)
print('ff_title')
print(ff.sys_title[0].lines) # title lines are found here
print('')

masstypes = """MASSATOMTYPECODE
#   NRMATY    NMATY
       21     80
#       Table 2.4.1
#       GROMOS mass atom type codes, masses and names
#
# mass atom type code
# mass
# mass atom name
#
#
# N     ATMAS   ATMASN
1       1.008   H
3       13.019  CH1
4       14.027  CH2
5       15.035  CH3
6       16.043  CH4
12      12.011  C
14      14.0067 N
16      15.9994 O
19      18.9984 F
23      22.9898 NA
24      24.305  MG
28      28.08   SI
31      30.9738 P
32      32.06   S
35      35.453  CL
39      39.948  AR
40      40.08   CA
56      55.847  FE
63      63.546  CU
65      65.37   ZN
80      79.904  BR
END
"""
gs = GromosString(masstypes)
ff.parse_gro(gs)
print('ff mass types')
print(ff.m_type) # mass type are stored in ff.m_type

ff.write_ff('mass_types_only.ifp')

# parsing force field (FF) from file
fd_gromos = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromos')

ifp_file = os.path.join(fd_gromos, '54a8.ifp')
ff = parse_ff(ifp_file)
ff.write_ifp('test.ifp') # write gromos FF ifp file

# convert to gromacs format
ff.gr2gm()
ff.write_ff_itp('gr2gm_ff.itp') # write gromacs itp file

# dump using pickle or yaml

ff.dump('ff_pickle_dump.p')
ff.dump('ff_yaml_dump.yaml', 'yaml')

ffp = SMArt.incl.DataDumping.load('ff_pickle_dump.p')
ffy = SMArt.incl.DataDumping.load('ff_yaml_dump.yaml', 'yaml')

# write again (check with diff to test.ifp)
ffp.write_ifp('testp.ifp')
ffy.write_ifp('testy.ifp')

# GROMOS mtb files
mtb_file = os.path.join(fd_gromos, '54a8.mtb')
mtb = BBdb(mtb_file) # read an mtb file - same as md.parse_mtb(fd_gromos + '54a8.mtb')
mtb.write_mtb('mtb1.mtb') # write out an mtb file

# dumping / writing again
mtb.dump('mtb_pickle_dump.p')
mtbp = SMArt.incl.DataDumping.load('mtb_pickle_dump.p')
mtbp.write_mtb('testp.mtb')


# mtb file together with the ifp file
mtb2 = BBdb(mtb_file, ifp_file)
mtb2.write_mtb('mtb2.mtb')

# dumping / writing again
mtb2.dump('mtb2_pickle_dump.p')
mtb2p = SMArt.incl.DataDumping.load('mtb2_pickle_dump.p')
mtb2p.write_mtb('test2p.mtb')

# GROMOS topology
top_file = os.path.join(fd_gromos, 'test_pept.top')
top = parse_top(top_file) # parsing
top.write_top('test.top') # writing


# dumping / writing again (not working with python2 and some platforms
try:
    top.dump('pickle_dump.p')
    topp = SMArt.incl.DataDumping.load('pickle_dump.p')
    topp.write_top('testp.top')
except:pass

top2_file = os.path.join(fd_gromos, 'top.top')
top = parse_top(top2_file, format_type='gr')
SMArt.incl.set_print_warnings()

top.write_top('lha_top.top')


# convert to gromacs format
top.gr2gm()

out_gr2gm_data_fd = os.path.join(out_data_fd, 'gr2gm_top')
if not os.path.isdir(out_gr2gm_data_fd):
    os.mkdir(out_gr2gm_data_fd)
os.chdir(out_gr2gm_data_fd)
top.write_top_gm('gromos2gomacs_top.top', from_str = 1, flag_segment_order=1, flag_if=1, flag_include = 1, sep_ff2itp='top_ff.itp', flag_defines_first_include = 1, sep_mol2itp=1)
# this produces separate itp files for each molecule - including separate files for each ion (mol_41 mol_42 ... )

# to avoid this, ions could be assigned to predefined ion molecule types using this code
for m in top.molecules:
    mt2del = set()
    if len(m.mol_type.atoms) == 1:
        for temp_mt_id in ('CU1', 'CU', 'ZN', 'MG', 'CA', 'NA', 'CL'):
            if temp_mt_id in top.molecule_types:
                temp_mt = top.molecule_types[temp_mt_id]
                at, at2 = list(m.mol_type.atoms.values())[0], temp_mt.atoms[1]
                if at.p_ch == at2.p_ch and at.m == at2.m and at.a_type == at2.a_type:
                    mt2del.add(m.mol_type)
                    m.mol_type = temp_mt
    for temp_mt in mt2del:
        top.molecule_types.pop(temp_mt.id)

out_gr2gm_data_fd = os.path.join(out_data_fd, 'gr2gm_top_ions')
if not os.path.isdir(out_gr2gm_data_fd):
    os.mkdir(out_gr2gm_data_fd)
os.chdir(out_gr2gm_data_fd)

top.write_top('gromos2gomacs_top2.top', from_str = 1, flag_segment_order=1, flag_if=1, flag_include = 1, sep_ff2itp='top_ff.itp', flag_defines_first_include = 1, sep_mol2itp=1, format_type='gm')


# some extra tricks with the topology object (getting molecules, adding a bond)
print('')
mols = top.get_molecules(flag_sort = True) # without flag_sort, order of atoms determined by the DFS algorithm
flag_check, last_at_index, mols2 = top.check_molecule_atoms_order(molecules = mols, flag_mol_atom_sort = 1)
print('molecule sorting\n', flag_check, last_at_index)

# adding a bond (cross link - between molecules)
print('')
b = top.Interaction(top.BondType)
b.add_atom(top.atoms['1'], top.atoms['400']) # atoms 1 and 400
b.add_state(top.ff.gr_bonds['1']) # bond type 1
top.add2container(b) # general way of adding information (atoms, interactions, etc.) to the topology object
top.write_top('lha_top2.top')

print('')
mols3 = top.get_molecules(flag_sort = 1)
flag_check, last_at_index, mols4 = top.check_molecule_atoms_order(molecules = mols3, flag_mol_atom_sort = 1)
print('')
print('molecule sorting after a cross link bond\n', flag_check, last_at_index)
## look at molecules in mols, mols2, mols3, mols4


print('\n')
print('data of any object is stored in _containers')
print('\ntop containers\n', top._containers)
print('\ntop.ff containers\n', top.ff._containers)

#SMArt.incl.set_print_warnings(False)
# gromacs
# parsing from string
print('\n\ngromacs part')
temp ="""[ atomtypes ]
;name  at.num   mass      charge  ptype       c6           c12
    O    8      0.000      0.000     A  0.0022619536       1e-06
   OM    8      0.000      0.000     A  0.0022619536  7.4149321e-07
   OA    8      0.000      0.000     A  0.0022619536  1.505529e-06
   OE    8      0.000      0.000     A  0.0022619536    1.21e-06
   OW    8      0.000      0.000     A  0.0026173456  2.634129e-06
    N    7      0.000      0.000     A  0.0024364096  2.319529e-06
   NT    7      0.000      0.000     A  0.0024364096  5.0625e-06
   NL    7      0.000      0.000     A  0.0024364096  2.319529e-06
   NR    7      0.000      0.000     A  0.0024364096  3.389281e-06
   NZ    7      0.000      0.000     A  0.0024364096  2.319529e-06
   NE    7      0.000      0.000     A  0.0024364096  2.319529e-06
    C    6      0.000      0.000     A  0.0023406244  4.937284e-06
  CH0    6      0.000      0.000     A  0.0023970816  0.0002053489
  CH1    6      0.000      0.000     A  0.00606841  9.70225e-05
  CH2    6      0.000      0.000     A  0.0074684164  3.3965584e-05
  CH3    6      0.000      0.000     A  0.0096138025  2.6646244e-05
  CH4    6      0.000      0.000     A  0.01317904  3.4363044e-05
 CH2r    6      0.000      0.000     A  0.0073342096  2.8058209e-05
  CR1    6      0.000      0.000     A  0.0055130625  1.5116544e-05
   HC    1      0.000      0.000     A   8.464e-05  1.5129e-08
    H    1      0.000      0.000     A           0           0
  DUM    0      0.000      0.000     A           0           0
    S   16      0.000      0.000     A  0.0099840064  1.3075456e-05
 CU1+   29      0.000      0.000     A  0.0004182025  5.1251281e-09
 CU2+   29      0.000      0.000     A  0.0004182025  5.1251281e-09
   FE   26      0.000      0.000     A           0           0
 ZN2+   30      0.000      0.000     A  0.0004182025  9.4400656e-09
 MG2+   12      0.000      0.000     A  6.52864e-05  3.4082244e-09
 CA2+   20      0.000      0.000     A  0.00100489  4.9801249e-07
    P   15      0.000      0.000     A  0.01473796  2.2193521e-05
   AR   18      0.000      0.000     A  0.0062647225  9.847044e-06
    F    9      0.000      0.000     A  0.0011778624  7.6073284e-07
   CL   17      0.000      0.000     A  0.0087647044  1.5295921e-05
   BR   35      0.000      0.000     A  0.02765569  6.5480464e-05
 CMet    6      0.000      0.000     A  0.0088755241   1.936e-05
 OMet    8      0.000      0.000     A  0.0022619536  2.325625e-06
  NA+   11      0.000      0.000     A  7.884019264e-05 7.290000e-08
  CL-   17      0.000      0.000     A  0.0128097124 6.0466176e-05
 CChl    6      0.000      0.000     A  0.0026308693  4.064256e-06
CLChl   17      0.000      0.000     A  0.0083066819  1.3764842e-05
 HChl    1      0.000      0.000     A  3.76996e-05  4.2999495e-09
SDmso   16      0.000      0.000     A  0.010561673  2.149806e-05
CDmso    6      0.000      0.000     A  0.0096138025  2.6646244e-05
ODmso    8      0.000      0.000     A  0.0022707131  7.5144626e-07
 CCl4    6      0.000      0.000     A  0.0026308693  7.5999462e-06
CLCl4   17      0.000      0.000     A  0.0076040144  1.2767758e-05
 FTFE    9      0.000      0.000     A  0.0011778624       1e-06
 CTFE    6      0.000      0.000     A  0.0023406244  3.374569e-06
CHTFE    6      0.000      0.000     A  0.0071048041  2.5775929e-05
 OTFE    8      0.000      0.000     A  0.0022619536  1.505529e-06
CUrea    6      0.000      0.000     A  0.0048868488  1.3589545e-05
OUrea    8      0.000      0.000     A  0.0023639044  1.5898688e-06
NUrea    7      0.000      0.000     A  0.0033527574  3.9509513e-06
   SI   14      0.000      0.000     A  0.01473796  2.2193521e-05
 MNH3    0      0.000      0.000     A  0.0           0.0
   MW    0      0.000      0.000     D  0.0           0.0
 CH3p    6      0.000      0.000     A  0.0096138025  2.6646244e-05

#define at140_1    0.28504386  0.29879737
"""

gms = GromacsString(temp)
ff = FF()
ff.parse_gm(gms)
print('atom types 1\n', ff.a_type)

# same as
ff = FF()
ff.parse_gm(temp, parse_from_file=False)
print('atom types 2\n', ff.a_type)

# parsing  force fields
fffs = """./oplsaa.ff/forcefield.itp
./amber99sb-ildn.ff/forcefield.itp
./amber03.ff/forcefield.itp
./gromos43a2.ff/forcefield.itp
./gromos43a1.ff/forcefield.itp
./gromos45a3.ff/forcefield.itp
./gromos53a5.ff/forcefield.itp
./gromos54a7.ff/forcefield.itp
./charmm27.ff/forcefield.itp
./gromos53a6.ff/forcefield.itp
./amber94.ff/forcefield.itp
./amber99.ff/forcefield.itp
./amberGS.ff/forcefield.itp
./amber99sb.ff/forcefield.itp
./amber96.ff/forcefield.itp"""

for i in fffs.split():
    ff = parse_ff(i, format_type = 'gm')
    print(i, '# atom types', len(ff.a_type))
    print('')

os.chdir(out_data_fd)
# convert to gromos format
ff = parse_ff('./gromos54a7.ff/forcefield.itp', format_type = 'gm')
ff.gm2gr()
ff.write_ifp('gromacs2gromos_gff.ifp')

fd_gromacs = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromacs')
top_file = os.path.join(fd_gromacs, 'gromos', 'topol.top')

# parse gromacs top file (gromos topology)
top = parse_top(top_file, format_type = 'gm')

# parse gromacs top file using MolSystem (allows for more options/functionality) class (amber topology)
top_file = os.path.join(fd_gromacs, 'amber', 'topol.top')
atop = MolSystem()
atop.parse_top(top_file, allow_multiple_matches = 1, format_type = 'gm')

# parse gromacs top file (charmm and opls topology)
top_file = os.path.join(fd_gromacs, 'charmm', 'topol.top')
ctop = parse_top(top_file, allow_multiple_matches = 1, format_type = 'gm')
top_file = os.path.join(fd_gromacs, 'opls', 'topol.top')
otop = parse_top(top_file, allow_multiple_matches = 1, format_type = 'gm', flag_molsys = 1)

# converting to gromos format
top.gm2gr()
top.write_top('gm2gr_top.top')

atop.gm2gr()
atop.write_top('gm2gr_atop.top')

import glob
for f_path in glob.glob('#bck*'):
    os.remove(f_path)

out_gr2gm_data_fd = os.path.join(out_data_fd, 'gr2gm_top')
os.chdir(out_gr2gm_data_fd)
for f_path in glob.glob('#bck*'):
    os.remove(f_path)

out_gr2gm_data_fd = os.path.join(out_data_fd, 'gr2gm_top_ions')
os.chdir(out_gr2gm_data_fd)
for f_path in glob.glob('#bck*'):
    os.remove(f_path)

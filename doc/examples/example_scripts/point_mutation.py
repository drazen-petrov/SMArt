import sys
import os

smart_fd = os.path.abspath(sys.argv[0])
smart_fd = os.path.split(smart_fd)[0]

for _ in range(3):
    smart_fd = os.path.split(smart_fd)[0]

sys.path.insert(0, smart_fd)

from SMArt.incl import set_print_warnings
set_print_warnings(False)

from SMArt import md
from SMArt import alchemy

# GROMACS
top_file = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromacs', 'gromos', 'topol.top')
top_gm = md.parse_top(top_file, format_type='gm')
top_file = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromacs', 'gromos', 'mut', 'topol.top')
top_gm_mut = md.parse_top(top_file, format_type='gm')
mt, mt_mut = top_gm.molecule_types['Protein_chain_A'], top_gm_mut.molecule_types['Protein_chain_A']

sol = alchemy.point_mutation(mt, mt_mut, ff_dumm = top_gm.get_DUM_type) # flag_top_prune = 'bond' (by default)
sol.toptp.write_itp('toptp.itp', flag_generate_excl_pairs = True)

sol = alchemy.point_mutation(mt, mt_mut, flag_top_prune = None, ff_dumm = top_gm.get_DUM_type)
sol.toptp.write_itp('toptp_bond.itp', flag_generate_excl_pairs = True)


#GROMOS
top_file = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromos', 'villin', 'wt.top')
top_gr = md.parse_top(top_file)
top_file = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromos', 'villin', 'mut.top')
top_gr_mut = md.parse_top(top_file)

sol = alchemy.point_mutation(top_gr, top_gr_mut) # flag_top_prune = 'bond' (by default)
sol.toptp.write_top('toptp.top')
sol.toptp.write_ptp('toptp.ptp')

sol = alchemy.point_mutation(top_gr, top_gr_mut, flag_top_prune = None)
sol.toptp.write_top('toptp_bond.top')
sol.toptp.write_ptp('toptp_bond.ptp')

import glob
for f_path in glob.glob('#bck*'):
    os.remove(f_path)

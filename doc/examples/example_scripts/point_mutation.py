# importing and preparing some variables
import sys
import os
import glob
import re
"""
SMArt_fd = os.path.abspath(os.getcwd())
for _ in range(3):
    SMArt_fd = os.path.split(smart_fd)[0]
sys.path.insert(0, smart_fd)
#"""
# if SMArt not in the path, uncomment the code above

import SMArt
data_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'some_data'))
out_fd = os.path.abspath(os.path.join(os.getcwd(), 'PTP_out_data'))
if not os.path.isdir(out_fd):os.mkdir(out_fd)

from SMArt.md import parse_top
from SMArt import alchemy


# GROMACS
# parse the topologies
top_wt_file = os.path.join(data_fd, 'gromacs', 'LYSH.top')
top_wt = parse_top(top_wt_file, format_type='gm')
top_PTM_file = os.path.join(data_fd, 'gromacs', 'K3C.top')
top_PTM = parse_top(top_PTM_file, format_type='gm')
mt_wt, mt_PTP = top_wt.molecule_types['mol_1'], top_PTM.molecule_types['mol_1']

# run the MCS algorithm
mcs = alchemy.point_mutation(mt_wt, mt_PTP)
sol = alchemy.generate_2state_top(mcs, ff_dumm = top_wt.get_DUM_type) # solution = 0 by default
out_itp = os.path.join(out_fd, 'toptp.itp')
sol.toptp.write_itp(out_itp, flag_generate_excl_pairs = True)

# GROMOS
# parse topologies
top_wt_file = os.path.join(data_fd, 'gromos', 'LYS_K3C_KAC', 'single_AA', 'LYS.top')
top_wt = parse_top(top_wt_file)
top_PTM_file = os.path.join(data_fd, 'gromos', 'LYS_K3C_KAC', 'single_AA', 'K3C.top')
top_PTM = parse_top(top_PTM_file)

# run the MCS algorithm
mcs = alchemy.point_mutation(top_wt, top_PTM)
sol = alchemy.generate_2state_top(mcs) # solution = 0 by default
sol.toptp.write_top(os.path.join(out_fd, 'toptp.top'))
sol.toptp.write_ptp(os.path.join(out_fd, 'toptp.ptp'))


# alternative MCS search - no bond perturbations allowed
mcs = alchemy.point_mutation(top_wt, top_PTM, flag_top_prune = 'bond')
sol = alchemy.generate_2state_top(mcs) # solution = 0 by default
sol.toptp.write_top(os.path.join(out_fd, 'toptp_nobond_ptp.top'))
sol.toptp.write_ptp(os.path.join(out_fd, 'toptp_nobond_ptp.ptp'))


import glob
for f_path in glob.glob(os.path.join(out_fd, '#bck*')):
    os.remove(f_path)

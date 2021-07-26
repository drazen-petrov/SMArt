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
out_fd = os.path.abspath(os.path.join(os.getcwd(), 'EDS_out_data'))
if not os.path.isdir(out_fd):os.mkdir(out_fd)

from SMArt.md import parse_top
from SMArt import alchemy

TOPs_path = os.path.join(data_fd, 'gromos', 'GRA2', 'tops', '*top')
TOPs = glob.glob(TOPs_path)
for fpath in TOPs:print(fpath)
print('\n\nMCS search...')

tops = [parse_top(fpath) for fpath in TOPs]

mcs = alchemy.get_EDS(*tops[:5]) # first 5 tops; faster than all 16 topologies (line below)
#mcs = alchemy.get_EDS(*tops)

sol_EDS, state_names = alchemy.generate_EDS_top(mcs)
sol_EDS.toptp.write_top(os.path.join(out_fd, 'EDS.top'))
sol_EDS.toptp.write_EDS(os.path.join(out_fd, 'EDS.ptp'), state_names=state_names)
print(sol_EDS)

for f_path in glob.glob('#bck*'):
    os.remove(f_path)

import sys
import os

smart_fd = os.path.abspath(sys.argv[0])
smart_fd = os.path.split(smart_fd)[0]

for _ in range(3):
    smart_fd = os.path.split(smart_fd)[0]

sys.path.insert(0, smart_fd)

import SMArt
from SMArt import md
from SMArt import alchemy

SMArt.incl.set_print_warnings(0)

import glob
TOPs_path = os.path.join(smart_fd, 'doc', 'examples', 'some_data', 'gromos', 'GRA2', 'tops', '*top')
print(TOPs_path)
TOPs = glob.glob(TOPs_path)

tops = [md.parse_top(fpath) for fpath in TOPs]
#sol, state_names = alchemy.get_EDS(*tops[:4])
sol, state_names = alchemy.get_EDS(*tops)

sol.toptp.write_top('EDS.top')
sol.toptp.write_EDS('EDS.ptp', state_names=state_names)

for f_path in glob.glob('#bck*'):
    os.remove(f_path)

import os
import glob
import SMArt
from SMArt import alchemy
from SMArt.md import parse_top
from SMArt.md.incl import Dummy

data_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'some_data'))

TOPs_path = os.path.join(data_fd, 'gromos', 'GRA2', 'tops', '*top')
TOPs = glob.glob(TOPs_path)
tops = [parse_top(fpath) for fpath in TOPs]

import unittest

class Test_EDS(unittest.TestCase):
    
    def test_get_EDS_GRA2_5states(self):
        mcs = alchemy.get_EDS(*tops[:5]) # first 5 tops; faster than all 16 topologies
        sol_EDS, state_names = alchemy.generate_EDS_top(mcs)
        assert sol_EDS._sol.shape == (24, 5)
        
    def test_get_EDS_GRA2_5states_core_common_atms(self):
        mcs = alchemy.get_EDS(*tops[:5], flag_get_core_common_atoms_top0=1) # first 5 tops; faster than all 16 topologies
        sol_EDS, state_names = alchemy.generate_EDS_top(mcs)
        assert sol_EDS._sol.shape == (24, 5)

    def test_get_EDS_GRA2_core_common_atms(self):
        mcs = alchemy.get_EDS(*tops, flag_get_core_common_atoms_top0=1) # first 5 tops; faster than all 16 topologies
        sol_EDS, state_names = alchemy.generate_EDS_top(mcs)
        assert sol_EDS._sol.shape == (24, 16)

if __name__ == '__main__':
    unittest.main()

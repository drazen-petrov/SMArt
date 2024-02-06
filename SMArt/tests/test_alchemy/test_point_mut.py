import os
import SMArt
from SMArt import alchemy
from SMArt.md import parse_top
from SMArt.md.incl import Dummy

data_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'some_data'))
# data2compare_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'examples', 'example_scripts', 'PTP_out_data'))
# out_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], 'tests', 'test_alchemy', 'PTP_out_data'))
# if not os.path.isdir(out_fd):os.mkdir(out_fd)

# GROMOS
# parse topologies
top_wt_file = os.path.join(data_fd, 'gromos', 'LYS_K3C_KAC', 'single_AA', 'LYS.top')
top_wt = parse_top(top_wt_file)
top_PTM_file = os.path.join(data_fd, 'gromos', 'LYS_K3C_KAC', 'single_AA', 'K3C.top')
top_PTM = parse_top(top_PTM_file)

# GROMACS
# parse the topologies
top_wt_file = os.path.join(data_fd, 'gromacs', 'LYSH.top')
gm_top_wt = parse_top(top_wt_file, format_type='gm')
top_PTM_file = os.path.join(data_fd, 'gromacs', 'K3C.top')
gm_top_PTM = parse_top(top_PTM_file, format_type='gm')
mt_wt, mt_PTP = gm_top_wt.molecule_types['mol_1'], gm_top_PTM.molecule_types['mol_1']

import unittest

class Test_PointMutation(unittest.TestCase):
    data_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'some_data'))
    # data2compare_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], '..', 'doc', 'examples', 'example_scripts', 'PTP_out_data'))
    # out_fd = os.path.abspath(os.path.join(os.path.split(SMArt.__file__)[0], 'tests', 'test_alchemy', 'PTP_out_data'))
    # if not os.path.isdir(out_fd):os.mkdir(out_fd)
    
    def test_GROMOS(self):
        mcs = alchemy.point_mutation(top_wt, top_PTM)
        sol = alchemy.generate_2state_top(mcs) # solution = 0 by default
        #atoms
        assert len(top_wt.atoms) == len(top_PTM.atoms) == len(sol.toptp.atoms)
        assert len(sol.toptp.ptp_int) == 2
        assert sol.toptp.ptp_int['bonds']==3
        assert sol.toptp.ptp_int['angles']==6 
        assert (sol.toptp.ptp_atom  == [3,4,3]).all()
        # bonds
        c_ptp_int = 0
        for sol_int in sol.toptp.bonds:
            if len(sol_int.ND_states)==2 and None not in sol_int.ND_states:
                c_ptp_int+=1
        assert c_ptp_int == 3
        # angles
        c_ptp_int = 0
        for sol_int in sol.toptp.angles:
            if len(sol_int.ND_states)==2 and None not in sol_int.ND_states:
                c_ptp_int+=1
        assert c_ptp_int == 6
        
    def test_GROMOS_no_ptpbond(self):#no bond perturbations allowed
        mcs = alchemy.point_mutation(top_wt, top_PTM, flag_top_prune = 'bond')
        sol = alchemy.generate_2state_top(mcs) # solution = 0 by default
        # atoms
        assert len(top_PTM.atoms) + 3 == len(sol.toptp.atoms)
        assert sol.toptp.ptp_int == {}
        assert (sol.toptp.ptp_atom == [0,7,6]).all()
        c_dummies = [0]*2
        for r_i, row in sol._sol.iterrows():
            for i in range(2):
                if row[i] is Dummy:
                    c_dummies[i]+=1
        assert c_dummies == [3,3]
        # bonds
        c_ptp_int = [0,0]
        for sol_int in sol.toptp.bonds:
            for i in range(2):
                if sol_int.int_states[i] is None:
                    c_ptp_int[i]+=1
        assert c_ptp_int == [3,3]
        # angles
        c_ptp_int = [0,0]
        for sol_int in sol.toptp.angles:
            for i in range(2):
                if sol_int.int_states[i] is None:
                    c_ptp_int[i]+=1
        assert c_ptp_int == [6,6]
    
    def test_GROMACS(self):
        mcs = alchemy.point_mutation(mt_wt, mt_PTP)
        sol = alchemy.generate_2state_top(mcs, ff_dumm = top_wt.get_DUM_type) # solution = 0 by default
        # atoms
        assert len(mt_wt.atoms) == len(mt_PTP.atoms) == len(sol.toptp.atoms)
        assert sol.toptp.ptp_int['bonds']==3
        assert sol.toptp.ptp_int['angles']==6 
        assert (sol.toptp.ptp_atom  == [3,4,3]).all()
        #bonds
        c_ptp_int = 0
        for sol_int in sol.toptp.bonds:
            if len(sol_int.ND_states)==2 and None not in sol_int.ND_states:
                c_ptp_int+=1
        assert c_ptp_int == 3
        # angles
        c_ptp_int = 0
        for sol_int in sol.toptp.angles:
            if len(sol_int.ND_states)==2 and None not in sol_int.ND_states:
                c_ptp_int+=1
        assert c_ptp_int == 6

if __name__ == '__main__':
    unittest.main()

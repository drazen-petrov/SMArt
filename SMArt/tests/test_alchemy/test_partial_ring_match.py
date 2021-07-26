from SMArt.incl import np
from SMArt.graph import Graph
from SMArt.alchemy.MCS import MCS

from SMArt.tests.test_graph.gen_test_graphs import *

import unittest

class Test_MCS_max_partial_ring_match(unittest.TestCase):
    G_622 = get_rc(6,2,2)
    G_621 = get_rc(6,2,1)
    G_522 = get_rc(5,2,2)
    G_521 = get_rc(5,2,1)
    G7 = get_c(7)
    G_56 = comb_G(get_r(6, 1), get_r(5, 5))
    g622 = Graph(G_622)
    g621 = Graph(G_621)
    g522 = Graph(G_522)
    g521 = Graph(G_521)
    g7 = Graph(G7)
    g56 = Graph(G_56)

    def test_chain_ring_max_PRM0(self):
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g622, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (15, 2), 'should have matched 2 atoms\n'+str(sol))
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 1
        mcs = MCS(self.g7, self.g621, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (15, 2), 'should have matched 2 atoms\n'+str(sol))
        # test chain_7 to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g522, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (14, 2), 'should have matched 2 atoms\n'+str(sol))

    def test_chain_ring_max_PRM1(self):
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g622, max_partial_ring_match=1)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (14, 2), 'should have matched 3 atoms\n'+str(sol))
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 1
        mcs = MCS(self.g7, self.g621, max_partial_ring_match=1)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (12, 2), 'should have matched 5 atoms\n'+str(sol))
        # test chain_7 to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g522, max_partial_ring_match=1)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (13, 2), 'should have matched 3 atoms\n'+str(sol))

    def test_chain_ring_max_PRM2(self):
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g622, max_partial_ring_match=2)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (11, 2), 'should have matched 6 atoms\n'+str(sol))
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 1
        mcs = MCS(self.g7, self.g621, max_partial_ring_match=2)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (12, 2), 'should have matched 5 atoms\n'+str(sol))
        # test chain_7 to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g7, self.g522, max_partial_ring_match=2)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (10, 2), 'should have matched 6 atoms\n'+str(sol))

    def test_ring_ring_max_PRM0(self):
        # test 6-ring to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g622, self.g522, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (17, 2), 'should have matched 2 atoms\n'+str(sol))
        # test 6-ring with 2 chains of len 2 and overlap 2 to 5-ring with overlap 1
        mcs = MCS(self.g622, self.g521, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (17, 2), 'should have matched 2 atoms\n'+str(sol))
        # test 6-ring to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g622, self.g622, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (10, 2), 'should have matched all atoms\n'+str(sol))

    def test_ring_ring_max_PRM1(self):
        # test 6-ring to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g622, self.g522, max_partial_ring_match=1)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (16, 2), 'should have matched 3 atoms\n'+str(sol))
        # test 6-ring with 2 chains of len 2 and overlap 2 to 5-ring with overlap 1
        mcs = MCS(self.g622, self.g521, max_partial_ring_match=1)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (16, 2), 'should have matched 3 atoms\n'+str(sol))

    def test_ring_ring_max_PRM2(self):
        # test chain_7 to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g622, self.g522, max_partial_ring_match=2)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (13, 2), 'should have matched 6 atoms\n'+str(sol))
        # test 6-ring with 2 chains of len 2 and overlap 2 to 5-ring with overlap 1
        mcs = MCS(self.g622, self.g521, max_partial_ring_match=2)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (14, 2), 'should have matched 5 atoms\n'+str(sol))

    def test_ring_ring_partial(self):
        # test 56-ring to 6-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g56, self.g622, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (13, 2), 'should have matched 6 atoms\n'+str(sol))
        # test 56-ring to 5-ring with 2 chains of len 2 and overlap 2
        mcs = MCS(self.g56, self.g522, max_partial_ring_match=0)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (13, 2), 'should have matched 5 atoms\n'+str(sol))
        # test 56-ring to 6-ring with 2 chains of len 2 and overlap 2 (no partial ring match)
        mcs = MCS(self.g56, self.g622, max_partial_ring_match=0, flag_partial_ring=False)
        mcs.enumerate_stepwise_sorted()
        #sol = mcs.solutions[0]
        self.assertEqual(len(mcs.solutions), 0, 'should be no solutions, but something found')
        # test 56-ring to 6-ring with 2 chains of len 2 and overlap 2 (no partial ring match)
        mcs = MCS(self.g56, self.g622, max_partial_ring_match=1, flag_partial_ring=False)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (18, 2), 'should have matched 1 atoms\n'+str(sol))
        # test 56-ring to 6-ring with 2 chains of len 2 and overlap 2 (no partial ring match)
        mcs = MCS(self.g56, self.g622, max_partial_ring_match=2, flag_partial_ring=False)
        mcs.enumerate_stepwise_sorted()
        sol = mcs.solutions[0]
        self.assertEqual(sol._sol.shape, (17, 2), 'should have matched 2 atoms\n'+str(sol))


if __name__ == '__main__':
    unittest.main()
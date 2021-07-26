from SMArt.incl import np
from SMArt.graph import Graph
from SMArt.alchemy.MCS import MCS

from SMArt.tests.test_graph.gen_test_graphs import *

import unittest


class Test_MCS_RMSD(unittest.TestCase):
    N = 6
    G = get_r(N)
    g1 = Graph(G)
    g2 = Graph(G)

    def __get_coord(self, offset=0):
        coord = np.arange(18).reshape(self.N,3)
        new_coord = []
        for i in range(coord.shape[0]):
            temp_i = (i + offset) % coord.shape[0]
            new_coord.append(coord[temp_i])
        return new_coord

    def test_RMDS_simpe(self):
        for offset in range(self.N):
            coords = [self.__get_coord(), self.__get_coord(offset)]
            m = MCS(self.g1, self.g2, add_RMSD='simple', coords = coords)
            m.enumerate_stepwise_sorted()
            self.assertTrue(((m.solutions[0]._sol[1] + offset) % self.N == m.solutions[0]._sol[0]).all())

    def test_RMDS_pairwise(self):
        for offset in range(self.N):
            coords = [self.__get_coord(), self.__get_coord(offset)]
            m = MCS(self.g1, self.g2, add_RMSD='pairwise', coords = coords)
            m.enumerate_stepwise_sorted()
            self.assertTrue(((m.solutions[0]._sol[1] + offset) % self.N == m.solutions[0]._sol[0]).all())


if __name__ == '__main__':
    unittest.main()
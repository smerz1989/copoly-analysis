import unittest
import os, sys
import subprocess as sb
import numpy as np

sys.path.append('../')

from dump import dump
import trajectory_class as tjc

class TestSimulationSnapshot(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestSimulationSnapshot, cls).setUpClass()
        cls.trajfileAABB = os.path.abspath('./test_files/simulation_snapshots/lts/poly.lammpstrj')
        cls.bondsfileAABB = os.path.abspath('./test_files/simulation_snapshots/lts/bonddump.dump')
        cls.dataAABB = os.path.abspath('./test_files/simulation_snapshots/lts/system.data')
        cls.trajAABB = dump(cls.trajfileAABB)
        cls.bondsAABB = dump(cls.bondsfileAABB)
        cls.snapshots_AABB = tjc.construct_molecule_trajectory(cls.dataAABB,cls.bondsAABB,cls.trajAABB)
        cls.timestep1, cls.snapshot1_AABB = next(cls.snapshots_AABB)
        cls.trajfileABAB = os.path.abspath('./test_files/simulation_snapshots/ABAB/poly.lammpstrj')
        cls.bondsfileABAB = os.path.abspath('./test_files/simulation_snapshots/ABAB/bonddump.dump')
        cls.dataABAB = os.path.abspath('./test_files/simulation_snapshots/ABAB/system.data')
        cls.trajABAB = dump(cls.trajfileABAB)
        cls.bondsABAB = dump(cls.bondsfileABAB)
        cls.snapshots_ABAB = tjc.construct_molecule_trajectory(cls.dataABAB,cls.bondsABAB,cls.trajABAB)
        cls.timestep1, cls.snapshot1_ABAB = next(cls.snapshots_ABAB)
   

    def test_get_sequences_returns_correct_sequences_for_all_A_and_all_B_polymers(self):
        sequences = self.snapshot1_AABB.get_sequences()
        chain = self.snapshot1_AABB.topology.bfsiter(self.snapshot1_AABB.topology.vs[0])
        str_sequence = self.snapshot1_AABB.get_chain_sequence(chain)
        expected_sequence = '8 2 6 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5 5 3 5' 
        self.assertEqual(str_sequence,expected_sequence)
        

    def test_get_sequence_gets_correct_sequence_probs_for_all_A_and_all_B_polymers(self):
        pAA,pBB,pAB = self.snapshot1_AABB.get_all_probs()
        np.testing.assert_almost_equal([pAA,pBB,pAB],[0.5,0.5,0.],err_msg="get_all_probs does not return corect pAA,pBB,pAB for uniform (all A and all B) polymers")

    def test_get_sequence_gets_correct_sequence_probs_for_ABAB_polymers(self):
        pAA,pBB,pAB = self.snapshot1_ABAB.get_all_probs()
        np.testing.assert_almost_equal([pAA,pBB,pAB],[0.0,0.0,1.],err_msg="get_all_probs does not return corect pAA,pBB,pAB for uniform (all A and all B) polymers")


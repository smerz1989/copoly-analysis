import unittest
import os, sys
import subprocess as sb
import numpy as np

#sys.path.append('../')
from copoly import dump_generator as dumpg
from copoly.dump import dump
import copoly.trajectory_class as tjc

class TestPotentials(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestPotentials, cls).setUpClass()
        test_location = os.path.dirname(os.path.abspath(__file__))
        cls.lj_expand_twopiece_simulation_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0/system.in')) 
        cls.lj_expand_twopiece_simulation_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0/pair_energy.data'))   
        cls.lj_expand_twopiece_simulation_AB_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0_AB/system.in')) 
        cls.lj_expand_twopiece_simulation_AB_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0_AB/pair_energy.data'))   
        cls.lj_expand_twopiece_simulation_AA_eps5_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0_AA_eps5/system.in')) 
        cls.lj_expand_twopiece_simulation_AA_eps5_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift0_AA_eps5/pair_energy.data'))   
        cls.lj_twopiece_simulation_AA_eps5_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_twopiece/shift0_AA_eps5/system.in')) 
        cls.lj_twopiece_simulation_AA_eps5_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_twopiece/shift0_AA_eps5/pair_energy.data'))   
        cls.lj_expand_twopiece_oneshift_simulation_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift1/system.in')) 
        cls.lj_expand_twopiece_oneshift_simulation_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_expand_twopiece/shift1/pair_energy.data'))   
        cls.lj_twopiece_simulation_input_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_twopiece/shift0_comparison/system.in')) 
        cls.lj_twopiece_simulation_data_file = os.path.abspath(os.path.join(test_location,'test_files/potential_simulations/lj_twopiece/shift0_comparison/pair_energy.data'))   


    @staticmethod
    def expected_lj_expand_twopiece(r,sigma,eps11,eps_rep,shift=0.):
        r0 = (2**(1/6.))*sigma
        energy = np.zeros(len(r))
        r_shifted = r-shift
        r_rep = r_shifted[r_shifted<r0]
        r_attr = r_shifted[r_shifted>=r0]
        energy[r_shifted<r0] = 4*eps_rep*((sigma/r_rep)**12-(sigma/r_rep)**6)+eps_rep-eps11
        energy[r_shifted>=r0] = 4*eps11*((sigma/r_attr)**12-(sigma/r_attr)**6)
        energy[np.isnan(energy)]=np.inf
        return(energy)
        
    @staticmethod
    def datafile_to_array(filepath):
        dg = dumpg.read_dump(filepath)
        num_timesteps = dumpg.get_number_of_timesteps(filepath) 
        energy_data = np.empty((num_timesteps,3))
        for entry,i in zip(dg,range(num_timesteps)):
            energy_data[i,:] = entry[1]
        return(energy_data)

    @staticmethod
    def run_lammps_sim(input_file_path,lammps_path='lmp'):
        current_place = os.path.abspath('.')
        os.chdir(os.path.abspath(os.path.dirname(input_file_path)))
        sb.call([lammps_path,"-i",input_file_path],stdout=open('tmp.out','w'))
        os.chdir(current_place)

    def test_lj_expand_twopiece_returns_correct_energy_vs_distance(self):
        self.run_lammps_sim(self.lj_expand_twopiece_simulation_input_file)
        energy_data = self.datafile_to_array(self.lj_expand_twopiece_simulation_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=1.0,eps_rep=1.0,shift=0.)
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,err_msg="lj/expand/twopiece does not produce expected analytical result for shift=0, sigma=0.5,eps11=1.0, and eps_rep=1.0")

    def test_lj_expand_twopiece_returns_correct_energy_vs_distance_with_two_different_types_of_monomers(self):
        self.run_lammps_sim(self.lj_expand_twopiece_simulation_AB_input_file)
        energy_data = self.datafile_to_array(self.lj_expand_twopiece_simulation_AB_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=1.0,eps_rep=1.0,shift=0.)
        test_data = np.empty((energy_data.shape[0],4))
        test_data[:,0:3]=energy_data
        test_data[:,3]=analytical_data
        np.savetxt('test_data_AB.txt',test_data) 
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,err_msg="lj/expand/twopiece does not produce expected analytical result for shift=0, sigma=0.5,eps11=1.0, and eps_rep=1.0")

    @unittest.skip("Skipping this test until lj/expand/twopiece potential is fixed")
    def test_lj_expand_twopiece_returns_correct_energy_vs_distance_with_different_attractive_and_repulsive_epsilon(self):
        self.run_lammps_sim(self.lj_expand_twopiece_simulation_AA_eps5_input_file)
        energy_data = self.datafile_to_array(self.lj_expand_twopiece_simulation_AA_eps5_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=5.0,eps_rep=1.0,shift=0.)
        test_data = np.empty((energy_data.shape[0],4))
        test_data[:,0:3]=energy_data
        test_data[:,3]=analytical_data
        np.savetxt('test_data_AA_eps5.txt',test_data) 
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,err_msg="lj/expand/twopiece does not produce expected analytical result for shift=0, sigma=0.5,eps11=1.0, and eps_rep=1.0")

    def test_lj_twopiece_returns_correct_energy_vs_distance_with_different_attractive_and_repulsive_epsilon(self):
        self.run_lammps_sim(self.lj_twopiece_simulation_AA_eps5_input_file)
        energy_data = self.datafile_to_array(self.lj_twopiece_simulation_AA_eps5_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=5.0,eps_rep=1.0,shift=0.)
        test_data = np.empty((energy_data.shape[0],4))
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,err_msg="lj/expand/twopiece does not produce expected analytical result for shift=0, sigma=0.5,eps11=1.0, and eps_rep=1.0")

    def test_lj_expand_twopiece_returns_correct_energy_vs_distance_with_shift(self):
        self.run_lammps_sim(self.lj_expand_twopiece_oneshift_simulation_input_file)
        energy_data = self.datafile_to_array(self.lj_expand_twopiece_oneshift_simulation_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=1.0,eps_rep=1.0,shift=1.)
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,verbose=True,err_msg="lj/expand/twopiece does not produce expected analytical result for shift=1, sigma=0.5,eps11=1.0, and eps_rep=1.0")

    def test_lj_twopiece_returns_correct_energy_vs_distance(self):
        self.run_lammps_sim(self.lj_twopiece_simulation_input_file)
        energy_data = self.datafile_to_array(self.lj_twopiece_simulation_data_file)
        analytical_data = self.expected_lj_expand_twopiece(energy_data[:,1],sigma=0.5,eps11=1.0,eps_rep=1.0,shift=0.)
        np.testing.assert_allclose(energy_data[:,2],analytical_data,rtol=1e-5,atol=1e-5,verbose=True,err_msg="lj/twopiece does not produce expected analytical result for shift=1, sigma=0.5,eps11=1.0, and eps_rep=1.0")


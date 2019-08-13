import os
import subprocess as sb

class Simulation(object):
    def __init__(self,total_monomers=3000,monomer_A_fraction=0.5,p=0.9,dump_frequency=1000,lt_dir = os.path.abspath('../../lt_files/'),xyz_dir = os.path.abspath('../../xyzs/')):
        self.total_monomers = total_monomers
        self.monomer_A_fraction = monomer_A_fraction
        self.p = p
        self.dump_frequency = dump_frequency
        self.lt_dir = lt_dir
        self.xyz_dir = xyz_dir

    def compile_simulation(self,packmol_path=os.path.abspath('/scratch/snm8xf/packmol/packmol')):
        self.change_monomer_count()
        os.chdir(self.xyz_dir)
        sb.call([packmol_path],stdin=open('np.inp'))
        os.chdir(self.lt_dir)
        sb.call(["moltemplate.sh","-xyz",self.xyz_dir+"/np.xyz","-atomstyle","angle","system.lt"])

    def change_monomer_count(self):
        num_monA = int(self.total_monomers*self.monomer_A_fraction)
        num_monB = self.total_monomers-num_monA
        os.chdir(self.xyz_dir)
        sb.call(["sed","-i",'/abead/,/bbead/ s/number\ [0-9]\+/number\ '+str(num_monA)+'/g',"np.inp"])
        sb.call(["sed","-i",'/bbead/, s/number\ [0-9]\+/number\ '+str(num_monB)+'/g',"np.inp"])


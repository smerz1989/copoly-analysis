import os
import subprocess as sb
import glob
import shutil

class Simulation(object):
    def __init__(self,total_monomers=3000,monomer_A_fraction=0.5,monomer_attractions=(1,1,1),p=0.9,
                    dump_frequency=1000,lt_dir = os.path.abspath('../../lt_files/'),
                    xyz_dir = os.path.abspath('../../xyzs/')):
        self.total_monomers = total_monomers
        self.monomer_A_fraction = monomer_A_fraction
        self.p = p
        self.dump_frequency = dump_frequency
        self.lt_dir = lt_dir
        self.xyz_dir = xyz_dir
        self.eAA,self.eBB,self.eAB = monomer_attractions

    def compile_simulation(self,packmol_path=os.path.abspath('/scratch/snm8xf/packmol/packmol')):
        self.change_monomer_count()
        self.change_monomer_attraction()
        os.chdir(self.xyz_dir)
        sb.call([packmol_path],stdin=open('np.inp'))
        os.chdir(self.lt_dir)
        sb.call(["moltemplate.sh","-xyz",self.xyz_dir+"/np.xyz","-atomstyle","angle","system.lt"])
        sb.call(["sed","-i",'s/a\\"/a\\"\ extra\/special\/per\/atom\ 4\ extra\/bond\/per\/atom\ 2\ extra\/angle\/per\/atom\ 2/g',"system.in"])
        sb.call(["sed","-i",'s/\!\(.*\)\!/\$\{\\1\}/g',"system.in.run"])
        sb.call(["sed","-i",'s/\!(\(.*\))/\$(\\1)/g',"system.in.run"])
    
    def change_monomer_attraction(self):
        cur_path = os.path.abspath('.')
        os.chdir(self.lt_dir)
        sb.call(["sed","-i",'/\@atom\:A\ \@atom\:A/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eAA)+'/g','copolyff.lt'])
        sb.call(["sed","-i",'/\@atom\:A\ \@atom\:B/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eAB)+'/g','copolyff.lt'])
        sb.call(["sed","-i",'/\@atom\:B\ \@atom\:B/ s/twopiece\ [0-9]\?\.\?[0-9]\?/twopiece\ '+str(self.eBB)+'/g','copolyff.lt'])
        os.chdir(cur_path)


    def change_monomer_count(self):
        curr_dir = os.path.abspath('.')
        num_monA = int(self.total_monomers*self.monomer_A_fraction)
        num_monB = self.total_monomers-num_monA
        os.chdir(self.xyz_dir)
        sb.call(["sed","-i",'/abead/,/bbead/ s/number\ [0-9]\+/number\ '+str(num_monA)+'/g',"np.inp"])
        sb.call(["sed","-i",'/bbead/,+2 s/number\ [0-9]\+/number\ '+str(num_monB)+'/g',"np.inp"])
        os.chdir(self.lt_dir)
        sb.call(["sed",'-i','s/ABEAD\ \[[0-9]\+\]/ABEAD\ \['+str(num_monA)+'\]/g',"system.lt"])
        sb.call(["sed",'-i','s/BBEAD\ \[[0-9]\+\]/BBEAD\ \['+str(num_monB)+'\]/g',"system.lt"])
        os.chdir(curr_dir)

    def move_simulation_files(self,dest_dir):
        dest_folder = os.path.abspath(dest_dir+'/copoly_{}monomers_{}percentA_{}epsAA_{}epsBB_{}epsAB'.format(self.total_monomers,
                                                                                                 int(100*self.monomer_A_fraction),
                                                                                                    self.eAA,self.eBB,self.eAB))
        self.dest_folder = dest_folder
        if not os.path.exists(dest_folder):
            os.mkdir(dest_folder)
        for simfile in glob.glob(r''+self.lt_dir+'/system.*'):
            shutil.copy(simfile,os.path.abspath(dest_folder))
        for simfile in glob.glob(r''+self.lt_dir+'/*.txt'):
            shutil.copy(simfile,os.path.abspath(dest_folder))
        shutil.copy("submit.sbatch",dest_folder)

    def analyze_simulation(self):
        print("placeholder")
    
    def start_simulation(self):
        os.chdir(self.dest_folder)
        sb.call(["sbatch","submit.sbatch"])

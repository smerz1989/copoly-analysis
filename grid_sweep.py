import subprocess as sb
geometry = 'no_nanoparticle'
geometry_name = geometry.capitalize() if not ("_" in geometry) else '_'.join([part.capitalize() for part in geometry.split('_')])
num_monomers=5000
fA=0.5
anglestrength=20
if anglestrength==20:
    anglefolder = "" 
elif anglestrength>20:
    anglefolder="/Increased_Strength_{}kbT/".format(int(anglestrength))
else:
    anglefolder="/Decreased_Strength_{}kbT/".format(int(anglestrength))
 
epsilons = [0./3.,1/3.,2/3.,3./3.,4./3.,5./3.]
dest_folder = '/scratch/snm8xf/NP-Copoly/{}/LJ_Twopiece/Temperature1kbT/Equilibration_Time/Decreased_Concentration{}'.format(geometry_name,anglefolder)
lt_files = '~/Template_Files/'+geometry+'/lt_files'
xyz_files = '~/Template_Files/'+geometry+'/xyzs'
p=0.9
for epsilon in epsilons:
  print("python ./copoly/create_simulation.py -N {} -fA {} --epsAA {} --epsBB {} --epsAB 0. --angle {} -p {} --xyz {} --lt {} -f {} --slurm --send_to_cluster".format(num_monomers,fA,epsilon,epsilon,anglestrength,p,xyz_files,lt_files,dest_folder))
  for trial in range(1,4):
      sb.call(["python",'./copoly/create_simulation.py','-N',str(num_monomers),
                                                        '-fA',str(fA),'--epsAA',str(epsilon),'--epsBB',str(epsilon),'--epsAB','0.',
                                                         '--angle',str(anglestrength),'-p',str(p),'--xyz',xyz_files,'--lt',lt_files,
                                                          '-f',dest_folder,'--slurm','--send_to_cluster','--suffix','_trial'+str(trial)])

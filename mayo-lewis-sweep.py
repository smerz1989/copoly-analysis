import subprocess as sb
geometry = 'icosahedron'
geometry_name = geometry.capitalize() if not ("_" in geometry) else '_'.join([part.capitalize() for part in geometry.split('_')])
num_monomers=5000
fAs=[0.15,0.35,0.65,0.85]
anglestrength=20
epsilon = 0.
if anglestrength==20:
    anglefolder = "" 
elif anglestrength>20:
    anglefolder="Increased_Strength_{}kbT/".format(int(anglestrength))
else:
    anglefolder="Decreased_Strength_{}kbT/".format(int(anglestrength))
 
dest_folder = '/scratch/snm8xf/NP-Copoly/{}/LJ_Twopiece/{}Mayo-Lewis-Analysis'.format(geometry_name,anglefolder)
lt_files = 'copoly/'+geometry+'/lt_files'
xyz_files = 'copoly/'+geometry+'/xyzs'
p=0.9
for fA in fAs:
  print("python ./copoly/create_simulation.py -N {} -fA {} --epsAA {} --epsBB {} --epsAB 0. --angle {} -p {} --xyz {} --lt {} -f {} --slurm --send_to_cluster".format(num_monomers,fA,epsilon,epsilon,anglestrength,p,xyz_files,lt_files,dest_folder))
  for trial in range(1,4):
      sb.call(["python",'./copoly/create_simulation.py','-N',str(num_monomers),
                                                        '-fA',str(fA),'--epsAA',str(epsilon),'--epsBB',str(epsilon),'--epsAB','0.',
                                                         '--angle',str(anglestrength),'-p',str(p),'--xyz',xyz_files,'--lt',lt_files,
                                                          '-f',dest_folder,'--slurm','--send_to_cluster','--suffix','_trial'+str(trial)])

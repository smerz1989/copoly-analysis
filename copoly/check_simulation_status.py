import simulation_class as sim
import argparse
import os
import glob
import subprocess as sb
import re


class CopolyJob(object):
    def __init__(self,directory):
        print("Placeholder")


def get_job_statuses(slurm_numbers):
    job_statuses = sb.check_output(["squeue","-j",",".join([jobid for jobid in slurm_numbers])]).decode('utf-8').split('\n')
    return(job_status) 


parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-f','--folder',dest='dest_folder',help="Top directory containing multiple simulation folders")

args = parser.parse_args()

folders = glob.glob(args.dest_folder+'/copoly_*')

slurm_numbers=[]
number_atoms=[]
timesteps=[]
conversions = []
for folder in folders:
    slurm_file_name = os.path.split(glob.glob(folder+'/slurm*')[0])[1]
    print("Analyzing slurm file {}".format(slurm_file_name))
    slurm_numbers.append(re.findall(r'[0-9]+',slurm_file_name)[0])
    traj_file = glob.glob(folder+'/*.lammpstrj')[0]
    num_atoms = sb.check_output(["sed",'-n',r'/NUMBER\ OF\ ATOMS/,+1 s/\([0-9]\+\)\+/\1/p',traj_file]).decode("utf-8")
    number_atoms = float(num_atoms.split('\n')[0])
    num_traj_lines = float(sb.check_output(["wc","-l",traj_file]).decode('utf-8').split()[0])
    timesteps.append(int(1000*num_traj_lines/(number_atoms+9)))
    bond_file = glob.glob(folder+'/*.dump')[0]
    num_bonds = [int(num_bonds) for num_bonds in sb.check_output(["sed",'-n',r'/NUMBER/,/BOX/ s/\([0-9]\+\)\+/\1/p',bond_file]).decode("utf-8").split('\n')[:-1]]
    conversions.append((num_bonds[-1]-num_bonds[0])/(num_bonds[0]/2))
    


output = sb.check_output(["squeue","-j",",".join([jobid for jobid in slurm_numbers])])

job_statuses = output.decode("utf-8")[:-1].split('\n')
job_statuses.sort()
new_statuses = ['{}\t{}\t{}'.format(status,timesteps[i],conversions[i],folders[i]) for i,status in enumerate(job_statuses[1:])]
new_statuses.insert(0,job_statuses[0]+'\t'+'Timestep\tConversion\tFolder')
print("\n".join(new_statuses))

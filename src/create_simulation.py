import simulation_class as smc
import argparse
import os

parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-N','--total_monomers',dest='total_monomers',type=int,default=3000)
parser.add_argument('-f','--folder',dest='dest_folder')

args = parser.parse_args()

project_path = os.path.abspath(args.dest_folder)

sim = smc.Simulation(total_monomers=args.total_monomers)
sim.compile_simulation()
sim.move_simulation_files(project_path)

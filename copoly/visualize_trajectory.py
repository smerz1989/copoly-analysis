import copoly
from copoly.dump import dump
#from copoly.dump_generator import read_dump
import copoly.trajectory_class as tjc
import os,sys
import matplotlib.pyplot as plt
import seaborn as sns
import re
sys.path.append('../')
from copoly.dump_generator import read_dump
from tqdm import tqdm
import itertools as it
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Script for downloading and visualizing LAMMPS trajectory using OVITO')
parser.add_argument('-f','--folder',dest='folder', help="Parent directory where simulations are kept")
parser.add_argument('-d','--dest_folder',dest='dest_folder',help="Location where results are placed.")
parser.add_argument('-s','--specific',default=None,dest='specific_file',help="Name of specific folder to analyze if only one folder is wanted.")
args = parser.parse_args()

data_folder = os.path.abspath(args.folder)
data_file = data_folder+'/system.data'
traj_file = data_folder+'/atom_trj.lammpstrj'
bond_file = data_folder+'/bonddump.dump'

snapshots = tjc.construct_molecule_trajectory_from_generators(data_file,bond_file,traj_file)

for i,snapshot in tqdm(enumerate(snapshots),total=1286):
    if i%100==0:
        snapshot.visualize_snapshot(args.dest_folder+'/snapshot_timestep{}.png'.format(i))

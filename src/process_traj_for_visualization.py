import trajectory_class as tc
import dump as dump
import mdtraj
import argparse
import pandas as pd
import nglview
from nglview.contrib.movie import MovieMaker

parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-f','--filename',dest='filename', help="Parent directory where simulations are kept")
parser.add_argument('-d','--dest_folder',dest='dest_folder',help="Location where results are placed.")
args = parser.parse_args()


def topology_from_LAMMPS_data_file(filename,type_to_element_mapping={'1':'Ag','2':'C','3':'C','4':'C','5':'H','6':'H','7':'H','8':'H'}):
    atoms,atom_array = tc.loadAtoms(filename) 
    bonds_array = tc.loadBonds(filename)
    snapshot = tc.SimulationSnapshot(atoms,bonds_array[:,1:4],atom_array)
    atom_dataframe = pd.DataFrame(atom_array[:,[0,2,2,1,1,1]],columns=['serial','name','element','resSeq','resName','chainID'])
    atom_dataframe = atom_dataframe.astype({'serial':'int64','name':'str','element':'str','resSeq':'int64','resName':'str','chainID':'int64'})
    atom_dataframe['element'] = atom_dataframe['element'].map(type_to_element_mapping)
    topology = mdtraj.Topology.from_dataframe(atom_dataframe,bonds_array[:,2:4]-1)
    return(topology)          


#import pdb;pdb.set_trace()
topology = topology_from_LAMMPS_data_file(args.filename)
traj = mdtraj.load('atom_trj.lammpstrj',top=topology)

traj_view = nglview.show_mdtraj(traj)
MovieMaker(traj_view,output='my.gif',in_memory=True).make()

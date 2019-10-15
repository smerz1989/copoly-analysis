import compile_simulation_results as csr
import server_class as svc
import argparse
import os
import pandas as pd
import re
import itertools as it
import json

parser = argparse.ArgumentParser(description='Script for modifying, compiling lt and xyz template files and copying the LAMMPS input script output into a standalone simulation folder')
parser.add_argument('-f','--folder',dest='folder', help="Parent directory where simulations are kept")
parser.add_argument('-d','--dest_folder',dest='dest_folder',help="Location where results are placed.")
parser.add_argument('-s','--specific',default=None,dest='specific_file',help="Name of specific folder to analyze if only one folder is wanted.")
parser.add_argument('-l','--local',dest='local',action='store_true',help='Set this to analyze simulation folders on the local machine')
parser.add_argument('-c','--chain-only',dest='only_chain',action='store_true',help='Analyze just chain block lengths.')
args = parser.parse_args()

if not args.local:
    server_connection = svc.ServerConnection()
    subfolders = server_connection.lsdir(args.folder)
else:
    subfolders = os.listdir(args.folder)
print("Found subfolders {}".format(subfolders))

if args.specific_file==None:
    simfolders = [folder for folder in subfolders if "copoly_" in folder]
    print("Found simulation folders {}".format(simfolders))
else:
    simfolders = [folder for folder in subfolders if args.specific_file in folder]
    print("Found simulation folders {}".format(simfolders))

dest_path= args.dest_folder

for i,folder in enumerate(simfolders):
    print("Analyzing {} simulation of a total of {}".format(i,len(simfolders)))
    os.makedirs(dest_path+'/'+folder,exist_ok=True)  
    if args.local:
        result = csr.SimulationResults(args.folder+'/'+folder,dest_path+'/'+folder,is_remote=False)
    else:
        print("Args.folder is {}".format(args.folder))
        result = csr.SimulationResults(args.folder+'/'+folder,dest_path+'/'+folder)
    result.get_trajectory(dest_path+'/'+folder)
    if not args.only_chain:
        print(dest_path+'/'+folder)
        data = result.analyze_trajectory(dest_path+'/'+folder)
        data.to_csv(dest_path+'/'+folder+'/traj_analysis.csv')
        result.plot_trajectory(data)
    else:
        data = result.analyze_trajectory_by_function(dest_path+'/'+folder)
        with open(dest_path+'/'+folder+'/seq_analysis.json','w') as jfile:
            json.dump(data,jfile,sort_keys=True,indent=2)


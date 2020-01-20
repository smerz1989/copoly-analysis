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
parser.add_argument('--skip',dest='skip',default=None,help="A comma separated list of folders to skip.")
parser.add_argument('-l','--local',dest='local',action='store_true',help='Set this to analyze simulation folders on the local machine')
parser.add_argument('-c','--chain-only',dest='only_chain',action='store_true',help='Analyze just chain block lengths.')
parser.add_argument('--nooverwrite',dest='no_overwrite',action='store_true',help='If flag is used script will not download or analyze simulation data of folders that are already present locally')
parser.add_argument('--no-overwrite-chains',dest='no_overwrite_chains',action='store_true',help='If flag is used script will not download or analyze simulation data of folders that have seq_analysis.json locally')

args = parser.parse_args()

def check_for_file(folder,filename):
    files = os.listdir(folder)
    return(filename in files)

if not args.local:
    server_connection = svc.ServerConnection()
    subfolders = server_connection.lsdir(args.folder)
else:
    subfolders = os.listdir(os.path.abspath(args.folder))
print("Found subfolders {}".format(subfolders))

if args.specific_file==None:
    skipfolders = [] if args.skip==None else args.skip.split(',')
    print("Skipping folders {}".format(skipfolders))
    simfolders = [folder for folder in subfolders if "copoly_" in folder and (not folder in skipfolders)]
    print("Found simulation folders {}".format(simfolders))
else:
    simfolders = [folder for folder in subfolders if (args.specific_file in folder)]
    print("Found simulation folders {}".format(simfolders))

if args.no_overwrite:
    local_folders = os.listdir(os.path.abspath(args.dest_folder))
    copoly_local_folders = [folder for folder in local_folders if "copoly_" in folder]
    simfolders = [folder for folder in simfolders if not (folder in copoly_local_folders)]

if args.no_overwrite_chains:
    local_folders = os.listdir(os.path.abspath(args.dest_folder))
    copoly_local_folders = [folder for folder in local_folders if check_for_file(folder,'seq_analysis.json')]
    simfolders = [folder for folder in simfolders if not (folder in copoly_local_folders)]

dest_path=os.path.abspath(args.dest_folder)

for i,folder in enumerate(simfolders):
    print("Analyzing {} simulation of a total of {}".format(i+1,len(simfolders)))
    try:
        os.makedirs(os.path.join(dest_path,folder),exist_ok=True)  
    except FileExistsError:
        print("Folder already exists continuing without creating folder") 
    if args.local:
        result = csr.SimulationResults(os.path.join(os.path.abspath(args.folder),folder),os.path.join(os.path.abspath(dest_path),folder),is_remote=False)
    else:
        print("Args.folder is {}".format(args.folder))
        result = csr.SimulationResults(args.folder+'/'+folder,dest_path+'/'+folder)
    result.get_trajectory(os.path.join(dest_path,folder))
    if not args.only_chain:
        print(dest_path+'/'+folder)
        data = result.analyze_trajectory(dest_path+'/'+folder)
        data.to_csv(dest_path+'/'+folder+'/traj_analysis.csv')
        result.plot_trajectory(data)
    else:
        print("Sending results to folder {}".format(os.path.join(dest_path,folder)))
        data = result.analyze_trajectory_by_function(os.path.join(dest_path,folder))
        with open(os.path.join(dest_path,folder,'seq_analysis.json'),'w') as jfile:
            json.dump(data,jfile,sort_keys=True,indent=2)


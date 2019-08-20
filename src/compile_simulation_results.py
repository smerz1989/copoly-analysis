import numpy as np
import pandas as pd
import server_class as svc
import trajectory_class as trj
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from tqdm import tqdm, tqdm_notebook
from dump import dump
import time
import os


class SimulationResults(object):
    def __init__(self,sim_path,local_path):
        self.sim_path = sim_path
        self.server_connection = svc.ServerConnection()
        self.local_path = os.path.abspath(local_path)

    def is_trajectory_downloaded(self,local_path):
        file_path = os.path.abspath(local_path)
        traj_existing = (os.path.exists(local_path+'/atom_trj.lammpstrj'),
                         os.path.exists(local_path+'/bonddump.dump'),
                         os.path.exists(local_path+'/system.data'))
        return(traj_existing)
 

    def get_trajectory(self,local_path,redownload=False):
        traj_existing = self.is_trajectory_downloaded(local_path)
        if redownload or not traj_existing[0]:
            self.server_connection.get_file(self.sim_path+'/atom_trj.lammpstrj', local_path+'/atom_trj.lammpstrj')
        if redownload or not traj_existing[1]:    
            self.server_connection.get_file(self.sim_path+'/bonddump.dump',local_path+'/bonddump.dump')
        if redownload or not traj_existing[2]:
            self.server_connection.get_file(self.sim_path+'/system.data',local_path+'/system.data')
        self.bdump_path = local_path+'/bonddump.dump'
        self.atomtrj_path = local_path+'/atom_trj.lammpstrj'
        self.data_file_path = local_path+'/system.data'

    def analyze_trajectory(self,local_path):
        bdump = dump(local_path+'/bonddump.dump')
        snapshots = trj.construct_molecule_trajectory(local_path+'/system.data',bdump)
        timesteps = bdump.time()
        simulation_data = pd.DataFrame(columns=['NMonomers','fA','p','DOP','PDI','pAA','pBB','pAB'],
                                        index=timesteps)
        with tqdm(total=len(timesteps)) as pbar:
            for i,(timestep,snapshot) in enumerate(snapshots):
                simulation_data.loc[timestep,'NMonomers'] = snapshot.get_number_monomers()
                simulation_data.loc[timestep,'PDI'] = snapshot.get_pdi()
                simulation_data.loc[timestep,'DOP'] = snapshot.get_dop()
                (newpAA,newpBB,newpAB) = snapshot.get_all_probs() 
                simulation_data.loc[timestep,'pAA'] = newpAA
                simulation_data.loc[timestep,'pBB'] = newpBB
                simulation_data.loc[timestep,'pAB'] = newpAB            
                if i==0:
                    No = snapshot.get_number_monomers()
                    simulation_data.loc[timestep,'p']=0.
                else:
                    N = simulation_data.loc[timestep,'NMonomers']+snapshot.get_number_chains()
                    simulation_data.loc[timestep,'p'] = ((No-N)/No)
                simulation_data.loc[timestep,'fA']= snapshot.get_monomer_type_fraction()
                pbar.update(1)
        return(simulation_data)

 

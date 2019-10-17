import numpy as np
import pandas as pd
import server_class as svc
import trajectory_class as trj
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from tqdm import tqdm, tqdm_notebook
from dump import dump
import time
import os
import shutil
import re
import itertools as it
import dump_generator

class SimulationResults(object):
    def __init__(self,sim_path,local_path,is_remote=True):
        self.sim_path = sim_path
        self.is_remote = is_remote
        if is_remote:
            self.server_connection = svc.ServerConnection()
        self.local_path = os.path.abspath(local_path)

    def is_trajectory_downloaded(self,local_path):
        file_path = os.path.abspath(local_path)
        traj_existing = (os.path.exists(local_path+'/atom_trj.lammpstrj'),
                         os.path.exists(local_path+'/bonddump.dump'),
                         os.path.exists(local_path+'/system.data'))
        return(traj_existing)
 

    def get_trajectory(self,local_path,redownload=False):
        if self.is_remote:
            traj_existing = self.is_trajectory_downloaded(local_path)
            if redownload or not traj_existing[0]:
                print("Getting file {}".format(self.sim_path+'/atom_trj.lammpstrj'))
                self.server_connection.get_file(self.sim_path+'/atom_trj.lammpstrj', local_path+'/atom_trj.lammpstrj')
            if redownload or not traj_existing[1]:    
                self.server_connection.get_file(self.sim_path+'/bonddump.dump',local_path+'/bonddump.dump')
            if redownload or not traj_existing[2]:
                self.server_connection.get_file(self.sim_path+'/system.data',local_path+'/system.data')
        else:
            if not self.sim_path==local_path:
                print("Sim path is: {}\n\nLocal path is {}".format(self.sim_path,local_path))
                shutil.copy2(self.sim_path+'/atom_trj.lammpstrj',local_path)
                shutil.copy2(self.sim_path+'/bonddump.dump',local_path)
                shutil.copy2(self.sim_path+'/system.data',local_path)
        self.bdump_path = local_path+'/bonddump.dump'
        self.atomtrj_path = local_path+'/atom_trj.lammpstrj'
        self.data_file_path = local_path+'/system.data'

    def analyze_trajectory(self,local_path):
        print(r''+local_path+'/bonddump.dump')
        #bdump = dump(r''+local_path+'/bonddump.dump')
        num_timesteps = dump_generator.get_number_of_timesteps(local_path+'/atom_trj.lammpstrj')
        snapshots = trj.construct_molecule_trajectory_from_generators(local_path+'/system.data',local_path+'/bonddump.dump',local_path+'/atom_trj.lammpstrj')
        #timesteps = bdump.time()
        simulation_data = pd.DataFrame(columns=['NMonomers','fA','fB','p','DOP','PDI','pAA','pBB','pAB','pBA'],
                                        dtype=float)
        i=0
        with tqdm(total=num_timesteps) as pbar:
            for timestep,snapshot in enumerate(snapshots):
                simulation_data.loc[timestep,'NMonomers'] = snapshot.get_number_monomers()
                simulation_data.loc[timestep,'PDI'] = snapshot.get_pdi()
                simulation_data.loc[timestep,'DOP'] = snapshot.get_dop()
                (newpAA,newpBB,newpAB,newPBA) = snapshot.get_all_probs() 
                simulation_data.loc[timestep,'pAA'] = newpAA
                simulation_data.loc[timestep,'pBB'] = newpBB
                simulation_data.loc[timestep,'pAB'] = newpAB
                simulation_data.loc[timestep,'pBA'] = newpAB
                if i==0:
                    No = snapshot.get_number_monomers()
                    simulation_data.loc[timestep,'p']=0.
                else:
                    N = simulation_data.loc[timestep,'NMonomers']+snapshot.get_number_chains()
                    simulation_data.loc[timestep,'p'] = ((No-N)/No)
                simulation_data.loc[timestep,'fA']= snapshot.get_monomer_type_fraction()
                simulation_data.loc[timestep,'fB']=1-simulation_data.loc[timestep,'fA']
                i+=1
                pbar.update(1)
        return(simulation_data)


    def analyze_trajectory_by_function(self,local_path):
        print(local_path+'/bonddump.dump')
        bdump = dump(r''+local_path+'/bonddump.dump')
        snapshots = trj.construct_molecule_trajectory(local_path+'/system.data',bdump)
        timesteps = bdump.time()
        simulation_data = pd.DataFrame(columns=['block_lengths'],
                                        index=timesteps,dtype=float)
        sequence_data={}
        with tqdm(total=len(timesteps)) as pbar:
            for i,(timestep,snapshot) in enumerate(snapshots):
                sequences = [snapshot.get_chain_sequence_filtered(chain) for chain in snapshot.get_sequences()]
                conditional_probs = snapshot.get_conditional_probs()
                conditional_probs_dict = {'pAA': conditional_probs[0],
                                          'pBB': conditional_probs[1],
                                          'pAB': conditional_probs[2],
                                          'pBA': conditional_probs[3]}
                block_lengths = [get_block_lengths(chain_string) for chain_string in sequences]
                a_lengths_by_chain, b_lengths_by_chain = zip(*block_lengths)
                a_lengths = list(it.chain.from_iterable(a_lengths_by_chain))
                b_lengths = list(it.chain.from_iterable(b_lengths_by_chain))
                sequence_data[timestep]={'a_lengths': a_lengths, 'b_lengths': b_lengths,'sequences': sequences,'conditional_probs': conditional_probs_dict}
                #simulation_data.loc[timestep,'a_lengths'] = a_lengths
                #simulation_data.loc[timestep,'b_lengths'] = b_lengths
                pbar.update(1)
        return(sequence_data)


    def plot_trajectory(self,traj_df):
        sns.set_style("ticks")
        fig = plt.figure(figsize=(18,24))
        gs1 = gridspec.GridSpec(nrows=11,ncols=6,left=0.2,right=0.9,bottom=0.5,wspace=0.75,hspace=0.75)
        sns.set(style='ticks',palette='muted')
        sns.relplot(kind="line",data=traj_df['NMonomers'],ax=fig.add_subplot(gs1[0:3,0:3]))
        sns.relplot(kind="line",data=traj_df['DOP'],ax=fig.add_subplot(gs1[0:3,3:6]))
        sns.relplot(kind="line",data=traj_df['PDI'],ax=fig.add_subplot(gs1[3:6,0:3]))
        paxis = fig.add_subplot(gs1[3:6,3:6])
        paxis.set_ylim([0.,1.])
        sns.relplot(kind="line",data=traj_df['p'],ax=paxis)
        prob_names = ['pAA','pBB','pAB','pBA']
        prob_max = traj_df[prob_names].max().max()

        for i,prob in enumerate(prob_names):
            col = i*2
            if i<3:
                axis = fig.add_subplot(gs1[6:8,col:(col+2)])
                axis.set_ylim(0.,prob_max)
            sns.lineplot(data=traj_df[prob],ax=axis)

        mapaxis = fig.add_subplot(gs1[8:11,0:3])
        mapaxis.set_ylim([0.,1.])
        mbpaxis = fig.add_subplot(gs1[8:11,3:6])
        mbpaxis.set_ylim([0.,1.])
        sns.relplot(kind="line",data=traj_df['fA'],ax=mapaxis)
        sns.relplot(kind="line",data=traj_df['fB'],ax=mbpaxis)
        fig.savefig(self.local_path+'/monomer_plot.png')
        plt.close(fig)


def get_block_lengths(chain_string,a_string='3',b_string='4'):
    A_blocks = [len(block) for block in re.findall(r''+a_string+'+',chain_string)]
    B_blocks = [len(block) for block in re.findall(r''+b_string+'+',chain_string)]
    return(A_blocks,B_blocks)



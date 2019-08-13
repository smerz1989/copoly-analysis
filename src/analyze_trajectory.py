import trajectory_class as trj
import seaborn as sns
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm, tqdm_notebook
from dump import dump
import time

bdump = dump('bonddump.dump')
snapshots = trj.construct_molecule_trajectory('system.data',bdump)

number_of_monomers = []
dop = []
pdi = []
extentr = []
pAA = []
pBB = []
pAB = []
timesteps = []
monomer_fraction = []

time.sleep(1)
with tqdm(total=500,ascii=True) as pbar:
    time.sleep(1)
    for i,(timestep,snapshot) in enumerate(snapshots):
        #print("Analyzing timestep: "+str(timestep))
        timesteps.append(timestep)
        number_of_monomers.append(snapshot.get_number_monomers())
        dop.append(snapshot.get_dop())
        pdi.append(snapshot.get_pdi())
        (newpAA,newpBB,newpAB) = snapshot.get_all_probs() 
        pAA.append(newpAA)
        pBB.append(newpBB)
        pAB.append(newpAB)
        if i==0:
            No = snapshot.get_number_monomers()
            extentr.append(0.)
        else:
            N = snapshot.get_number_monomers()+snapshot.get_number_chains()
            extentr.append((No-N)/No)
        monomer_fraction.append(snapshot.get_monomer_type_fraction())
        pbar.update(1)


mon_df = pd.DataFrame(list(zip(timesteps,number_of_monomers,dop,pdi,extentr,pAA,pBB,pAB,monomer_fraction)),columns=['timestep','nmonomers','dop','pdi','p','pAA','pBB','pAB','mfraction'])


sns.set_style("ticks")
fig = plt.figure(figsize=(18,24))
gs1 = gridspec.GridSpec(nrows=8,ncols=6,left=0.2,right=0.9,bottom=0.5,wspace=0.75,hspace=0.75)

#fig.set_size_inches(24.5,18)

sns.set(style='ticks',palette='muted')
sns.relplot(x='timestep',y='nmonomers',kind="line",data=mon_df,ax=fig.add_subplot(gs1[0:3,0:3]))
sns.relplot(x='timestep',y='dop',kind="line",data=mon_df,ax=fig.add_subplot(gs1[0:3,3:6]))
sns.relplot(x='timestep',y='pdi',kind="line",data=mon_df,ax=fig.add_subplot(gs1[3:6,0:3]))
paxis = fig.add_subplot(gs1[3:6,3:6])
paxis.set_ylim([0.,1.])
sns.relplot(x='timestep',y='p',kind="line",data=mon_df,ax=paxis)

prob_names = ['pAA','pBB','pAB']

for i,prob in enumerate(prob_names):
    col = i*2
    axis = fig.add_subplot(gs1[6:8,col:(col+2)])
    axis.set_ylim(0.,0.5)
    sns.relplot(x='timestep',y=prob,kind="line",data=mon_df,ax=axis)

#mapaxis = fig.add_subplot(gs1[6:9,0:3])
#mapaxis.set_ylim([0.,1.])
#mbpaxis = fig.add_subplot(gs1[6:9,3:6])
#mbpaxis.set_ylim([0.,1.])
#sns.relplot(x='timestep',y='mfraction',kind="line",data=mon_df,ax=mapaxis)
#sns.relplot(x='timestep',y=1-'mfraction',kind="line",data=mon_df,ax=mbpaxis)


fig.savefig('monomer_plot.png')


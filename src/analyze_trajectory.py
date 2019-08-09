import trajectory_class as trj
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

snapshots = trj.construct_molecule_trajectory('system.data','bonddump.dump')

number_of_monomers = np.empty(len(snapshots))
dop = np.empty(len(snapshots))
pdi = np.empty(len(snapshots))
extentr = np.empty(len(snapshots))
pAA = np.empty(len(snapshots))
pBB = np.empty(len(snapshots))
pAB = np.empty(len(snapshots))

#length_dict = {timestep: np.array(snapshots[timestep].get_chain_lengths()) for timestep in snapshots}

for i,timestep in enumerate(snapshots):
    print("Analyzing timestep: "+str(timestep))
    number_of_monomers[i] = snapshots[timestep].get_number_monomers()
    dop[i] = snapshots[timestep].get_dop()
    pdi[i] = snapshots[timestep].get_pdi()
    pAA[i], pBB[i], pAB[i] = snapshots[timestep].get_all_probs() 
    if i==0:
        No = snapshots[timestep].get_number_monomers()
        extentr[i] = 0
    else:
        N = snapshots[timestep].get_number_monomers()+snapshots[timestep].get_number_chains()
        extentr[i] = (No-N)/No
    del snapshots[timestep]


mon_df = pd.DataFrame(list(zip(snapshots.keys(),number_of_monomers,dop,pdi,extentr,pAA,pBB,pAB)),columns=['timestep','nmonomers','dop','pdi','p','pAA','pBB','pAB'])
#mon_df['timestep'] = mon_df.index 

fig, axises = plt.subplots(3,3)
#fig.tight_layout()
fig.subplots_adjust(left=0.15,right=0.9,bottom=0.15,wspace=0.5,hspace=0.5)
fig.set_size_inches(12.5,6.5)

sns.set(style='ticks',palette='muted')
axises[0][1].axis('off')
axises[1][1].axis('off')
sns.relplot(x='timestep',y='nmonomers',data=mon_df,ax=axises[0][0])
sns.relplot(x='timestep',y='dop',data=mon_df,ax=axises[0][2])
sns.relplot(x='timestep',y='pdi',data=mon_df,ax=axises[1][0])
axises[1][2].set_ylim(0.,1.)
sns.relplot(x='timestep',y='p',data=mon_df,ax=axises[1][2])
for axis in axises[2]:
    axis.set_ylim(0.,0.5)
sns.relplot(x='timestep',y='pAA',data=mon_df,ax=axises[2][0])
sns.relplot(x='timestep',y='pBB',data=mon_df,ax=axises[2][1])
sns.relplot(x='timestep',y='pAB',data=mon_df,ax=axises[2][2])

fig.savefig('monomer_plot.png')


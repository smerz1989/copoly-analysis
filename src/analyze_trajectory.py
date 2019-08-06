import trajectory_class as trj
import seaborn as sns
import numpy as np
import pandas as pd

snapshots = trj.construct_molecule_trajectory('system.data','bonddump.dump')

number_of_monomers = np.empty(len(snapshots))
dop = np.empty(len(snapshots))
length_dict = {timestep: np.array(snapshots[timestep].get_chain_lengths()) for timestep in snapshots}
for i,timestep in enumerate(snapshots):
    number_of_monomers[i] = len(np.extract(length_dict[timestep]==3,length_dict[timestep]))
    dop[i] = np.mean([length for length in length_dict[timestep] if length>1])
#    mave = mean([length**2 for length in length_dict[timestep] if length>1])


mon_df = pd.DataFrame(list(zip(snapshots.keys(),number_of_monomers,dop)),columns=['timestep','nmonomers','dop'])
#mon_df['timestep'] = mon_df.index 

sns.set(style='ticks',palette='muted')
g = sns.relplot(x='timestep',y='dop',data=mon_df)
g.fig.savefig('monomer_plot.png')


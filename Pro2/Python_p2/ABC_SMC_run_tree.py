import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from ABC_SMC_DVmodel import SMC_ABC, single_trait_sim
from DV_model_sim_along_phy import DVtraitsim_tree
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])
# Observation generated
#  load the data for a given tree
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
simresult = DVtraitsim_tree(file = file, gamma1 = par_obs[0], a = par_obs[1],scalor = 1000)
evo_time, total_species = simresult[0].shape
evo_time = evo_time-1
trait_RI_dr = simresult[0]
population_RI_dr = simresult[1]
# empirical data for trait and population
trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]
# observation data
obs = np.array([trait_dr_tips,population_tips])


epsilon = 30
timestep = 10
particlesize = 1000
prior = [0.03,0.5,0.2,0.5]
SMC_ABC_result = SMC_ABC(timestep = timestep, particlesize = particlesize, obs = obs, prior = prior, epsilon = epsilon,file = file)
pickle.dump( SMC_ABC_result, open( "c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result_DV.p", "wb" ) )

# read result
# SMC_ABC_result = pickle.load( open( "c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result_DV.p", "rb" ) )

# plot parameters distribution for each time step
gamma = SMC_ABC_result['gamma'].flatten()
a = SMC_ABC_result['a'].flatten()
par = np.array([gamma,a]).flatten()
# column and row labels
t_str = ['t=%i' % x for x in range(1,timestep+1)]
par_str = ['$\gamma$','a']

tlabel = np.repeat(t_str, particlesize)
tlabel = np.tile(tlabel,2)
parlabel = np.repeat(['gamma','a'],timestep*particlesize)

df_ori = pd.DataFrame(dict(par=par,tlabel=tlabel,parlabel = parlabel))

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
g_ori = sns.FacetGrid(df_ori, col = "parlabel"  ,row="tlabel", hue="tlabel", aspect=15, size=.5,
                      palette=pal,
                      xlim=[-1,5])

# Draw the densities in a few steps
g_ori.map(sns.kdeplot, "par", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g_ori.map(sns.kdeplot, "par", clip_on=False, color="w", lw=2, bw=.2)

g_ori.map(plt.axhline, y=0, lw=2, clip_on=False)

# column title
for i in range(len(par_str)):
    g_ori.axes[0,i].text(0.5, 0.8, s = par_str[i], fontweight="bold", color = pal[0],
                          ha="right", va="center", transform=g_ori.axes[0,i].transAxes)
# row legend
for i in range(len(t_str)):
    g_ori.axes[i, 1].text(0.9, .2, s = t_str[i], fontweight="bold", color = pal[i],
                          ha="right", va="center", transform=g_ori.axes[i, 1].transAxes)
# Set the subplots to overlap
g_ori.fig.subplots_adjust(hspace=-.25)

# Remove axes details that don't play will with overlap
g_ori.set_titles("")
g_ori.set(yticks=[])
g_ori.despine(bottom=True, left=True)
g_ori.savefig("c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result_DV.png")
# g_ori.set_xlabels("Distance")

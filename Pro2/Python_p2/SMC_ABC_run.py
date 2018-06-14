import sys, os
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from ABC_MCMC import SMC_ABC, single_trait_sim
import numpy as np
from Trait_sim_in_branches_stat import traitsim
import scipy.stats
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=0)

epsilon = 30
timestep = 8
particlesize = 100
prior = [0.3,0.5,0.2,0.5]
SMC_ABC_result = SMC_ABC(timestep = timestep, particlesize = particlesize, obs = obs, prior = prior, epsilon = epsilon)
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result.npy",SMC_ABC_result)
pickle.dump( SMC_ABC_result, open( "c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result.p", "wb" ) )


SMC_ABC_result = pickle.load( open( "c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result.p", "rb" ) )

gamma = SMC_ABC_result['gamma'].flatten()
a = SMC_ABC_result['a'].flatten()
par = np.array([gamma,a]).flatten()
tlabel = np.repeat(['t=1','t=2','t=3','t=4','t=5','t=6','t=7','t=8'], particlesize)
tlabel = np.tile(tlabel,2)
parlabel = np.repeat(['gamma','a'],timestep*particlesize)

df_ori = pd.DataFrame(dict(par=par,tlabel=tlabel,parlabel = parlabel))

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#
# def kde_plot(x, color_pal, **kwargs):
#     cmap = sns.cubehelix_palette(color_pal, rot=-.25, light=.7)
#     sns.kdeplot(x, color_pal=cmap, **kwargs)


g_ori = sns.FacetGrid(df_ori, col = "parlabel"  ,row="tlabel", hue="tlabel", aspect=15, size=.5,
                      palette=pal,
                      xlim=[0,5])

# Draw the densities in a few steps
g_ori.map(sns.kdeplot, "par", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g_ori.map(sns.kdeplot, "par", clip_on=False, color="w", lw=2, bw=.2)

g_ori.map(plt.axhline, y=0, lw=2, clip_on=False)

# # Define and use a simple function to label the plot in axes coordinates
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color=color,
#             ha="left", va="center", transform=ax.transAxes)
# #
# # g_ori.axes[0,5].text(0, .2, label, fontweight="bold", color=pal,
# #             ha="left", va="center", transform=g_ori.axes[0,5].transAxes)
# g_ori.map(label, "dis_gamma_sort")

t_str = ['t=1','t=2','t=3','t=4','t=5','t=6','t=7','t=8']
par_str = ['$\gamma$','a']


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
# g_ori.set_xlabels("Distance")

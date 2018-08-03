import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from testABC_sim2 import DVtraitsim_tree
from testABCfunc import SMC_ABC_MS
import timeit

tic_test = timeit.default_timer()

# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.5])

# Observation generated
# Load the data for a given tree
# The directory of the tree data
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
# Simulate data; model = 0: normal distribution; model=1: uniform distribution
simresult = DVtraitsim_tree(file = file,replicate = 3,model=0, gamma1 = par_obs[0], a = par_obs[1],scalar = 1000)
# We only need the data at tips.
evo_time, total_species = simresult[0].shape
# evo_time = 1
# total_species = len(simresult[0])
evo_time = evo_time-1
trait_RI_dr = simresult[0]
population_RI_dr = simresult[1]
traitvar = simresult[3]
trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]
traitvar = traitvar[evo_time,:][~np.isnan(traitvar[evo_time,:])]
# observation data
obs = np.array([trait_dr_tips,population_tips,traitvar])

#SMC test for model selection
timestep= 30
particlesize = 1000
epsilon = 100
prior = [0,1,0,1]
testsmc = SMC_ABC_MS(timestep, particlesize, obs, epsilon, prior, file, sort = 0)




# possessing data
total = np.arange(timestep)
total.fill(1)
total = pd.Series(total)
modeldata = np.sum(testsmc['model'],axis=1)/particlesize
modeldata = pd.Series(modeldata)
population = np.arange(timestep)
pop_strs = ["Pop %d" % (d+1) for d in range(timestep)]
population_str = pd.Series(pop_strs)

smcdata = {'total':total,'modeldata':modeldata,'population':population_str}
smcdf = pd.DataFrame(data=smcdata)

#plot distribution for model selection
# Initialize the matplotlib figure
sns.set(style="whitegrid")
f, ax = plt.subplots(figsize=(6, 15))

# Plot the total crashes
sns.set_color_codes("pastel")
sns.barplot(x="total", y="population", data=smcdf,
            label="M0", color="b")

# Plot the crashes where alcohol was involved
sns.set_color_codes("muted")
sns.barplot(x="modeldata", y="population", data=smcdf,
            label="M1", color="b")

# Add a legend and informative axis label
ax.legend(ncol=2, loc="lower right", frameon=True)
ax.set(xlim=(0, 1), ylabel="Population",
       xlabel="Model selection")
sns.despine(left=True, bottom=True)




gamma = testsmc['gamma'].flatten()
a = testsmc['a'].flatten()
par = np.array([gamma,a]).flatten()
tlabel = np.repeat(pop_strs, particlesize)
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
                      xlim=[0,1])

# Draw the densities in a few steps
g_ori.map(sns.kdeplot, "par", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g_ori.map(sns.kdeplot, "par", clip_on=False, color="w", lw=1, bw=.2)

g_ori.map(plt.axhline, y=0, lw=1, clip_on=False)

# # Define and use a simple function to label the plot in axes coordinates
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color=color,
#             ha="left", va="center", transform=ax.transAxes)
# #
# # g_ori.axes[0,5].text(0, .2, label, fontweight="bold", color=pal,
# #             ha="left", va="center", transform=g_ori.axes[0,5].transAxes)
# g_ori.map(label, "dis_gamma_sort")

t_str = pop_strs
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
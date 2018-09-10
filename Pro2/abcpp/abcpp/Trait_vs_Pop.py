import os
import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kendalltau
sns.set(style="white")
#
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec

dir_path = 'c:/Liang/Code/Pro2/abcpp'
files = dir_path + '/tree_data/example1/'

td = DVTreeData(path=files, scalar=10000)
K = 10e8
nu=1/(100*K)
num = 10
trait_w = []
trait_v = []
pop_w = []

for gamma in gamma_vec:
    for a in a_vec:
        # let's try to find a true simulation:
        obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                            split_stddev=0.2)
        trait_data = ()
        population_data = ()
        traitvar_data = ()
        for loop in range(1,num):
            str = 'gamma = %.3f; a = %.3f; loop = %d' % (gamma,a,loop)
            print(str)
            par_obs = np.array([gamma, a])
            simresult = dvcpp.DVSim(td, obs_param)
            if simresult['sim_time'] == td.evo_time:
                trait_tips = simresult['Z']
                population_tips = simresult['N']
                traitvar_tips = simresult['V']
                # empirical data for trait and population
                trait_tips = trait_tips[trait_tips != np.nan]
                population_tips = population_tips[population_tips != np.nan]
                traitvar_tips = traitvar_tips[traitvar_tips != np.nan]

                trait_data = np.append(trait_data, trait_tips)
                population_data = np.append(population_data, population_tips)
                traitvar_data = np.append(traitvar_data, traitvar_tips)
            if simresult['sim_time'] < td.evo_time:
                print('Jump to the next loop')

        trait_w.append(trait_data)
        trait_v.append(traitvar_data)
        pop_w.append(population_data)

num_tips = len(trait_tips)

normed_trait = []
normed_traitvar = []
normed_pop = []

for i in range(0,36):
    if len(trait_w[i])==0:
        normed_trait.append([0])
        normed_traitvar.append([0])
        normed_pop.append([0])
    else:
        normed_trait.append((trait_w[i]- np.min(trait_w[i])) / (np.max(trait_w[i]) - np.min(trait_w[i])))
        normed_traitvar.append((trait_v[i]- np.min(trait_v[i])) / (np.max(trait_v[i]) - np.min(trait_v[i])))
        normed_pop.append((pop_w[i] - np.min(pop_w[i])) / (np.max(pop_w[i]) - np.min(pop_w[i])))


count = 0
label_a = (['a=0','a=.001','a=.01','a=.01','a=.1','a=.5'])
label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5'])

# Set up the matplotlib figure
f, axes = plt.subplots(6, 6, figsize=(9, 9), sharex=True, sharey=True) #, sharex=True, sharey=True

# Rotate the starting point around the cubehelix hue circle
for ax, s in zip(axes.flat, np.linspace(0, 3,36)):
    # Create a cubehelix colormap to use with kdeplot
    cmap = sns.cubehelix_palette(start=s, light=1, as_cmap=True)
    if len(normed_trait[count]) == 1:
        ax.plot()

    else:
        # trait = trait_w[count]
        # traitvar = trait_v[count]
        # pop = pop_w[count]
        trait = normed_trait[count]
        traitvar = normed_traitvar[count]
        pop = normed_pop[count]
        ax.set_xlim([0,1])
        # Generate and plot a random bivariate dataset
        # sns.kdeplot(trait, pop, cmap=cmap, shade=True, cut=5, ax=ax)
        sns.scatterplot(trait, pop,ax=ax)

    if count in range(0,6):
        ax.title.set_text(label_a[count])

    if count in ([5, 11, 17, 23, 29, 35]):
        ax.set_ylabel(label_gamma[int(count/6)])
        ax.yaxis.set_label_position("right")
    count += 1
f.text(0.5, 0, 'Trait', ha='center')
f.text(0.01, 0.5, 'Population', va='center', rotation='vertical')
f.tight_layout()

# f.savefig('C:/Liang/Googlebox/Research/Project2/DVmodel/1stClusterStudy/traitvsvar.png')

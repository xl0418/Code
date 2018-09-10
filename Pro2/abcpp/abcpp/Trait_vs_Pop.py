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
sns.set(style="ticks")
#
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec

dir_path = 'c:/Liang/Code/Pro2/abcpp'
files = dir_path + '/tree_data/example1/'

td = DVTreeData(path=files, scalar=10000)
K = 10e8
nu=1/(100*K)
num = 6
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
    # trait = trait_w[count]
    # traitvar = trait_v[count]
    # pop = pop_w[count]
    else:
        trait = normed_trait[count]
        traitvar = normed_traitvar[count]
        pop = normed_pop[count]
        ax.set_xlim([0,1])
        # Generate and plot a random bivariate dataset
        sns.kdeplot(trait, pop, cmap=cmap, shade=True, cut=5, ax=ax)
    if count in range(0,6):
        ax.title.set_text(label_a[count])

    if count in ([5, 11, 17, 23, 29, 35]):
        ax.set_ylabel(label_gamma[int(count/6)])
        ax.yaxis.set_label_position("right")
    count += 1
f.text(0.5, 0, 'Trait', ha='center')
f.text(0.01, 0.5, 'Population', va='center', rotation='vertical')
f.tight_layout()

f.savefig('C:/Liang/Googlebox/Research/Project2/DVmodel/1stClusterStudy/traitvsvar.png')


sns.jointplot(trait_w[5], trait_v[5], kind="kde", stat_func=None, color="#4CB391")


#
# gamma_label = np.repeat(['gamma=.001','gamma=.01'], (num-1)*num_tips*len(a_vec)) #,'$\gamma$=.01','$\gamma$=.1','$\gamma$=.5'
# label_a = np.repeat(['a=.001','a=.01'], (num-1)*num_tips) #,'a=.01','a=.1','a=.5'
# a_label =np.tile(label_a,len(gamma_vec))
# df_ori = pd.DataFrame(dict(trait=trait_data,pop=population_data,traitvar=traitvar_data,
#                            a_label=a_label,gamma_label = gamma_label))
#
# x = [df_ori[(df_ori['gamma_label']== 'gamma=.001') & (df_ori['a_label']== 'a=.001')]['trait'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.01') & (df_ori['a_label']== 'a=.001')]['trait'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.001') & (df_ori['a_label']== 'a=.01')]['trait'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.01') & (df_ori['a_label']== 'a=.01')]['trait'].values]
# y = [df_ori[(df_ori['gamma_label']== 'gamma=.001') & (df_ori['a_label']== 'a=.001')]['pop'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.01') & (df_ori['a_label']== 'a=.001')]['pop'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.001') & (df_ori['a_label']== 'a=.01')]['pop'].values,
# df_ori[(df_ori['gamma_label']== 'gamma=.01') & (df_ori['a_label']== 'a=.01')]['pop'].values]
# np.max(x)
# np.min(x)


tvmap = sns.FacetGrid(df_ori, col = "gamma_label"  ,row="label", hue="label", aspect=15, size=.5,
                      palette=pal)

# Draw the densities in a few steps
tvmap.map(sns.kdeplot, "dis_gamma_sort", clip_on=False, shade=True, alpha=1, lw=1, bw=2)


gamma_str = ['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5']
# gamma_str = ['$\gamma$=0','$\gamma$=.001']

for i in range(len(gamma_str)):
    tvmap.axes[0,i].text(1, 0.8, s = gamma_str[i], fontweight="bold", color = pal[0],
                          ha="right", va="center", transform=tvmap.axes[0,i].transAxes)
a_str = ['a=0','a=.001','a=.01','a=.1','a=.5']
for i in range(len(a_str)):
    tvmap.axes[i, 4].text(1, .2, s = a_str[i], fontweight="bold", color = pal[i],
                          ha="right", va="center", transform=tvmap.axes[i, 4].transAxes)
# Set the subplots to overlap
tvmap.fig.subplots_adjust(hspace=-.25)

# Remove axes details that don't play will with overlap
tvmap.set_titles("")
tvmap.set(yticks=[])
tvmap.despine(bottom=True, left=True)
# g_ori.set_xlabels("Distance")
import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from ABC_SMC_DVmodel import SMC_ABC, single_trait_sim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kendalltau
sns.set(style="ticks")
#
gamma_vec = np.array([0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec
# Observation parameters [gamma,a]
# trait_data = ()
# population_data = ()
# traitvar_data = ()
# Observation generated
#  load the data for a given tree
# file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
file = 'D:\\Googlebox\\Python\\Project2\\R-tree_sim\\'

num = 6
trait_w = []
trait_v = []
pop_w = []

for gamma in gamma_vec:
    for a in a_vec:
        trait_data = ()
        population_data = ()
        traitvar_data = ()
        for loop in range(1,num):
            str = 'gamma = %.3f; a = %.3f; loop = %d' % (gamma,a,loop)
            print(str)
            par_obs = np.array([gamma, a])
            simresult = single_trait_sim(file = file,replicate = 0, par = par_obs)
            trait_tips = simresult[0]
            population_tips = simresult[1]
            traitvar_tips = simresult[3]
            # empirical data for trait and population
            trait_tips = trait_tips[trait_tips != 0]
            population_tips =  population_tips[population_tips != 0]
            traitvar_tips = traitvar_tips[traitvar_tips != 0]

            trait_data = np.append(trait_data,trait_tips)
            population_data = np.append(population_data,population_tips)
            traitvar_data = np.append(traitvar_data,traitvar_tips)
        trait_w.append(trait_data)
        trait_v.append(traitvar_data)
        pop_w.append(population_data)

num_tips = len(trait_tips)

normed_trait = []
normed_traitvar = []
normed_pop = []

for i in range(0,25):
    normed_trait.append((trait_w[i]- np.min(trait_w[i])) / (np.max(trait_w[i]) - np.min(trait_w[i])))
    normed_traitvar.append((trait_v[i]- np.min(trait_v[i])) / (np.max(trait_v[i]) - np.min(trait_v[i])))
    normed_pop.append((pop_w[i] - np.min(pop_w[i])) / (np.max(pop_w[i]) - np.min(pop_w[i])))

np.max(normed_pop)
np.min(normed_pop)
count = 0
label_a = (['a=.001','a=.01','a=.01','a=.1','a=.5'])
label_gamma = (['$\gamma$=.001','$\gamma$=.01','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5'])

# Set up the matplotlib figure
f, axes = plt.subplots(5, 5, figsize=(9, 9)) #, sharex=True, sharey=True

# Rotate the starting point around the cubehelix hue circle
for ax, s in zip(axes.flat, np.linspace(0, 3, 25)):
    # Create a cubehelix colormap to use with kdeplot
    cmap = sns.cubehelix_palette(start=s, light=1, as_cmap=True)
    # trait = trait_w[count]
    # traitvar = trait_v[count]
    # pop = pop_w[count]
    trait = normed_trait[count]
    traitvar = normed_traitvar[count]
    pop = normed_pop[count]
    ax.set_xlim([0,1])
    # Generate and plot a random bivariate dataset
    sns.kdeplot(trait, traitvar, cmap=cmap, shade=True, cut=5, ax=ax)
    if count in range(0,5):
        ax.title.set_text(label_a[count])

    if count in ([4,9,14,19,24]):
        ax.set_ylabel(label_gamma[int(count/5)])
        ax.yaxis.set_label_position("right")
    count += 1
f.text(0.5, 0, 'Trait', ha='center')
f.text(0.02, 0.5, 'Trait variance', va='center', rotation='vertical')
f.tight_layout()

f.savefig('C:/Liang/Googlebox/Research/Project2/DVmodel/1stClusterStudy/traitvsvar.png')


sns.jointplot(x, y, kind="kde", stat_func=None, color="#4CB391")


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
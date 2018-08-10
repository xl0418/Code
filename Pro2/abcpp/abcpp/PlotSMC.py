import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


generation = len(para_data['gamma'])
population = len(para_data['gamma'][0])

pop_strs = ["G %d" % (d+1) for d in range(generation)]

gamma = np.concatenate(para_data['gamma']).ravel()
a = np.concatenate(para_data['a']).ravel()
par = np.array([gamma,a]).flatten()
tlabel = np.repeat(pop_strs, population)
tlabel = np.tile(tlabel,2)
parlabel = np.repeat(['gamma','a'],generation*population)

df_ori = pd.DataFrame(dict(par=par,tlabel=tlabel,parlabel = parlabel))

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#
# def kde_plot(x, color_pal, **kwargs):
#     cmap = sns.cubehelix_palette(color_pal, rot=-.25, light=.7)
#     sns.kdeplot(x, color_pal=cmap, **kwargs)


g_ori = sns.FacetGrid(df_ori, col = "parlabel"  ,row="tlabel", hue="tlabel", aspect=15, height=.5,
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



import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="dark")
p=2
gamma_vec = [0,0.001,0.01,0.1,0.5,1]
a_vec =  [0,0.001,0.01,0.1,0.5,1]
gamma_list = []
a_list = []
fit_list = []
# Create the data
if platform.system()=='Windows':
    filedir = 'C:/Liang/Googlebox/Research/Project2/smc_newdata/test3/'
elif platform.system()=='Darwin':
    filedir = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/abcpp/smcndata/tree9+1q/'
for gamma_index in range(len(gamma_vec)):
    for a_index in range(len(a_vec)):
        file = filedir + 'smc%dg%da.npy' % (gamma_index, a_index)
        if os.path.isfile(file):
            para_data = np.load(file).item()
            generation = len(para_data['gamma'])
            population = len(para_data['gamma'][0])
            gamma_list.append(para_data['gamma'][generation - 1])
            a_list.append(para_data['a'][generation - 1])
            fit_list.append(para_data['fitness'][generation - 1])

        else:
            gamma_list.append([0])
            a_list.append([0])
            fit_list.append([0])


label_a = (['a=0','a=.001','a=.01','a=.1','a=.5','a=1'])
label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5','$\gamma$=1'])
row_a = len(label_a)
row_gamma = len(label_gamma)
pos=[1,2,3,4]
pos1 = [1,2]
pos2 = [3,4]
pos_label=['$\gamma$','a','$\gamma$','a']
# Set up the matplotlib figure
f, axes = plt.subplots(row_a, row_gamma, figsize=(9, 9), sharex=True, sharey=True) #
gamma_vec_point = np.repeat(gamma_vec,len(a_vec))
a_vec_point = np.tile(a_vec,len(gamma_vec))
cmap = sns.cubehelix_palette(p, rot=-.5, dark=.3)

count = 0
# Rotate the starting point around the cubehelix hue circle
for ax in axes.flat:
    gamma = gamma_list[count]
    a = a_list[count]
    fitness = fit_list[count]
    # d=np.column_stack((gamma,a))
    q5 = np.argsort(fitness)[-population // 20] # best 5%
    fit_index = np.where(fitness > fitness[q5])[0]
    d=[gamma,a]
    gamma5th=gamma_list[count][fit_index]
    a5th = a_list[count][fit_index]
    d5th=[gamma5th,a5th]
    dg=[gamma,a,gamma5th,a5th]
    if len(gamma)==1:
        ax.text(0.45,0.45,"X")

    # Create a cubehelix colormap to use with kdeplot
    else:
        # Generate and plot a random bivariate dataset
        # sns.kdeplot(a, gamma, cmap=cmap, shade=True, cut=5, ax=ax)
        # sns.violinplot(data=d, palette=cmap, inner="points",ax=ax)
        ax.violinplot(dataset=dg, positions=pos, points=20, widths=0.8,
                    showmeans=True, showextrema=True, showmedians=False)
        # ax.violinplot(dataset=d5th, positions=pos2, points=20, widths=0.3,
        #               showmeans=True, showextrema=True, showmedians=False)
        ax.scatter(x=1,y=gamma_vec_point[count],color='r',s=5,alpha=1)
        ax.scatter(x=2,y=a_vec_point[count],color='r',s=5,alpha=1)
        ax.scatter(x=3,y=gamma_vec_point[count],color='g',s=5,alpha=1)
        ax.scatter(x=4,y=a_vec_point[count],color='g',s=5,alpha=1)
    if count in range(0,row_a):
        ax.title.set_text(label_a[count])

    if count in ([5,11,17,23,29,35]):
        ax.set_ylabel(label_gamma[int(count/row_a)])
        ax.yaxis.set_label_position("right")
    # ax.yaxis.set_major_locator(plt.NullLocator())
    # ax.xaxis.set_major_locator(plt.NullLocator())
    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.set(ylim=(-.2,1.5))
    ax.set_xticklabels(pos_label)
    count += 1

f.text(0.5, 0, '', ha='center')
f.text(0.01, 0.5, '', va='center', rotation='vertical')
f.tight_layout()
f.show()
# f.savefig('C:/Liang/Code/Pro2/abcpp/abcpp/smcdata/Tree2plot.png')

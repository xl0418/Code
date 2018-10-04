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
nv_list = []
truenv=1e-11
lowylim=truenv*(-0.4)
upperylim=2*truenv*6/5
# Create the data
if platform.system()=='Windows':
    filedir = 'C:/Liang/Googlebox/Research/Project2/smcvdata/tree5+nvka/'
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
            nv_list.append(para_data['nu'][generation - 1])

        else:
            gamma_list.append([0])
            a_list.append([0])
            nv_list.append([0])

label_a = (['a=0','a=.001','a=.01','a=.1','a=.5','a=1'])
label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5','$\gamma$=1'])
row_a = len(label_a)
row_gamma = len(label_gamma)
pos = [0,1]
posnv = [2]
pos_label=['$\gamma$','a','$\\nu$']
# Set up the matplotlib figure
f, axes = plt.subplots(row_a, row_gamma, figsize=(9, 9)) #, sharex=True, sharey=True
gamma_vec_point = np.repeat(gamma_vec,len(a_vec))
a_vec_point = np.tile(a_vec,len(gamma_vec))
cmap = sns.cubehelix_palette(p, rot=-.5, dark=.3)

count = 0
# Rotate the starting point around the cubehelix hue circle
for ax in axes.flat:
    gamma = gamma_list[count]
    a = a_list[count]
    nv = nv_list[count]
    # d=np.column_stack((gamma,a))
    d=[gamma,a]
    axnv = ax.twinx()  # instantiate a second axes that shares the same x-axis
    if len(gamma)==1:
        axnv.plot()
        if count in ([5,11,17,23,29,35]):
            axnv.set_ylabel(label_gamma[int(count/row_a)])
    # Create a cubehelix colormap to use with kdeplot
    else:
        # Generate and plot a random bivariate dataset
        ax.violinplot(dataset=d, positions=pos, points=20, widths=0.3,
                    showmeans=True, showextrema=True, showmedians=False)
        ax.scatter(x=0,y=gamma_vec_point[count],color='r',s=10,alpha=1)
        ax.scatter(x=1,y=a_vec_point[count],color='r',s=10,alpha=1)
        axnv.violinplot(dataset=nv, positions=posnv, points=20, widths=0.3,
                    showmeans=True, showextrema=True, showmedians=False)
        # axnv.get_yaxis().set_ticks([])
        # ax.get_yaxis().set_ticks([])
        axnv.scatter(x=2,y=truenv,color='r',s=10,alpha=1)
        if count in ([5, 11, 17, 23, 29, 35]):
            axnv.set_ylabel(label_gamma[int(count / row_a)])
    if count in range(0,row_a):
        ax.title.set_text(label_a[count])
        # ax.yaxis.set_label_position("right")
        # ax.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    if count in ([0,1,2,3,4,6,7,8,9,10,12,13,14,15,16,18,19,20,21,22,24,25,26,27,28,30,31,32,33,34]):
        axnv.get_yaxis().set_ticks([])
    if count in ([1, 2, 3, 4,5, 7, 8, 9, 10,11, 13, 14, 15, 16,17, 19, 20, 21, 22,23, 25, 26, 27, 28,29, 31, 32, 33, 34,35]):
        ax.get_yaxis().set_ticks([])
    if count in range(30):
        ax.get_xaxis().set_ticks([])
    if count in range(30,36):
        plt.xticks([0, 1, 2], pos_label)

    #  ax.xaxis.set_major_locator(plt.NullLocator())
    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.set(ylim=(-.2,1.2))
    axnv.set_ylim(lowylim, upperylim)
    count += 1

f.text(0.5, 0, '', ha='center')
f.text(0.01, 0.5, '', va='center', rotation='vertical')
# f.tight_layout()
plt.show(f)

# f.savefig('C:/Liang/Code/Pro2/abcpp/abcpp/smcdata/Tree2plot.png')
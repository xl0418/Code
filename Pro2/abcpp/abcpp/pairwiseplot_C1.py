import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="white")
treeno_vec = [9,11]
gno_vec = [0,0,0,0,1,1,1,2,2,3]
ano_vec = [2,3,4,5,3,4,5,4,5,5]
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec =  np.array([0,0.001,0.01,0.1,0.5,1])
gamma_list = []
a_list = []
nv_list = []
fit_list = []

truenv=1e-11
lowylim=truenv*(-0.4)
upperylim=2*truenv*6/5
for indicator in range(10):
    for tree in treeno_vec:
        g=gno_vec[indicator]
        a=ano_vec[indicator]
        if platform.system() == 'Windows':
            file = 'C:/Liang/Googlebox/Research/Project2/smcvdata/tree%d+nv/smc%dg%da.npy' % (tree,g,a)
        elif platform.system() == 'Darwin':
            file = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/abcpp/smcndata/tree%d+nv/smc%dg%da.npy' % (tree,g,a)
        if os.path.isfile(file):
            para_data = np.load(file).item()
            generation = len(para_data['gamma'])
            population = len(para_data['gamma'][0])
            gamma_list.append(para_data['gamma'][generation - 1])
            a_list.append(para_data['a'][generation - 1])
            nv_list.append(para_data['nu'][generation - 1])
            fit_list.append(para_data['fitness'][generation - 1])
        else:
            gamma_list.append([0])
            a_list.append([0])
            nv_list.append([0])
            fit_list.append([0])
label_tree = (['Test %d' % i for i in treeno_vec])
label_gamma = (['$\gamma$=%s \n $\\alpha$=%s' % (str(gamma_vec[int(gno_vec[i])]),str(a_vec[int(ano_vec[i])])) for i in range(len(gno_vec)) ])
row_a = len(label_tree)
row_gamma = len(label_gamma)
pos = [0,1]
posnv = [2]
pos_label=['$\gamma$','a','$\\nu$']
# Set up the matplotlib figure
f, axes = plt.subplots(row_gamma,row_a, figsize=(9, 9)) #, sharex=True, sharey=True
gamma_vec_point = np.repeat(gamma_vec[gno_vec],len(label_tree))
a_vec_point = np.repeat(a_vec[ano_vec],len(label_tree))

count = 0
# Rotate the starting point around the cubehelix hue circle
for ax in axes.flat:
    gamma = gamma_list[count]
    a = a_list[count]
    nv = nv_list[count]

    axnv = ax.twinx()  # instantiate a second axes that shares the same x-axis
    if len(gamma)==1:
        axnv.plot()

    else:
        fitness = fit_list[count]
        q5 = np.argsort(fitness)[-population // 20]  # best 5%
        fit_index = np.where(fitness > fitness[q5])[0]
        # d=np.column_stack((gamma,a))
        d = [gamma, a]
        gamma5th = gamma_list[count][fit_index]
        a5th = a_list[count][fit_index]
        nv5th = nv_list[count][fit_index]
        d5th = [gamma5th, a5th]
        dg = [gamma, a, gamma5th, a5th]
        # Generate and plot a random bivariate dataset
        ax.violinplot(dataset=d5th, positions=pos, points=20, widths=0.3,
                    showmeans=False, showextrema=True, showmedians=True)
        ax.scatter(x=0,y=gamma_vec_point[count],color='r',s=10,alpha=1)
        ax.scatter(x=1,y=a_vec_point[count],color='r',s=10,alpha=1)
        axnv.violinplot(dataset=nv5th, positions=posnv, points=20, widths=0.3,
                    showmeans=False, showextrema=True, showmedians=True)
        axnv.scatter(x=2,y=truenv,color='r',s=10,alpha=1)

    if (count+1)%len(label_tree)==0:
        axnv.set_ylabel(label_gamma[int((count+1)/len(label_tree))-1])
        ax.get_yaxis().set_ticks([])

    if count in range(0,row_a):
        ax.title.set_text(label_tree[count])
        # ax.yaxis.set_label_position("right")
        # ax.tick_params(axis='y', which='both', labelleft='off', labelright='on')
    if  (count)%len(label_tree)==0:
        axnv.get_yaxis().set_ticks([])
    if count in range(20):
        ax.get_xaxis().set_ticks([])
    if count in range(18,20):
        plt.xticks([0, 1, 2], pos_label)
    if count >1:
        axnv.yaxis.get_offset_text().set_visible(False)
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

import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="white")
treeno_vec = [1,2,9,3,4,10,5,6,11]
gno_vec = [0,0,1,1,2,2]
ano_vec = [0,1,0,1,0,1]
gamma_vec = np.array([0,0.001,0.01])
a_vec =  np.array([0.5,1])
gamma_list = []
a_list = []
nv_list = []
fit_list = []
gamma_5thmedian_list = []
gamma_5thvar_list = []
a_5thmedian_list = []
a_5thvar_list = []
nu_5thmedian_list = []
nu_5thvar_list = []

truenv=1e-11
lowylim=truenv*(-0.4)
upperylim=2*truenv*6/5
for tree in treeno_vec:
    for indicator_g in range(len(gamma_vec)):
        # g=gamma_vec[indicator_g]
        for indicator_a in [4,5]:
            # a=a_vec[indicator_a]
            if platform.system() == 'Windows':
                file = 'C:/Liang/Googlebox/Research/Project2/smcvdata/tree%d+nv/smc%dg%da.npy' % (tree,indicator_g,indicator_a)
            elif platform.system() == 'Darwin':
                file = '/Users/dudupig/GoogleDrive/Research/Project2/smcvdata/tree%d+nv/smc%dg%da.npy' % (tree,indicator_g,indicator_a)
            if os.path.isfile(file):
                para_data = np.load(file).item()
                generation = len(para_data['gamma'])
                population = len(para_data['gamma'][0])
                gamma = para_data['gamma'][generation - 1]
                a = para_data['a'][generation - 1]
                nu = para_data['nu'][generation - 1]
                fitness = para_data['fitness'][generation - 1]
                q5 = np.argsort(fitness)[-population // 20]  # best 5%
                fit_index = np.where(fitness > fitness[q5])[0]
                gamma5thmedian = np.median(gamma[fit_index])
                a5thmedian = np.median(a[fit_index])
                nu5thmedian = np.median(nu[fit_index])
                gamma5thvar = np.var(gamma[fit_index])
                a5thvar = np.var(a[fit_index])
                nu5thvar = np.var(nu[fit_index])

                gamma_5thmedian_list.append(gamma5thmedian)
                a_5thmedian_list.append(a5thmedian)
                nu_5thmedian_list.append(nu5thmedian)
                gamma_5thvar_list.append(gamma5thvar)
                a_5thvar_list.append(a5thvar)
                nu_5thvar_list.append(nu5thvar)


            else:
                print('Tree = %d; g = %d; a = %d...' % (tree,indicator_g,indicator_a))
                break
#
# label_tree = (['Test %d' % i for i in treeno_vec])
# label_gamma = (['$\gamma$=%s \n $\\alpha$=%s' % (str(gamma_vec[int(gno_vec[i])]),str(a_vec[int(ano_vec[i])])) for i in range(len(gno_vec)) ])

testgroups_label = ['Test group %d' %i for i in range(1,4)]
testgrouppos = [0,1,2]

paracombpos = [0,1,2,3,4,5]

sizegrouppos = [0,1,2]
whole_sizegroups = np.repeat(np.tile(sizegrouppos,len(testgrouppos)),len(paracombpos))

whole_paracombpos = np.tile(paracombpos,len(sizegrouppos)*len(testgrouppos))

whole_testgroup = np.repeat(testgrouppos,len(paracombpos)*len(sizegrouppos))

true_gamma = np.tile(np.tile(gamma_vec,len(a_vec)),len(testgrouppos))
true_a = np.tile(np.repeat(a_vec,len(gamma_vec)),len(testgrouppos))

figrow1 = np.tile(range(6),3)
figrow2 = np.tile(range(6,12),3)
figrow3 = np.tile(range(12,18),3)
whole_figno = np.concatenate((figrow1,figrow2,figrow3))

est_list = {"gammamean":gamma_5thmedian_list,"amean":a_5thmedian_list, "numean":nu_5thmedian_list,
            'gammavar':gamma_5thvar_list,'avar':a_5thvar_list,'nuvar':nu_5thvar_list,
            "size":whole_sizegroups,'paracombo':whole_paracombpos,'testgroup':whole_testgroup,
            'figno':whole_figno}
est_data = pd.DataFrame(data=est_list)

norest_data = est_data

for i in range(len(paracombpos)*len(testgrouppos)):
    gammaest = norest_data[norest_data.figno == i]['gammamean']
    norest_data['figno'][i]=(gammaest-true_gamma[i])/abs(np.max(gammaest-true_gamma[i]))

# Initialize a grid of plots with an Axes for each walk
grid = sns.FacetGrid(est_data, col="figno", hue="paracombo", palette="tab20c",
                     col_wrap=6, height=1.5)

# Draw a horizontal line to show the starting point
grid.map(plt.axhline, y=0, ls=":", c=".5")

# Draw a line plot to show the trajectory of each random walk
grid.map(plt.plot, "size", "gammavar", marker="o")

# Adjust the tick positions and labels
grid.set(xticks=np.arange(5), yticks=[-3, 3],
         xlim=(-.5, 2.5), ylim=(-.5, 1.2))

# Adjust the arrangement of the plots
grid.fig.tight_layout(w_pad=1)



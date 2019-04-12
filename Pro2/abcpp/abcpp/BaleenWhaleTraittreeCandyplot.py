import sys, os
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from PhyloDiff_model_sim import Candimodels
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam
import numpy as np
import matplotlib.pyplot as plt

scalar = 20000
sigma2 = 0.5  # Brownian Motion variance
meantrait = 0.0
# the generating params for models

generating = 'DR'
if generating == 'DR':
    # a = 0.5
    m = 1
elif generating == 'NH':
    a = 0
else:
    print('Please specify the candidate models: DR or NH.')

timegap = 1

# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'



# trait evolution plot

gamma_vec = np.array([0, 0.001, 0.01, 0.1, 0.5, 1])
m_vec = gamma_vec
a_vec = gamma_vec
row_gamma = len(gamma_vec)
count = 0

td = DVTreeData(path=files, scalar=scalar)


f1, axes1 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #

# label_a = (['m=0','m=.001','m=.01','m=.1','m=.5','m=1'])
label_a = (['a=0','a=.001','a=.01','a=.1','a=.5','a=1'])

label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5','$\gamma$=1'])

for index_g in range(len(gamma_vec)):
    gamma1=gamma_vec[index_g]
    for index_a in range(len(m_vec)):
        a=a_vec[index_a]
        print( count)
        candiparam = np.array([gamma1, a, meantrait, m, meantrait, sigma2])

        simresult = Candimodels(td,candiparam,mode = 'f')
        # if pic==0:
        evo_time, total_species = td.sim_evo_time,td.total_species
        trait_RI_dr = simresult['Z']
        trait_RI_dr[trait_RI_dr==0] = np.nan
        num_lines = total_species
        x = np.arange(evo_time/timegap+1)

        labels = []
        for i in range( num_lines ):
            axes1[index_g,index_a].plot(x, trait_RI_dr[::timegap, i])

        # axes[index_g, index_a].yaxis.set_major_locator(plt.NullLocator())
        axes1[index_g, index_a].xaxis.set_major_locator(plt.NullLocator())

        # axes[index_g, index_a].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        # axes[index_g, index_a].set_yscale('log')
        axes1[index_g, index_a].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


        if count in range(0, row_gamma):
            axes1[index_g, index_a].title.set_text(label_a[count])


        if count in ([5, 11, 17, 23, 29, 35]):
            axes1[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
            axes1[index_g, index_a].yaxis.set_label_position("right")

        count += 1


dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree/'
f1.savefig(dir_fig+'TP2w.png')
plt.close(f1)

plt.close('all')


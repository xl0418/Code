import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam
from dvtraitsim_py import DVSim
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from PhyloDiff_model_sim import Candimodels

scalar = 1000
sigma2 = 0.02  # Brownian Motion variance
meantrait = 32.50540571860487

# paras for DR
gamma_DR_mean = 0.638579798972643
a_DR_mean =  0.7627260725192937
m_DR_mean = 18.421144252900813

# paras for NH
gamma_NH_mean = 0.231243
m_NH_mean = 0.014966

# paras for TP
gamma = .000236
a = .617900
nu = 2.411e-11
K = 10e8


timegap = 10

# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)


f1, axes1 = plt.subplots(3, 1, figsize=(9, 9),sharey=True,sharex=True) #

label_model = ['TP','DR','NH']

count = 0
for model in label_model:
    if model == 'TP':
        param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500,
                        initpop_sigma=10.0, break_on_mu=False)
        simresult = DVSim(td, param)

    elif model == 'DR':
        candiparam = np.array([gamma_DR_mean, a_DR_mean, meantrait, m_DR_mean, meantrait, sigma2])
        simresult = Candimodels(td, candiparam,mode = 'f')
    elif model == 'NH':
        candiparam = np.array([gamma_NH_mean, 0, meantrait, m_NH_mean, meantrait, sigma2])
        simresult = Candimodels(td, candiparam,mode = 'f')
    else:
        print('Pls specify the models...')
        break

    # if pic==0:
    evo_time, total_species = td.sim_evo_time,td.total_species
    trait_RI_dr = simresult['Z']
    trait_RI_dr[trait_RI_dr==0] = np.nan
    num_lines = total_species
    x = np.arange(evo_time/timegap)

    labels = []
    for i in range( num_lines ):
        axes1[count].plot(x, trait_RI_dr[::timegap, i])

    # axes[index_g, index_a].yaxis.set_major_locator(plt.NullLocator())

    # axes1[count].xaxis.set_major_locator(plt.NullLocator())

    # axes[index_g, index_a].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    # axes[index_g, index_a].set_yscale('log')
    #
    # axes1[count].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    #
    axes1[count].set_ylabel(label_model[count])
    axes1[count].yaxis.set_label_position("right")

    count += 1

#
# dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'
# f1.savefig(dir_fig+'TP1q.png')
# plt.close(f1)
#
# plt.close('all')


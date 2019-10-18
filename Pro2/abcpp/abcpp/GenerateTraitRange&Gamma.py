import sys, os
import platform

if platform.system() == 'Windows':
    sys.path.append('C:/Liang/abcpp_master8/abcpp')
elif platform.system() == 'Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 10000
K = 10e7
nu = 1 / (100 * K)
timegap = 100
replicate = 1000
# let's try to find a true simulation:


# trait evolution plot
for no_tree in range(1, 23):
    gamma_vec = np.array([0, 0.001, 0.01, 0.1, 0.5, 1])
    a_vec = gamma_vec
    row_gamma = len(gamma_vec)
    count = 0
    tree = 'tree' + '%d' % no_tree
    example = 'example' + '%d' % no_tree
    if platform.system() == 'Windows':
        dir_path = 'c:/Liang/Googlebox/Research/Project2'
        files = dir_path + '/treesim_newexp/' + example + '/'
        td = DVTreeData(path=files, scalar=scalar)
    elif platform.system() == 'Darwin':
        file = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/tree_data/' + example + '/'

    f1, axes1 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9), sharey=True, sharex=True)  #

    label_a = (['$\\alpha$=0', '$\\alpha$=.001', '$\\alpha$=.01', '$\\alpha$=.1', '$\\alpha$=.5', '$\\alpha$=1'])
    label_gamma = (['$\gamma$=0', '$\gamma$=.001', '$\gamma$=.01', '$\gamma$=.1', '$\gamma$=.5', '$\gamma$=1'])

    xticks = (0, td.evo_time * scalar / timegap)
    xlabels = ['0', '150K']
    for index_g in range(len(gamma_vec)):
        gamma1 = gamma_vec[index_g]
        for index_a in range(len(a_vec)):
            a = a_vec[index_a]
            print(count)
            if index_a >= index_g:
                obs_param_single = DVParamLiang(gamma=gamma1, a=a, K=K,h=1, nu=nu, r=r,
                                             theta=theta,
                                         V00=.01, V01=.01,Vmax=1, inittrait=0, initpop=500,
                                    initpop_sigma=10.0, break_on_mu=False)
                obs_param = np.tile(obs_param_single,(replicate, 1))
                simresult = dvcpp.DVSimTVP(td, obs_param)
                valid = np.where(simresult['sim_time'] == td.sim_evo_time)[0]
                if len(valid)>0:
                    pic = 0
                    trait_RI_dr = simresult['Z']
                    trait_range = np.max(trait_RI_dr[valid], axis=1) - np.min(trait_RI_dr[valid], axis=1)
                    sns.distplot(trait_range, kde=False,ax=axes1[index_g, index_a])

                else:
                    pic = 1
                    axes1[index_g, index_a].plot([])

                axes1[index_g, index_a].set_xticks(xticks)
                axes1[index_g, index_a].set_xticklabels(xlabels, minor=False)

                if count in range(0, row_gamma):
                    axes1[index_g, index_a].title.set_text(label_a[count])

                if count in ([5, 11, 17, 23, 29, 35]):
                    axes1[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
                    axes1[index_g, index_a].yaxis.set_label_position("right")

            else:
                axes1[index_g, index_a].plot([])

                axes1[index_g, index_a].axis('off')
            count += 1


    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_treeuppertri/' + tree
    f1.text(0.84, 0.04, 'Trait', ha='center', fontsize=15)
    f1.text(0.04, 0.84, 'Frequency', va='center', rotation='vertical', fontsize=15)

    f1.savefig(dir_fig + 'TP.png')
    plt.close(f1)

    plt.close('all')

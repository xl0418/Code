import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import numpy as np
import matplotlib.pyplot as plt

theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 100000
K=10e8
nu=1/(100*K)
meantrait = 0
# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'



# trait evolution plot

gamma_vec = np.array([0, 0.001, 0.01, 0.1, 0.5, 1])
a_vec = gamma_vec
row_gamma = len(gamma_vec)
count = 0
population = 100
td = DVTreeData(path=files, scalar=scalar)
obs_param = DVParam(gamma=0.01, a=0.5, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500000,
                    initpop_sigma=10.0, break_on_mu=False)
params_TP = np.tile(obs_param, (population, 1))

f1, axes1 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #

label_a = (['a=0','a=.001','a=.01','a=.1','a=.5','a=1'])
label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5','$\gamma$=1'])

for index_g in range(len(gamma_vec)):
    gamma1=gamma_vec[index_g]
    for index_a in range(len(a_vec)):
        a=a_vec[index_a]
        if index_g > index_a:
            axes1[index_g, index_a].plot([])
            axes1[index_g, index_a].axis('off')
        else:
            params_TP[:,0] = gamma1
            params_TP[:,1] = a

            print( count)
            simmodelTP = dvcpp.DVSim(td, params_TP)
            valid_TP = np.where(simmodelTP['sim_time'] == td.sim_evo_time)[0]
            if len(valid_TP)>0:
                Z_modelTP = simmodelTP['Z'][valid_TP]
                axes1[index_g,index_a].hist(Z_modelTP.flatten())
                axes1[index_g,index_a].set_ylim([0, 300])
            else:
                axes1[index_g, index_a].text(0.45,0.45,"X")


            if count in range(0, row_gamma):
                axes1[index_g, index_a].title.set_text(label_a[count])


            if count in ([5, 11, 17, 23, 29, 35]):
                axes1[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
                axes1[index_g, index_a].yaxis.set_label_position("right")

        count += 1


# dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'
#
#
# f1.savefig(dir_fig+'TP2w.png')
# plt.close(f1)


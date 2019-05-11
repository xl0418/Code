import os,sys
sys.path.append('C:/Liang/abcpp_ms2/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
from dvtraitsim_py import DVSimLiang
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from Dvtraitsim_without_population import DVSimLiang_nopop
from Dvtraitsim_metabolism import DVSimMetabolism
from Dvtraitsim_patiallywithoutpop import DVSimLiang_pnopop
from multiprocessing import Pool
from itertools import repeat
num_cores = Pool(8)  # the number of cores

#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'
td = DVTreeData(path=files, scalar=1000)

K = 1
nu = 1e-4
meantrait = 1300
population =100
obs_param = DVParamLiang(gamma=1e-13, a=1e-2, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=10, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
params_TV = np.tile(obs_param, (population, 1))  # duplicate
params_TV[:, 0] = np.random.uniform(0, 1e-13, params_TV.shape[0])  # randomize 'gamma'
params_TV[:, 1] = np.random.uniform(0, 1e-2, params_TV.shape[0])  # randomize 'a'
params_TV[:, 4] = np.random.uniform(0, 1e-4, params_TV.shape[0])  # randomize 'nu'
params_TV[:, 9] = np.random.uniform(0, 10, params_TV.shape[0])  # randomize 'Vm'

for replicate in range(100):
    print(replicate)

    simresult = DVSimLiang_nopop(td, obs_param)

    if simresult['sim_time'] == td.sim_evo_time:
        pic = 0
        break
    else:
        pic = 1

simmodeltv_list = num_cores.starmap(DVSimLiang_pnopop, zip(repeat(td), params_TV))
valid_TV = np.where([simmodeltv_list[i]['sim_time'] == td.sim_evo_time for i in range(population)])[0]
if len(valid_TV) == 0:
    print('No complete results from TV model ')
else:
    for valid_TV_Z in valid_TV:
        Z_modeltv = simmodeltv_list[valid_TV_Z]['Z']
        Z = np.vstack([Z, sorted(Z_modeltv)])



K = 1e13
nu = 1e-4
meantrait = 1300

param_metabolism = DVParamLiang(gamma=1e-6, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e7,
                                initpop_sigma=10.0, break_on_mu=False)

for replicate in range(100):
    print(replicate)

    simresult = DVSimMetabolism(td, param_metabolism)

    if simresult['sim_time'] == td.sim_evo_time:
        pic = 0
        break
    else:
        pic = 1




K = 1e8
nu = 1e-4
meantrait = 1300

param_metabolism = DVParamLiang(gamma=1e-6, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e7,
                                initpop_sigma=10.0, break_on_mu=False)

for replicate in range(100):
    print(replicate)

    simresult = DVSimLiang(td, param_metabolism)

    if simresult['sim_time'] == td.sim_evo_time:
        pic = 0
        break
    else:
        pic = 1

ratio31 = np.array(simresult['Vcomponent3'])/np.array(simresult['Vcomponent1'])
ratio32 = np.array(simresult['Vcomponent3'])/np.array(simresult['Vcomponent2'])
x = np.arange(0,len(ratio31))
import matplotlib.pyplot as plt
plt.scatter(x,ratio31,marker='+',label='V_sele/V_rep')
plt.scatter(x,ratio32,marker='o',label='V_sele/V_mut',s=0.1)
plt.legend()

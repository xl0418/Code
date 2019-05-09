import os,sys
sys.path.append('C:/Liang/abcpp_ms2/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from Dvtraitsim_without_population import DVSimLiang_nopop
from Dvtraitsim_metabolism import DVSimMetabolism


#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'
td = DVTreeData(path=files, scalar=1000)

K = 1e6
nu = 1e-4
meantrait = 1300

obs_param = DVParamLiang(gamma=1e-7, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)

for replicate in range(100):
    print(replicate)

    simresult = DVSimLiang_nopop(td, obs_param)

    if simresult['sim_time'] == td.sim_evo_time:
        pic = 0
        break
    else:
        pic = 1


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



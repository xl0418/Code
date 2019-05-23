import os,sys
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
from Dvtraitsim_TVP import DVSimTVP
from Dvtraitsim_TV import DVSimTV
from Dvtraitsim_TVM import DVSimTVM

#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'
td = DVTreeData(path=files, scalar=2000)

K = 1e6
nu = 1e-4
meantrait = 1300
gamma_TVM_est =  7.539611977498784e-05 #np.mean(gamma_TVM[fit_index_TVM])
a_TVM_est = 0.0021762021948271803 #np.mean(a_TVM[fit_index_TVM])
nu_TVM_est =0.0036982608897495085 # np.mean(nu_TVM[fit_index_TVM])
vm_TVM_est = 42.941# np.mean(vm_TVM[fit_index_TVM])

# obs_param = DVParamLiang(gamma=1e-7, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e5,
#                                 initpop_sigma=10.0, break_on_mu=False)
obs_param = DVParamLiang(gamma=gamma_TVM_est, a=a_TVM_est, K=1e15, h=1, nu=nu_TVM_est, r=1, theta=meantrait,
                         V00=.5,V01=.5, Vmax=vm_TVM_est, inittrait=meantrait, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)

simresult = DVSimTVM(td, obs_param)
simresult['Z'][-1]
#
#
#
#
# K = 1e13
# nu = 1e-4
# meantrait = 1300
#
# param_metabolism = DVParamLiang(gamma=1e-6, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e7,
#                                 initpop_sigma=10.0, break_on_mu=False)
#
# for replicate in range(100):
#     print(replicate)
#
#     simresult = DVSimMetabolism(td, param_metabolism)
#
#     if simresult['sim_time'] == td.sim_evo_time:
#         pic = 0
#         break
#     else:
#         pic = 1
#
#
#
#
# K = 1e8
# nu = 1e-4
# meantrait = 1300
#
# param_metabolism = DVParamLiang(gamma=1e-6, a=1e-4, K=K,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e7,
#                                 initpop_sigma=10.0, break_on_mu=False)
#
# for replicate in range(100):
#     print(replicate)
#
#     simresult = DVSimLiang(td, param_metabolism)
#
#     if simresult['sim_time'] == td.sim_evo_time:
#         pic = 0
#         break
#     else:
#         pic = 1
#
# ratio31 = np.array(simresult['Vcomponent3'])/np.array(simresult['Vcomponent1'])
# ratio32 = np.array(simresult['Vcomponent3'])/np.array(simresult['Vcomponent2'])
# x = np.arange(0,len(ratio31))
# import matplotlib.pyplot as plt
# plt.scatter(x,ratio31,marker='+',label='V_sele/V_rep')
# plt.scatter(x,ratio32,marker='o',label='V_sele/V_mut',s=0.1)
# plt.legend()

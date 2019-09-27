import sys, os
sys.path.append('C:/Liang/abcpp_ms7/abcpp')
from Dvtraitsim_TVM import DVSimTVMLog10
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
from matplotlib.pylab import *
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
theta = 3  # optimum of natural selection
gamma = 1e-7 # intensity of natural selection
r = 1  # growth rate
a = 5e-6 # intensity of competition
K = 10e8  # carrying capacity
nu=1.1e-03
Vmax = 0.01
scalar = 20000


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)

# parameter settings
obs_param = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.1,V01=.1, Vmax=Vmax, inittrait=theta, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)


population = 1000

obs_param_TVMlog10 = np.tile(obs_param, (population, 1))
simmodelTVM = dvcpp.DVSimTVMLog10(td, obs_param_TVMlog10)
valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
len(simmodelTVM['Z'][valid_TVM])
simmodelTVM['Z'][valid_TVM].max()
simmodelTVM['Z'][valid_TVM].min()


for rep in range(100):
    simresult = DVSimTVMLog10(td,obs_param)
    if simresult['sim_time'] == td.sim_evo_time:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)

simresult['Z']
np.var(simresult['Z'])


params_TVM = np.tile(obs_param, (population, 1))  # duplicate
params_TVM[:, 0] = np.random.uniform(0, 1e-7, params_TVM.shape[0])  # randomize 'gamma'
params_TVM[:, 1] = np.random.uniform(1e-7, 1e-5, params_TVM.shape[0])  # randomize 'a'
params_TVM[:, 4] = np.random.uniform(1e-5, 1e-3, params_TVM.shape[0])  # randomize 'nu'
params_TVM[:, 6] = np.random.uniform(2.5, 3.3, params_TVM.shape[0])  # randomize 'theta'
params_TVM[:, 9] = np.random.uniform(0, 1, params_TVM.shape[0])  # randomize 'Vm'

simmodelTVM = dvcpp.DVSimTVMLog10(td, params_TVM)
valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
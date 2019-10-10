import sys, os
sys.path.append('C:/Liang/abcpp_ms8/abcpp')
from Dvtraitsim_TVM import DVSimTVMLog10, DVSimTVM
from Dvtraitsim_TVP import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import numpy as np
theta = 2.8  # optimum of natural selection
gamma =1 # intensity of natural selection
r = 1  # growth rate
a =1200# intensity of competition   1100-1700
K = 10e11  # carrying capacity
nu=6.957e-03
Vmax = 0.836507e-4
scalar = 20000


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)

# parameter settings
obs_param = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.001,V01=.001,
                         Vmax=Vmax, inittrait=theta, initpop=1000,
                initpop_sigma = 10.0, break_on_mu=False)


population = 100

obs_param_TVMlog10 = np.tile(obs_param, (population, 1))
simmodelTVMlog10 = dvcpp.DVSimTVMLog10(td, obs_param_TVMlog10)
valid_TVMlog10 = np.where(simmodelTVMlog10['sim_time'] == td.sim_evo_time)[0]
if len(simmodelTVMlog10['Z'][valid_TVMlog10])>0:
    simmodelTVMlog10['Z'][valid_TVMlog10].max()
    simmodelTVMlog10['Z'][valid_TVMlog10].min()
else:
    print('No valid simulations')



theta = 1300  # optimum of natural selection
gamma =4.901e-9 # intensity of natural selection
r = 1  # growth rate
a = 3.531e-4# intensity of competition
K = 10e9  # carrying capacity
nu=6.957e-03
Vmax = 0.836507e1
scalar = 20000

obs_param = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.001,V01=.001,
                         Vmax=Vmax, inittrait=theta, initpop=1000,
                initpop_sigma = 10.0, break_on_mu=False)

population = 1000

obs_param_TVMlog10 = np.tile(obs_param, (population, 1))
simmodelTVM = dvcpp.DVSimTVM(td, obs_param_TVMlog10)
valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
if len(simmodelTVM['Z'][valid_TVM])>0:
    simmodelTVM['Z'][valid_TVM].max()
    simmodelTVM['Z'][valid_TVM].min()
else:
    print('No valid simulations')


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
params_TVM[:, 4] = np.random.uniform(1e-5, 1e-2, params_TVM.shape[0])  # randomize 'nu'
params_TVM[:, 6] = np.random.uniform(2.5, 3.3, params_TVM.shape[0])  # randomize 'theta'
params_TVM[:, 9] = np.random.uniform(0, 1, params_TVM.shape[0])  # randomize 'Vm'

simmodelTVM = dvcpp.DVSimTVMLog10(td, params_TVM)
valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
simmodelTVM['Z'][valid_TVM]


# Test TVP
theta = 3  # optimum of natural selection
gamma =4.629e-01 # intensity of natural selection
r = 1  # growth rate
a = 7.077e-01 # intensity of competition
K = 10e5  # carrying capacity
nu=3.694e-02
Vmax = 0.067450
scalar = 20000


# parameter settings
obs_param_TVP = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.1,V01=.1,
                          Vmax=Vmax, inittrait=theta, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)

for rep in range(10):
    simresult = DVSimTVP(td,obs_param_TVP)
    if simresult['sim_time'] == td.sim_evo_time:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)

simresult['Z']
simresult['N']

np.var(simresult['Z'])




# Test TVP
theta = 3.05549  # optimum of natural selection

gamma =8.665e+00# intensity of natural selection
r = 1  # growth rate
a = 2.554e+02   # intensity of competition
K = 10e5  # carrying capacity
nu=1.204e-03
Vmax = 0.000099
scalar = 20000

for gamma in range(1,20):
    gamma=gamma*0.1
    for a in range(20,200,10):
        # parameter settings
        obs_param_TVP = DVParamLiang(gamma=gamma, a=a, K=K,h=0.5, nu=nu, r=1, theta=theta,V00=1e-5,
                                     V01=1e-5,
                                  Vmax=Vmax, inittrait=theta, initpop=500,
                        initpop_sigma = 10.0, break_on_mu=False)

        population = 100

        obs_param_TVP_batch = np.tile(obs_param_TVP, (population, 1))
        simmodelTVP = dvcpp.DVSimTVP(td, obs_param_TVP_batch)
        valid_TVP = np.where(simmodelTVP['sim_time'] == td.sim_evo_time)[0]
        if len(simmodelTVP['Z'][valid_TVP])>0:
            print(gamma,a)
            simmodelTVP['Z'][valid_TVP].max()
            simmodelTVP['Z'][valid_TVP].min()
        # else:
        #     print('No valid simulations')


# Test tvp simulation on cluster. V00 V01 should be small. Otherwise, no complete simulations.
K=1e6
lefttrait=2.7
righttrait=3.3
prior = [3, 5.1, 80, 270, 0.0, nu*100, 0, 1e-3, lefttrait, righttrait]
sampleparam_TVP = DVParamLiang(gamma=1, a=1, K=K,h=1, nu=nu, r=1, theta=0,V00=.001,V01=.001,
                               Vmax=100, inittrait=3, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
params_TVP = np.tile(sampleparam_TVP, (population, 1))  # duplicate
params_TVP[:, 0] = np.random.uniform(prior[0], prior[1], params_TVP.shape[0])  # randomize 'gamma'
params_TVP[:, 1] = np.random.uniform(prior[2], prior[3], params_TVP.shape[0])  # randomize 'a'
params_TVP[:, 4] = np.random.uniform(prior[4], prior[5], params_TVP.shape[0])  # randomize 'nu'
params_TVP[:, 6] = np.random.uniform(prior[8], prior[9], params_TVP.shape[0])  # randomize 'theta'
params_TVP[:, 9] = np.random.uniform(prior[6], prior[7], params_TVP.shape[0])  # randomize 'Vm'

simmodelTVP = dvcpp.DVSimTVP(td, params_TVP)
valid_TVP = np.where(simmodelTVP['sim_time'] == td.sim_evo_time)[0]
if len(simmodelTVP['Z'][valid_TVP]) > 0:
    print(gamma, a)
    simmodelTVP['Z'][valid_TVP].max()
    simmodelTVP['Z'][valid_TVP].min()
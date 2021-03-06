import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
sys.path.append('C:/Liang/abcpp_ms10/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import csv
population = 40000
timescaling = 20000
total_population = population * 3
data_dir = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'
data_name = data_dir + 'results_ms_1107_continue40/pics_continue_to_40.npy'
obs_file = data_dir + 'results_ms_1107_continue40/pics_con_sim/'
tree_data_file = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'

iteration_of_est = -3
est_data = np.load(data_name,allow_pickle=True).item()
fitness = est_data['fitness'][iteration_of_est]
Z = est_data['Z']
gamma_TVP = est_data['gamma_data_TVP'][iteration_of_est]
a_TVP = est_data['a_data_TVP'][iteration_of_est]
nu_TVP = est_data['nu_data_TVP'][iteration_of_est]
vm_TVP = est_data['vm_data_TVP'][iteration_of_est]
theta_TVP = est_data['theta_data_TVP'][iteration_of_est]

gamma_TV = est_data['gamma_data_TV'][iteration_of_est]
a_TV = est_data['a_data_TV'][iteration_of_est]
nu_TV = est_data['nu_data_TV'][iteration_of_est]
vm_TV = est_data['vm_data_TV'][iteration_of_est]
theta_TV = est_data['theta_data_TV'][iteration_of_est]

gamma_TVM = est_data['gamma_data_TVM'][iteration_of_est]
a_TVM = est_data['a_data_TVM'][iteration_of_est]
nu_TVM = est_data['nu_data_TVM'][iteration_of_est]
vm_TVM = est_data['vm_data_TVM'][iteration_of_est]
theta_TVM = est_data['theta_data_TVM'][iteration_of_est]

q5_TVP = np.argsort(fitness[ :population ])[-int(population // 200+1)]  # best 5%
q5_TV = np.argsort(fitness[population:2 * population ])[-int(population // 200+1)] + population
# best 5%
q5_TVM = np.argsort(fitness[2 * population:3 * population])[-int(population // 200+1)] + 2 * \
         population  # best 5%

fit_index_TVP = np.where(fitness[ :population] > fitness[ q5_TVP])[0]
fit_index_TV = np.where(fitness[ population:2 * population] > fitness[ q5_TV])[0]
fit_index_TVM = np.where(fitness[ 2 * population:] > fitness[ q5_TVM])[0]

gamma_TVP_est = np.mean(gamma_TVP[fit_index_TVP])
a_TVP_est = np.mean(a_TVP[fit_index_TVP])
nu_TVP_est = np.mean(nu_TVP[fit_index_TVP])
vm_TVP_est = np.mean(vm_TVP[fit_index_TVP])
theta_TVP_est = np.mean(theta_TVP[fit_index_TVP])

gamma_TV_est = np.mean(gamma_TV[fit_index_TV])
a_TV_est = np.mean(a_TV[fit_index_TV])
nu_TV_est = np.mean(nu_TV[fit_index_TV])
vm_TV_est = np.mean(vm_TV[fit_index_TV])
theta_TV_est = np.mean(theta_TV[fit_index_TV])

gamma_TVM_est =  np.mean(gamma_TVM[fit_index_TVM])
a_TVM_est = np.mean(a_TVM[fit_index_TVM])
nu_TVM_est =np.mean(nu_TVM[fit_index_TVM])
vm_TVM_est =  np.mean(vm_TVM[fit_index_TVM])
theta_TVM_est = np.mean(theta_TVM[fit_index_TVM])

td = DVTreeData(path=tree_data_file, scalar=timescaling)

def simtraits(param, replicates,obs_dir,tree_data_file,timescaling,mode):
    td = DVTreeData(path=tree_data_file, scalar=timescaling)
    params = np.tile(param, (replicates, 1))  # duplicate
    if mode == 'TVP':
        predictsim = dvcpp.DVSimTVP(td, params)
    elif mode == 'TV':
        predictsim = dvcpp.DVSimTV(td, params)
    elif mode == 'TVM':
        predictsim = dvcpp.DVSimTVMLog10(td, params)
    else:
        return print('Please sepcify mode...')
    valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]
    print(len(valid))
    if len(valid)>0:
        Z = predictsim['Z'][valid]
        Z = np.nan_to_num(Z)
        Z_df = pd.DataFrame(Z)
        savefilename = obs_dir + 'predictsim%s.csv' % mode
        Z_df.to_csv(savefilename, sep=',', index=False)
        V = predictsim['V'][valid]
        V = np.nan_to_num(V)
        V_df = pd.DataFrame(V)
        savefilename_V = obs_dir + 'predictsimV%s.csv' % mode
        V_df.to_csv(savefilename_V, sep=',', index=False)
        if mode != 'TV':
            N = predictsim['N'][valid]
            N = np.nan_to_num(N)
            N_df = pd.DataFrame(N)
            savefilename = obs_dir + 'predictsimN%s.csv' % mode
            N_df.to_csv(savefilename, sep=',', index=False)
    else:
        return print('No valid simulation...')


with open(tree_data_file+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(tree_data_file+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = logTL
obsZ = length
s = np.argsort(obsZ)
obsZ = obsZ[s]
obsZ = obsZ.astype(np.float)
meantrait = np.mean(obsZ)
# sim TVP
# gamma_TVP_est=1e-2
# a_TVP_est=5
# nu_TVP_est=1.439e-03
# theta_TVP_est=3.078976
# vm_TVP_est=1e-2

param_TVP = DVParamLiang(gamma=gamma_TVP_est, a=a_TVP_est, K=1e6, h=1, nu=nu_TVP_est, r=1,
                         theta=theta_TVP_est,
                         V00=.0001,V01=.0001, Vmax=vm_TVP_est, inittrait=meantrait, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)

simtraits(param = param_TVP,replicates=1000,obs_dir= obs_file,
          tree_data_file=tree_data_file,timescaling=timescaling,mode = 'TVP')
#
# gamma_TV_est=1e-2
# a_TV_est=1e-2
# nu_TV_est=1.439e-03
# theta_TV_est=3.078976
# vm_TV_est=10.195081

# sim TV
param_TV = DVParamLiang(gamma=gamma_TV_est, a=a_TV_est, K=1e6, h=1, nu=nu_TV_est, r=1, theta=theta_TV_est,
                         V00=.0001,V01=.0001, Vmax=vm_TV_est, inittrait=meantrait, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)

simtraits(param = param_TV,replicates=1000,obs_dir= obs_file,tree_data_file=tree_data_file,
          timescaling=timescaling,mode = 'TV')
#
# gamma_TVM_est=4.901e-15
# a_TVM_est=3.531e-13
# nu_TVM_est=6.957e-04
# theta_TVM_est=3.078976
# vm_TVM_est=0.836507

# sim TVM
param_TVM = DVParamLiang(gamma=gamma_TVM_est, a=a_TVM_est, K=1e12, h=1, nu=nu_TVM_est, r=1, theta=theta_TVM_est,
                         V00=.0001,V01=.0001, Vmax=vm_TVM_est, inittrait=meantrait, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)

simtraits(param = param_TVM,replicates=1000,obs_dir= obs_file,tree_data_file=tree_data_file,\
timescaling=timescaling,mode = 'TVM')


model_index = np.array([0,1,2])
propose_model = np.repeat(model_index,repeats = population)

q5 = np.argsort(fitness)[-int(total_population // 3)]  # best 25%
fit_index = np.where(fitness > fitness[q5])[0]

modelTVPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
modelTVperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)
modelTVMperc = len(np.where(propose_model[fit_index] == 2)[0]) / len(fit_index)

print('Iteration 25th Model TVP: %.1f%% ;  Model TV: %.1f%% ; Model TVM: %.1f%%...'
      % ( modelTVPperc * 100, modelTVperc * 100, modelTVMperc * 100))
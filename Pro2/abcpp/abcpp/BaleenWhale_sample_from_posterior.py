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
sstats = 'pics'
data_name = data_dir + 'results_ms_posterior/'+sstats+'_con40.npy'
obs_file = data_dir + 'results_ms_posterior/'+sstats+'_con_sim/'
tree_data_file = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'
td = DVTreeData(path=tree_data_file, scalar=timescaling)

with open(tree_data_file + 'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(tree_data_file + 'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:, 1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:, 0] == label)[0][0])

logTL = lengthdata_array[length_index, 1].astype(np.float)
length = logTL

# reorder the traits according to the simulation order
sim_species_label = ["Balaena_mysticetus", "Balaenoptera_acutorostrata", "Caperea_marginata",
                     "Balaenoptera_borealis",
                     "Balaenoptera_physalus", "Eschrichtius_robustus", "Balaenoptera_musculus",
                     "Balaenoptera_omurai",
                     "Eubalaena_australis", "Megaptera_novaeangliae", "Balaenoptera_bonaerensis",
                     "Balaenoptera_brydei",
                     "Balaenoptera_edeni", "Eubalaena_glacialis", "Eubalaena_japonica"]

obsZ_ordered_sim = length[
    [np.where(sim_species_label[i] == extantlabels_array)[0][0] for i in range(15)]]
obsZ = obsZ_ordered_sim
meantrait = np.mean(obsZ)

iteration_of_est = -2
est_data = np.load(data_name,allow_pickle=True).item()
fitness = est_data['fitness'][-1]
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


K_TVP=1e6
K_TV = 1e6
K_TVM = 1e12
nu=1e-4

# let's try to find a true simulation:
sampleparam_TVP = DVParamLiang(gamma=1, a=1, K=K_TVP, h=1, nu=nu, r=1, theta=0, V00=.0001,
                               V01=.0001, Vmax=100, inittrait=meantrait, initpop=1e5,
                               initpop_sigma=10.0, break_on_mu=False)
sampleparam_TV = DVParamLiang(gamma=1, a=1, K=K_TV, h=1, nu=nu, r=1, theta=0, V00=.0001, V01=.0001,
                              Vmax=100, inittrait=meantrait, initpop=1e5,
                              initpop_sigma=10.0, break_on_mu=False)
sampleparam_TVM = DVParamLiang(gamma=1, a=1, K=K_TVM, h=1, nu=nu, r=1, theta=0, V00=.0001,
                               V01=.0001, Vmax=100.0, inittrait=meantrait, initpop=1e5,
                               initpop_sigma=10.0, break_on_mu=False)

fitness_TVP = fitness[:population]
fitness_TV = fitness[population:2*population]
fitness_TVM = fitness[2*population:]
sample_size = 1000
index_s = [i for i in range(population)]
sample_index_TVP = np.random.choice(a = index_s,size = sample_size,
                                    p = fitness_TVP/np.sum(fitness_TVP))
sample_index_TV = np.random.choice(a = index_s,size = sample_size,
                                    p = fitness_TV/np.sum(fitness_TV))
sample_index_TVM = np.random.choice(a = index_s,size = sample_size,
                                    p = fitness_TVM/np.sum(fitness_TVM))
params_TVP = np.tile(sampleparam_TVP, (sample_size, 1))  # duplicate
params_TVP[:, 0] = gamma_TVP[sample_index_TVP] # randomize 'gamma'
params_TVP[:, 1] = a_TVP[sample_index_TVP]  # randomize 'a'
params_TVP[:, 4] = nu_TVP[sample_index_TVP]  # randomize 'nu'
params_TVP[:, 6] = theta_TVP[sample_index_TVP]  # randomize 'theta'
params_TVP[:, 9] = vm_TVP[sample_index_TVP]  # randomize 'Vm'

params_TV = np.tile(sampleparam_TV, (sample_size, 1))  # duplicate
params_TV[:, 0] = gamma_TV[sample_index_TV]  # randomize 'gamma'
params_TV[:, 1] = a_TV[sample_index_TV]  # randomize 'a'
params_TV[:, 4] = nu_TV[sample_index_TV]  # randomize 'nu'
params_TV[:, 6] = theta_TV[sample_index_TV]  # randomize 'theta'
params_TV[:, 9] = vm_TV[sample_index_TV]  # randomize 'Vm'

params_TVM = np.tile(sampleparam_TVM, (sample_size, 1))  # duplicate
params_TVM[:, 0] = gamma_TVM[sample_index_TVM]  # randomize 'gamma'
params_TVM[:, 1] = a_TVM[sample_index_TVM]  # randomize 'a'
params_TVM[:, 4] = nu_TVM[sample_index_TVM]  # randomize 'nu'
params_TVM[:, 6] = theta_TVM[sample_index_TVM]   # randomize 'theta'
params_TVM[:, 9] = vm_TVM[sample_index_TVM] # randomize 'Vm'


def simtraits(params,obs_dir,tree_data_file,timescaling,mode):
    td = DVTreeData(path=tree_data_file, scalar=timescaling)
    if mode == 'TVP':
        predictsim = dvcpp.DVSimTVP(td, params)
    elif mode == 'TV':
        predictsim = dvcpp.DVSimTV(td, params)
    elif mode == 'TVM':
        predictsim = dvcpp.DVSimTVMLog10(td, params)
    else:
        return print('Please specify mode...')
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


simtraits(params = params_TVP,obs_dir= obs_file,
          tree_data_file=tree_data_file,timescaling=timescaling,mode = 'TVP')
#
# sim TV

simtraits(params = params_TV,obs_dir= obs_file,tree_data_file=tree_data_file,
          timescaling=timescaling,mode = 'TV')
#
# sim TVM

simtraits(params = params_TVM,obs_dir= obs_file,tree_data_file=tree_data_file,
timescaling=timescaling,mode = 'TVM')

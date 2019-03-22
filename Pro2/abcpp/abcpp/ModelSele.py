import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
from PhyloDiff_model_sim import Candimodels
from tp_update import tp_update
from bm_update import bm_update
from ou_update import ou_update
from dr_update import dr_update


import csv

# gamma = 0.001
# a = 0.1
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err

full =0

#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/'
if full==1:
    files = dir_path + 'full_data/'
else:
    files = dir_path + 'pruend_data/'

td = DVTreeData(path=files, scalar=1000)

prior = [0.5, 0.5, 0.5, 0.5,1e-11,1e-11]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
nu_prior_mean = prior[4]
nu_prior_var = prior[5]
K=10e8
nu=1/(100*K)
num_models = 4

with open(dir_path+'diameter.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    diameterdata = list(csv_reader)

ddarray = np.array(diameterdata)
tipslabel = ddarray[1:,0]
obsZ = ddarray[1:,1]
obsV = ddarray[1:,2]
obsN = ddarray[1:,3]
print('trying to estimate the parameters','...')
if full==1:
    missingdata = np.array([ 57, 160, 165, 163, 182])-1  # missing index in full data
else:
    missingdata = np.array([ 1, 15, 19, 18, 29])-1  # missing index in reconstructive data

naindex = np.where(obsZ == 'NA')[0]
obsZ = np.delete(obsZ,naindex)
s = np.argsort(obsZ)
obsZ = obsZ[s]
obsZ = obsZ.astype(np.float)

meantrait = np.mean(obsZ)
# let's try to find a true simulation:
obs_param = DVParam(gamma=0.01, a=0.5, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500000,
             initpop_sigma = 10.0, break_on_mu=False)
candiparam = np.array([0, 0, meantrait, 1, meantrait, 1])
# pop = dvcpp.DVSim(td, obs_param)

modelnum = 4
population = 5
total_population = population * modelnum
generations = 20
params_TP = np.tile(obs_param, (population, 1))  # duplicate
params_TP[:, 0] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'gamma'
params_TP[:, 1] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'a'
params_TP[:, 3] = np.random.uniform(0.0, 1e-10, params_TP.shape[0])  # randomize 'nu'

params_BM = np.tile(candiparam,(population,1))
params_BM[:,3] = 0
params_BM[:, 5] = np.random.uniform(0.0, 10.0, params_BM.shape[0])  # randomize delta


params_OU = np.tile(candiparam,(population,1))
params_OU[:,0]= np.random.uniform(0.0, 1.0, params_OU.shape[0])
params_OU[:,3]= 0
params_OU[:, 5] = np.random.uniform(0.0, 10.0, params_OU.shape[0])  # randomize delta


params_drury =  np.tile(candiparam,(population,1))
params_drury[:,0]= np.random.uniform(0.0, 1.0, params_drury.shape[0])
params_drury[:,1]= np.random.uniform(0.0, 1.0, params_drury.shape[0])
params_drury[:, 5] = np.random.uniform(0.0, 10.0, params_drury.shape[0])  # randomize delta

# model choice
model_index = np.array(range(num_models))
model_params = np.repeat(model_index,repeats = population)
model_data = np.zeros(shape=(generations, total_population))

gamma_data = np.zeros(shape=(generations, population))
a_data = np.zeros(shape=(generations, population))
nu_data = np.zeros(shape=(generations, population))

fitness= np.zeros(shape=(generations, total_population))
# Initialize the weights.
weight_model = np.zeros(total_population)
weight_model.fill(1 / total_population)

# weights for paras of TP
weight_gamma_TP = np.zeros(population)
weight_gamma_TP.fill(1 / population)
weight_a_TP = np.zeros(population)
weight_a_TP.fill(1 / population)
weight_nu_TP = np.zeros(population)
weight_nu_TP.fill(1 / population)

# weights for paras of BM
weight_del_BM = np.zeros(population)
weight_del_BM.fill(1 / population)


# weights for paras of OU
weight_gamma_OU = np.zeros(population)
weight_gamma_OU.fill(1 / population)
weight_del_OU = np.zeros(population)
weight_del_OU.fill(1 / population)


# weights for paras of dr
weight_gamma_dr = np.zeros(population)
weight_gamma_dr.fill(1 / population)
weight_a_dr = np.zeros(population)
weight_a_dr.fill(1 / population)
weight_del_dr = np.zeros(population)
weight_del_dr.fill(1 / population)

model_data[0, :] = model_params

for g in range(generations):
    # model 0
    simmodelTP = dvcpp.DVSim(td, params_TP)
    # access fitness
    # fitness = np.zeros(population)
    valid = np.where(simmodelTP['sim_time'] == td.sim_evo_time)[0]
    if len(valid)<20:
        print("WARNING:Valid simulations are too scarce!")
    if valid.size > 0:
        Z_modelTP = simmodelTP['Z'][valid]
        Z_modelTP = np.delete(Z_modelTP,np.s_[missingdata],axis=1)
        i, j = argsort2D(Z_modelTP)
        Z_modelTP = Z_modelTP[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTP = np.nan_to_num(Z_modelTP)
        # V = np.nan_to_num(V)
        #GOF: Goodness of fit
    Z = Z_modelTP
    # model 1
    for param_BM in params_BM:
        simmodelBM = Candimodels(td, param_BM)
    # access fitness
    # fitness = np.zeros(population)
        Z_modelBM = simmodelBM['Z'][simmodelBM['Z'].shape[0]-1,:]
        Z_modelBM = np.delete(Z_modelBM,np.s_[missingdata],axis=0)

        Z = np.vstack([Z,Z_modelBM])


    # model 2
    for param_OU in params_OU:
        simmodelOU = Candimodels(td, param_OU)
    # access fitness
    # fitness = np.zeros(population)
        Z_modelOU = simmodelOU['Z'][simmodelOU['Z'].shape[0]-1,:]
        Z_modelOU = np.delete(Z_modelOU,np.s_[missingdata],axis=0)

        Z = np.vstack([Z,Z_modelOU])

    # model 3
    for param_drury in params_drury:
        simmodeldr= Candimodels(td, param_drury)
    # access fitness
    # fitness = np.zeros(population)
        Z_modeldr = simmodeldr['Z'][simmodeldr['Z'].shape[0]-1,:]
        Z_modeldr = np.delete(Z_modeldr,np.s_[missingdata],axis=0)

        Z = np.vstack([Z,Z_modeldr])

    valid = np.concatenate([valid, range(len(np.where(model_params==0)[0]), total_population)])
    fitness[g,valid] += 1.0 - normalized_norm(Z, obsZ)
        # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g,:])[-total_population// 2]  # best 50%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]

    modelTPperc = len(np.where(model_params[fit_index]==0)[0])/len(fit_index)
    modelBMperc = len(np.where(model_params[fit_index]==1)[0])/len(fit_index)
    modelOUperc = len(np.where(model_params[fit_index]==2)[0])/len(fit_index)
    modeldrperc = len(np.where(model_params[fit_index]==3)[0])/len(fit_index)

    print('Iteration = %d 5th Model TP: %.1f%% ; Model TP: %.1f%% ; Model TP: %.1f%% ; Model TP: %.1f%%...'
          % (g,modelTPperc*100 , modelBMperc*100,modelOUperc*100 ,modeldrperc*100))

    # reevaluate the weight of the best fitted  models
    weight_model = weight_model[fit_index]/sum(weight_model[fit_index])
    # sample new models from the fitness of previous best fitted models
    sample_model_index = sorted(np.random.choice(fit_index, total_population, p=weight_model))
    propose_model = model_params[sample_model_index]

    previous_bestfitted_model = model_params[fit_index]
    model_data[g,:] = propose_model

    # update TP paras and weights
    weight_gamma_TP,weight_a_TP,weight_nu_TP,propose_gamma_TP,propose_a_TP,propose_nu_TP=\
        tp_update(previous_bestfitted_model, propose_model, params_TP, weight_gamma_TP, weight_a_TP, weight_nu_TP)
    modelTP = np.where(propose_model==0)
    params_TP = np.tile(obs_param, (len(modelTP[0]), 1))
    params_TP[:, 0] = propose_gamma_TP
    params_TP[:, 1] = propose_a_TP
    params_TP[:, 3] = propose_nu_TP

    # update BM paras and weights
    weight_del_BM, propose_del_BM = \
        bm_update(previous_bestfitted_model, propose_model, params_BM, weight_del_BM)
    modelBM = np.where(propose_model == 1)
    params_BM = np.tile(candiparam, (len(modelBM[0]), 1))
    params_BM[:, 5] = propose_del_BM

    # update OU paras and weights
    weight_gamma_OU, weight_del_OU, propose_gamma_OU,propose_del_OU = \
        ou_update(previous_bestfitted_model, propose_model, params_OU, weight_gamma_OU,weight_del_OU)
    modelOU = np.where(propose_model == 2)
    params_OU = np.tile(candiparam, (len(modelOU[0]), 1))
    params_OU[:, 0] = propose_gamma_OU
    params_OU[:, 5] = propose_del_OU
#


para_data = {'gamma': gamma_data, 'a': a_data, 'nu': nu_data,'fitness': fitness}
# file='C:/Liang/Code/Pro2/abcpp/abcpp/smcndata/'
# # file = '/home/p274981/abcpp/abcpp/'
# file2 = file + 'paraT11.npy'
# np.save(file2,para_data)
#
# smc = np.load(file2).item()
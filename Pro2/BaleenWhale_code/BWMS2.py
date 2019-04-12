import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang,DVParamDrury
import dvtraitsim_cpp as dvcpp
from tp_update import tp_update
from dr_update import dr_update
import csv
import sys

#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


#full tree
dir_path = '/home/p274981/abcpp/'

files = dir_path + 'BaleenWhales/treedata/'

time_scaling = float(sys.argv[1]) # time scaling parameter
dividing = float(sys.argv[2]) # the number dividing the total length
fix_m = float(sys.argv[3]) # fix m (1) or not (0)

savedir = dir_path + 'BaleenWhales/BWMSt%i_d%i_f%i.npy' % (int(time_scaling),int(dividing),int(fix_m))

td = DVTreeData(path=files, scalar=time_scaling)
allowmodeldie = 'off'
K=10e8
nu=1/(100*K)
del_mute ='on'
prior = [0.5, 0.5, 0.5, 0.5,nu,nu]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
nu_prior_mean = prior[4]
nu_prior_var = prior[5]

with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(files+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = 10**logTL/dividing
obsZ = length
print('trying to estimate the parameters','...')

s = np.argsort(obsZ)
obsZ = obsZ[s]
obsZ = obsZ.astype(np.float)
meantrait = np.mean(obsZ)
modelnum = 2
sigma2 = 0.5
population = 10000  # number of the particles for each iteration
total_population = population * modelnum
generations = 50  # number of the iterations
# let's try to find a true simulation:
# The generating paras

print('Try to find the generating model...' )
# let's try to find a true simulation:
param_Liang = DVParamLiang(gamma=0.01, a=0.5,K=K, h=1.0,nu=nu, r=1, theta=meantrait, V00 = .5, V01=.5,
                           Vmax=1, inittrait=meantrait, initpop=500000,
                    initpop_sigma=10.0, break_on_mu=False)
param_Drury = DVParamDrury(gamma = 0.0, a=0.0, theta=meantrait, m=1.0, inittrait=meantrait, var_trait=sigma2)
# pop = dvcpp.DVSim(td, obs_param)


params_TP = np.tile(param_Liang, (population, 1))  # duplicate
params_TP[:, 0] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'gamma'
params_TP[:, 1] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'a'
params_TP[:, 4] = np.random.uniform(0.0, 2 * nu, params_TP.shape[0])  # randomize 'nu'

params_DR = np.tile(param_Drury, (population, 1))
params_DR[:, 0] = np.random.uniform(0.0, 1.0, params_DR.shape[0])  # randomize 'gamma'
params_DR[:, 1] = np.random.uniform(0.0, 1.0, params_DR.shape[0])  # randomize 'a'
params_DR[:, 3] = 1.0 #np.random.uniform(0.0, 5.0, params_DR.shape[0])  # randomize 'm'
params_DR[:, 5] = np.random.uniform(0.0, sigma2, params_DR.shape[0])  # randomize delta

# model choice
model_index = np.array([0, 1])
model_params = np.repeat(model_index, repeats=population)
model_data = np.zeros(shape=(generations + 1, total_population))
model_data[0, :] = model_params
propose_model = model_params

# store parameters used
# TP
gamma_data_TP = np.zeros(shape=(generations + 1, population))
a_data_TP = np.zeros(shape=(generations + 1, population))
nu_data_TP = np.zeros(shape=(generations + 1, population))
# DR
gamma_data_DR = np.zeros(shape=(generations + 1, population))
a_data_DR = np.zeros(shape=(generations + 1, population))
m_data_DR = np.zeros(shape=(generations + 1, population))

gamma_data_TP[0, :] = params_TP[:, 0]
a_data_TP[0, :] = params_TP[:, 1]
nu_data_TP[0, :] = params_TP[:, 4]

gamma_data_DR[0, :] = params_DR[:, 0]
a_data_DR[0, :] = params_DR[:, 1]
m_data_DR[0, :] = params_DR[:, 3]

fitness = np.zeros(shape=(generations, total_population))
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

# weights for paras of dr
weight_gamma_dr = np.zeros(population)
weight_gamma_dr.fill(1 / population)
weight_a_dr = np.zeros(population)
weight_a_dr.fill(1 / population)
weight_m_dr = np.zeros(population)
weight_m_dr.fill(1 / population)
weight_del_dr = np.zeros(population)
weight_del_dr.fill(1 / population)

for g in range(generations):
    # model 0
    TP_sample_length = len(np.where(model_data[g, :] == 0)[0])
    DR_sample_length = len(np.where(model_data[g, :] == 1)[0])
    Z = np.zeros((1, td.total_species))

    if TP_sample_length > 0:
        print('TP simulations start...')

        simmodelTP = dvcpp.DVSimLiang(td, params_TP)
        valid_TP = np.where(simmodelTP['sim_time'] == td.sim_evo_time)[0]
        Z_modelTP = simmodelTP['Z'][valid_TP]
        i, j = argsort2D(Z_modelTP)
        Z_modelTP = Z_modelTP[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTP = np.nan_to_num(Z_modelTP)
        Z = np.vstack([Z, Z_modelTP])
        # V = np.nan_to_num(V)
        # GOF: Goodness of fit
    else:
        valid_TP = []

    if DR_sample_length > 0:
        print('DR simulations start...')
        # model 1
        # for param_drury in params_DR:
        simmodelDR = dvcpp.DVSimDrury(td, params_DR)
        valid_DR = np.where(simmodelDR['sim_time'] == td.sim_evo_time)[0]
        Z_modelDR = simmodelDR['Z'][valid_DR]
        i, j = argsort2D(Z_modelDR)
        Z_modelDR = Z_modelDR[i, j]
        Z = np.vstack([Z, Z_modelDR])
    else:
        valid_DR = []


    Z = Z[1:, ]

    valid = np.concatenate([valid_TP, np.array(valid_DR) + TP_sample_length
                            ]).astype(int)

    eudis = normalized_norm(Z, obsZ)
    # eudis = eudistance(Z, obsZ)

    fitness[g, valid] += 1 - eudis

    # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g, :])[-int(total_population // 2)]  # best 25%
    fit_index = np.where(fitness[g, :] > fitness[g, q5])[0]

    modelTPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
    modeldrperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)

    print('Iteration = %d 50th Model TP: %.1f%% ;  Model DR: %.1f%%...'
          % (g, modelTPperc * 100, modeldrperc * 100))
    print('Average fitness: %f' % np.mean(fitness[g, fit_index]))
    # reevaluate the weight of the best fitted  models
    weight_model_bestfitted = weight_model[fit_index] * fitness[g, fit_index] / sum(
        weight_model[fit_index] * fitness[g, fit_index])
    # sample new models from the fitness of previous best fitted models
    sample_model_index = sorted(np.random.choice(fit_index, total_population, p=weight_model_bestfitted))

    if allowmodeldie == 'on':
        propose_model = propose_model[sample_model_index]
        previous_bestfitted_index_TP = fit_index[np.where(fit_index < population)]
        previous_bestfitted_index_DR = fit_index[np.where(np.logical_and(fit_index >= population, \
                                                                         fit_index < 2 * population))] - \
                                       population


    else:
        propose_model = model_params

        q5_TP = np.argsort(fitness[g, :population - 1])[-int(population // 20)]  # best 25%
        q5_DR = np.argsort(fitness[g, population:2 * population - 1])[-int(population // 20)] + population  # best 25%

        fit_index_TP = np.where(fitness[g, :population - 1] > fitness[g, q5_TP])[0]
        fit_index_DR = np.where(fitness[g, population:2 * population - 1] > fitness[g, q5_DR])[0] + population

        previous_bestfitted_index_TP = fit_index_TP
        previous_bestfitted_index_DR = fit_index_DR - population

        chosengamma_TP, chosena_TP, chosennu_TP = np.mean(params_TP[previous_bestfitted_index_TP, 0]), \
                                                  np.mean(params_TP[previous_bestfitted_index_TP, 1]), \
                                                  np.mean(params_TP[previous_bestfitted_index_TP, 4])

        chosengamma_DR, chosena_DR, chosenm_DR = np.mean(params_DR[previous_bestfitted_index_DR, 0]), \
                                                 np.mean(params_DR[previous_bestfitted_index_DR, 1]), \
                                                 np.mean(params_DR[previous_bestfitted_index_DR, 3])


        print('Mean estimates: TP gamma: %f ; a: %f ; nu: %.3e' % (chosengamma_TP, chosena_TP, chosennu_TP))
        print('Mean estimates: DR gamma: %f ; a: %f ; m: %f' % (chosengamma_DR, chosena_DR, chosenm_DR))

    model_data[g + 1, :] = propose_model
    gamma_data_TP[g + 1, :] = params_TP[:, 0]
    a_data_TP[g + 1, :] = params_TP[:, 1]
    nu_data_TP[g + 1, :] = params_TP[:, 4]

    gamma_data_DR[g + 1, :] = params_DR[:, 0]
    a_data_DR[g + 1, :] = params_DR[:, 1]
    m_data_DR[g + 1, :] = params_DR[:, 3]

    if len(np.where(propose_model == 0)[0]) > 0:
        # update TP paras and weights
        weight_gamma_TP, weight_a_TP, weight_nu_TP, propose_gamma_TP, propose_a_TP, propose_nu_TP = \
            tp_update(previous_bestfitted_index_TP, propose_model, params_TP, weight_gamma_TP,
                      weight_a_TP, weight_nu_TP)
        modelTP = np.where(propose_model == 0)
        params_TP = np.tile(param_Liang, (len(modelTP[0]), 1))
        params_TP[:, 0] = propose_gamma_TP
        params_TP[:, 1] = propose_a_TP
        params_TP[:, 4] = propose_nu_TP

    if len(np.where(propose_model == 1)[0]) > 0:

        # update DR paras and weights
        weight_gamma_dr, weight_a_dr, weight_m_dr, weight_del_dr, propose_gamma_dr, propose_a_dr, propose_m_dr, propose_del_dr = \
            dr_update(previous_bestfitted_index_DR, propose_model, params_DR, weight_gamma_dr,
                      weight_a_dr, weight_m_dr, weight_del_dr)
        modelDR = np.where(propose_model == 1)
        params_DR = np.tile(param_Drury, (len(modelDR[0]), 1))
        params_DR[:, 0] = propose_gamma_dr
        params_DR[:, 1] = propose_a_dr
        if fix_m == 1:
            params_DR[:, 3] = 1.0  # propose_m_dr
        else:
            params_DR[:, 3] = propose_m_dr

        if del_mute == 'on':
            params_DR[:, 5] = sigma2
            weight_m_dr.fill(1 / len(weight_m_dr))

        else:
            params_DR[:, 5] = propose_del_dr
    #
#

para_data = {'model_data': model_data, 'fitness': fitness, 'gamma_data_TP': gamma_data_TP,
             'a_data_TP': a_data_TP, 'nu_data_TP': nu_data_TP, 'gamma_data_DR': gamma_data_DR,
             'a_data_DR': a_data_DR, 'm_data_DR': m_data_DR,'Z':Z,'valid':valid}


np.save(savedir, para_data)

# smc = np.load(file2).item()
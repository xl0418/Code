import sys
sys.path.append('C:/Liang/abcpp_ms/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang,DVParamDrury
import dvtraitsim_cpp as dvcpp
sys.path.append('C:/Liang/Code/Pro2/BaleenWhale_code/')
from dr_update import dr_update
import matplotlib.pyplot as plt
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err

# prep data
# dir_path = os.path.dirname(os.path.realpath(__file__))
# files = dir_path + '/../tree_data/example1/'
dir_path = 'C:\\Liang\\Googlebox\\Research\\Project2\\'
files = dir_path + 'treesim_newexp/example1/'
savedir = dir_path+'modelsele/'
file2 = savedir + 'example1/MSdr.npy'

td = DVTreeData(path=files, scalar=1000)
print(td.total_species)

meantrait=0.0
sigma2 = 0.001  # Brownian Motion variance
del_mute = 'on'
allowmodeldie = 'off'
fix_m =1
# the generating params for models


# DR
gene_gamma_dr = 0.00001
gene_a_dr=0.001
gene_m_dr=gene_a_dr


modelnum = 1
population = 10000
total_population = population * modelnum
generations = 30
# let's try to find a true simulation:
# The generating paras

gene_candiparam = DVParamDrury(gamma = 0.0, a=0.0, theta=meantrait, m=gene_m_dr, inittrait=meantrait, var_trait=sigma2)



gene_candiparam[ 0] = gene_gamma_dr  # randomize 'gamma'
gene_candiparam[ 1] = gene_a_dr  # randomize 'a'
gene_candiparam[ 3] = gene_m_dr  # randomize 'm'
gene_candiparam[ 5] = sigma2  # randomize delta
gene_modelDR = dvcpp.DVSimDrury(td, gene_candiparam)
s = np.argsort(gene_modelDR['Z'])
obsZ = gene_modelDR['Z'][s]


meantrait = np.mean(obsZ)

# The generating paras

print('With the generating params: gamma:%f; a:%f; m: %f...' % (gene_gamma_dr,gene_a_dr,gene_m_dr) )

# let's try to find a true simulation:

param_Drury = DVParamDrury(gamma = 0.0, a=0.0, theta=meantrait, m=1.0, inittrait=meantrait, var_trait=sigma2)
# pop = dvcpp.DVSim(td, obs_param)
fitness = np.zeros(shape=(generations, total_population))

params_DR = np.tile(param_Drury, (population, 1))
params_DR[:, 0] = np.random.uniform(0.0, 1.0, params_DR.shape[0])  # randomize 'gamma'
params_DR[:, 1] = np.random.uniform(0.0, 1.0, params_DR.shape[0])  # randomize 'a'
params_DR[:, 3] = np.random.uniform(0.0, 5.0, params_DR.shape[0])  # randomize 'm'
params_DR[:, 5] = sigma2  # randomize delta

propose_model = np.zeros(shape=(1, total_population))
propose_model.fill(1.0)
# store parameters used

# DR
gamma_data_DR = np.zeros(shape=(generations + 1, population))
a_data_DR = np.zeros(shape=(generations + 1, population))
m_data_DR = np.zeros(shape=(generations + 1, population))


gamma_data_DR[0, :] = params_DR[:, 0]
a_data_DR[0, :] = params_DR[:, 1]
m_data_DR[0, :] = params_DR[:, 3]


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
    print('DR simulations start...')
    # model 1
    # for param_drury in params_DR:
    simmodelDR = dvcpp.DVSimDrury(td, params_DR)
    valid_DR = np.where(simmodelDR['sim_time'] == td.sim_evo_time)[0]
    Z_modelDR = simmodelDR['Z'][valid_DR]
    i, j = argsort2D(Z_modelDR)
    Z_modelDR = Z_modelDR[i, j]

    valid = valid_DR
    eudis = normalized_norm(Z_modelDR, obsZ)
    # eudis = eudistance(Z, obsZ)

    fitness[g, valid] += 1 - eudis

    # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g, :])[-int(total_population // 100)]  # best 25%
    fit_index = np.where(fitness[g, :] > fitness[g, q5])[0]

    # sample new models from the fitness of previous best fitted models

    q5_DR = np.argsort(fitness[g, :])[-int(population // 100)] # best 25%

    fit_index_DR = np.where(fitness[g, :] > fitness[g, q5_DR])[0]

    previous_bestfitted_index_DR = fit_index_DR

    chosengamma_DR, chosena_DR, chosenm_DR = np.mean(params_DR[previous_bestfitted_index_DR, 0]), \
                                             np.mean(params_DR[previous_bestfitted_index_DR, 1]), \
                                             np.mean(params_DR[previous_bestfitted_index_DR, 3])


    print('%i th iteration, mean estimates: DR gamma: %f ; a: %f ; m: %f' % (g,chosengamma_DR, chosena_DR, chosenm_DR))

    gamma_data_DR[g + 1, :] = params_DR[:, 0]
    a_data_DR[g + 1, :] = params_DR[:, 1]
    m_data_DR[g + 1, :] = params_DR[:, 3]

    # update DR paras and weights
    weight_gamma_dr, weight_a_dr, weight_m_dr, weight_del_dr, propose_gamma_dr, propose_a_dr, propose_m_dr, propose_del_dr = \
        dr_update(previous_bestfitted_index_DR, propose_model, params_DR, weight_gamma_dr,
                  weight_a_dr, weight_m_dr, weight_del_dr)
    modelDR = np.where(propose_model == 1)
    params_DR = np.tile(param_Drury, (len(modelDR[0]), 1))
    params_DR[:, 0] = propose_gamma_dr
    params_DR[:, 1] = propose_a_dr
    if fix_m == 1.0:
        params_DR[:, 3] = propose_a_dr  # propose_m_dr
    else:
        params_DR[:, 3] = propose_m_dr

    if del_mute == 'on':
        params_DR[:, 5] = sigma2
        weight_m_dr.fill(1 / len(weight_m_dr))
    else:
        params_DR[:, 5] = propose_del_dr
    #
#

para_data = {'fitness': fitness,  'gamma_data_DR': gamma_data_DR,
             'a_data_DR': a_data_DR, 'm_data_DR': m_data_DR,'Z':Z_modelDR,'valid':valid}


# np.save(savedir, para_data)
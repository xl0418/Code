import os
import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm

gamma = 0
a = 0.001

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
dir_path = 'c:/Liang/Code/Pro2/abcpp'
files = dir_path + '/tree_data/example1/'

td = DVTreeData(path=files, scalar=10000)

prior = [0.5, 0.5, 0.5, 0.5]
gamma_prior_mean = prior[0]
gamma_prior_var = prior[1]
a_prior_mean = prior[2]
a_prior_var = prior[3]
K=1000000000
nu=0
# let's try to find a true simulation:
obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                    split_stddev=0.2)
print('try to find a completed true simulation with gamma =', obs_param[0], 'and a =', obs_param[1], '...')
for r in range(100):
    obs = dvcpp.DVSim(td, obs_param)
    if obs['sim_time'] == td.evo_time:
        break
if obs['sim_time'] < td.evo_time:
    print('hopeless, does not compute.')
    sys.exit(-1)
s = np.argsort(obs['Z'])
obsN = obs['N'][s]
obsZ = obs['Z'][s]
obsV = obs['V'][s]
obsN = np.nan_to_num(obsN)
obsZ = np.nan_to_num(obsZ)
obsV = np.nan_to_num(obsV)


population = 1000
generations = 10
params = np.tile(obs_param, (population, 1))  # duplicate
params[:, 0] = np.random.uniform(0.0, 1.0, params.shape[0])  # randomize 'gamma'
params[:, 1] = np.random.uniform(0.0, 1.0, params.shape[0])  # randomize 'a'
gamma_data = np.zeros(shape=(generations, population))
a_data = np.zeros(shape=(generations, population))
fitness= np.zeros(shape=(generations, population))
# Initialize the weights.
weight_gamma = np.zeros(population)
weight_gamma.fill(1 / population)
weight_a = np.zeros(population)
weight_a.fill(1 / population)
for g in range(generations):
    gamma_data[g, :] = params[:, 0]
    a_data[g, :] = params[:, 1]
    pop = dvcpp.DVSim(td, params)

    # access fitness
    # fitness = np.zeros(population)
    valid = np.where(pop['sim_time'] == td.evo_time)[0]
    if len(valid)<20:
        print("WARNING:Valid simulations are too scarce!")
    if valid.size > 0:
        Z = pop['Z'][valid]
        i, j = argsort2D(Z)
        Z = Z[i, j]
        V = pop['V'][valid][i, j]
        Z = np.nan_to_num(Z)
        V = np.nan_to_num(V)
        #GOF: Goodness of fit
        fitness[g,valid] += 1.0 - normalized_norm(Z, obsZ)
        fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g,:])[-population// 20]  # best 5%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]

    print('Iteration = %d 5th gamma = %f  a = %f  fitness = %f' % (g, np.mean(params[fit_index, 0]),
                                                                 np.mean(params[fit_index, 1]), np.mean(fitness[g,fit_index])))
    print('Iteration = %d all gamma = %f  a = %f  fitness = %f' % (g, np.mean(params[:, 0]),
                                                                 np.mean(params[:, 1]),
                                                                 np.mean(fitness[g,:])))

    weight_gamma = weight_gamma[fit_index]/sum(weight_gamma[fit_index])
    weight_a = weight_a[fit_index]/sum(weight_a[fit_index])

    gamma_pre_mean = np.sum(params[fit_index, 0] * weight_gamma)
    gamma_pre_var = np.sum((params[fit_index, 0] - gamma_pre_mean) ** 2 * weight_gamma)
    a_pre_mean = np.sum(params[fit_index, 1] * weight_a)
    a_pre_var = np.sum((params[fit_index, 1] - a_pre_mean) ** 2 * weight_a)

    # sample parameters by the weights computed in last loop.
    sample_gamma_index = np.random.choice(fit_index, population, p=weight_gamma)
    sample_a_index = np.random.choice(fit_index, population, p=weight_a)
    # mean of the sample for gamma
    propose_gamma0 = params[sample_gamma_index, 0]
    # draw new gamma with mean and variance
    propose_gamma = abs(np.random.normal(propose_gamma0, np.sqrt(2 * gamma_pre_var)))
    # mean of the sample for a
    propose_a0 = params[sample_a_index, 1]
    # draw new a with mean and variance
    propose_a = abs(np.random.normal(propose_a0, np.sqrt(2 * a_pre_var)))

    extend_weight_gamma = weight_gamma[fit_index.searchsorted(sample_gamma_index)]
    extend_weight_a = weight_a[fit_index.searchsorted(sample_a_index)]

    # compute new weights for gamma and a
    weight_gamma_denominator = np.sum(extend_weight_gamma * norm.pdf(propose_gamma, params[:, 0],
                                                                     np.sqrt(2 * gamma_pre_var)))
    weight_gamma_numerator = norm.pdf(propose_gamma, gamma_prior_mean, gamma_prior_var)
    weight_gamma = weight_gamma_numerator / weight_gamma_denominator

    weight_a_denominator = np.sum(extend_weight_a * norm.pdf(propose_a, params[:, 1],
                                                             np.sqrt(2 * a_pre_var)))
    weight_a_numerator = norm.pdf(propose_a, a_prior_mean, a_prior_var)
    weight_a = weight_a_numerator / weight_a_denominator
    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma = weight_gamma / sum(weight_gamma)
    weight_a = weight_a / sum(weight_a)
    params[:, 0] = propose_gamma
    params[:, 1] = propose_a
#
para_data = {'gamma': gamma_data, 'a': a_data,'fitness': fitness}
# file='C:/Liang/Code/Pro2/abcpp/abcpp/smcndata/'
# # file = '/home/p274981/abcpp/abcpp/'
# file2 = file + 'paraT11.npy'
# np.save(file2,para_data)
#
# smc = np.load(file2).item()
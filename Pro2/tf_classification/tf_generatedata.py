import os
import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
from scipy.stats import norm


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
dir_path = 'c:/Liang/Code/Pro2/abcpp'
files = dir_path + '/tree_data/example1/'

td = DVTreeData(path=files, scalar=10000)

# prior = [0.5, 0.5, 0.5, 0.5,1e-13,1e-12]
# gamma_prior_mean = prior[0]
# gamma_prior_var = prior[1]
# a_prior_mean = prior[2]
# a_prior_var = prior[3]
# nu_prior_mean = prior[4]
# nu_prior_var = prior[5]

gamma = 0.001
a_vec = [0.0,0.001,0.01,0.1,0.5,1]
# a_labels = [i for i in range(len(a_vec))]
Z_train = []
label_train = []
K=10e8
nu=1/(100*K)
# let's try to find a true simulation:
datasize_batch = 1000

for count in range(len(a_vec)):
    a = a_vec[count]

    obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                        initpop_sigma=10.0, break_on_mu=False)
    print('try to find a completed true simulation with gamma =', obs_param[0], 'and a =', obs_param[1], 'and nu =', obs_param[3],'...')
    for r in range(datasize_batch):
        print(r)
        obs = dvcpp.DVSim(td, obs_param)
        if obs['sim_time'] == td.sim_evo_time:
            s = np.argsort(obs['Z'])
            obsZ = obs['Z'][s]
            Z_train.append(list(obsZ))
            label_train = np.concatenate((label_train,[count]))
    if obs['sim_time'] < td.evo_time:
        print('hopeless, does not compute.')
        sys.exit(-1)

Z_train_ndarray = np.array(Z_train)

#
para_data = {'Z_train': Z_train_ndarray, 'Z_labels': label_train}
file='C:/Liang/Code/Pro2/tf_classification/'
# # file = '/home/p274981/abcpp/abcpp/'
filename = file + 'tf_c1_batch2.npy'
np.save(filename,para_data)
#
# smc = np.load(filename).item()
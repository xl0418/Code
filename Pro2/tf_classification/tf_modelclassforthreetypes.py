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
files = dir_path + '/tree_data/example2/'

td = DVTreeData(path=files, scalar=10000)

Z_train = []
label_train = []
K=10e8
nu=1/(100*K)
# let's try to find a true simulation:
datasize_batch = 200000
mode = ['BM','OU','TP','antiOU']
for count in range(len(mode)):
    model = mode[count]
    gamma = 0
    a = 0

    obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                        initpop_sigma=10.0, break_on_mu=False)

    params = np.tile(obs_param, (datasize_batch, 1))  # duplicate

    if model == 'OU':
        params[:, 0] = np.random.uniform(1e-10, 1e-2, params.shape[0])  # randomize 'gamma'
    elif model == 'TP':
        params[:, 0] = np.random.uniform(1e-10, 1e-2, params.shape[0])  # randomize 'gamma'
        params[:, 1] = np.random.uniform(1e-10, 1.0, params.shape[0])  # randomize 'a'
    elif model == 'antiOU':
        params[:, 1] = np.random.uniform(1e-10, 1.0, params.shape[0])  # randomize 'a'

    print('Batch ', count,': try to generate ', datasize_batch,' simulations with gamma =', obs_param[0], 'and a =', obs_param[1], 'and nu =', obs_param[3],'...')

    pop = dvcpp.DVSim(td, params)
    valid = np.where(pop['sim_time'] == td.sim_evo_time)[0]
    if len(valid) < 20:
        print("WARNING:Valid simulations are too scarce!")
    if valid.size > 0:
        Z = pop['Z'][valid]
        i, j = argsort2D(Z)
        Z = Z[i, j]
        # V = pop['V'][valid][i, j]
        Z = np.nan_to_num(Z)
        label_train_batch = np.zeros(shape=(1,Z.shape[0]))
        label_train_batch.fill(count)
        label_train = np.concatenate((label_train,label_train_batch),axis=None)
        if count == 0:
            Z_train = Z
        else:
            Z_train = np.concatenate((Z_train,Z))

        para_data = {'Z_train': Z_train, 'Z_labels': label_train}
        file = 'C:/Liang/Code/Pro2/tf_classification/%s' % model
        # # file = '/home/p274981/abcpp/abcpp/'
        filename = file + 'tf_tree2modelseletrain2b.npy'
        np.save(filename, para_data)
#

#
# smc = np.load(filename).item()
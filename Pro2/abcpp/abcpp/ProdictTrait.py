import os,sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import pandas as pd

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'
gamma = .000236
a = .617900
nu = 2.411e-11
K = 10e8
meantrait = 32.50540571860487

td = DVTreeData(path=files, scalar=1000)
param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)
params = np.tile(param, (1000, 1))       # duplicate


predictsim = dvcpp.DVSim(td, params)

valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]

Z = predictsim['Z'][valid]
i, j = argsort2D(Z)
Z = Z[i, j]
# V = pop['V'][valid][i, j]
Z = np.nan_to_num(Z)

Z_df = pd.DataFrame(Z)
Z_df.to_csv(files+'predictsim.csv',sep=',',index=False)





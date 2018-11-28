#%%
import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import numpy as np
import pandas as pd
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


theta = 0.0  # optimum of natural selection
r = 1.0  # growth rate
Vmax = 1.0
scalar = 1000
K=10e8
nu=1/(100*K)

# let's try to find a true simulation:


full =0

#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/'
if full==1:
    files = dir_path + 'full_data/'
else:
    files = dir_path + 'pruend_data/'


# trait evolution plot
td = DVTreeData(path=files, scalar=scalar)
gamma1 = 0.0001
a = 0.5

replicate = 1000
obs_param = DVParam(gamma=gamma1, a=a, K=K, nu=nu, r=r, theta=theta, Vmax=1, inittrait=0, initpop=500,
                    initpop_sigma=10.0, break_on_mu=False)
params = np.tile(obs_param, (replicate, 1))
simresult = dvcpp.DVSim(td, params)

obsZ_unsort = simresult['Z']

i, j = argsort2D(obsZ_unsort)
obsZ_sort = simresult['Z'][i,j]


Z_unsort_df = pd.DataFrame(obsZ_unsort)
Z_sort_df = pd.DataFrame(obsZ_sort)


fileunsort_name = os.getcwd()+ '\\Pro2\\data\\unsort.csv'
filesort_name = os.getcwd()+ '\\Pro2\\data\\sort.csv'

Z_unsort_df.to_csv(fileunsort_name, encoding='utf-8', index=False)
Z_sort_df.to_csv(filesort_name, encoding='utf-8', index=False)



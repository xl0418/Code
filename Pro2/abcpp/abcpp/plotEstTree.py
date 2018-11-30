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
import csv

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


with open(dir_path+'diameter.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    diameterdata = list(csv_reader)

ddarray = np.array(diameterdata)
tipslabel = ddarray[1:,0]
obsZ = ddarray[1:,1]
obsV = ddarray[1:,2]
obsN = ddarray[1:,3]
if full==1:
    missingdata = np.array([ 57, 160, 165, 163, 182])-1  # missing index in full data
else:
    missingdata = np.array([ 1, 15, 19, 18, 29])-1  # missing index in reconstructive data

naindex = np.where(obsZ == 'NA')[0]
obsZ_delna = np.delete(obsZ,naindex)
remove_nalabel = np.delete(tipslabel,naindex)
s = np.argsort(obsZ_delna)
obsZ_delna = obsZ_delna[s]
obsZ_delna = obsZ_delna.astype(np.float)
sorted_label = remove_nalabel[s]

sorted_label = np.concatenate([sorted_label,tipslabel[naindex]])
obsZ = np.concatenate([obsZ_delna,np.array([-1,-1,-1,-1,-1])])
emplist = {'specieslabel':sorted_label, 'trait':obsZ}
empiricaldata_sorted = pd.DataFrame(emplist)

fileemp_name = os.getcwd()+ '\\Pro2\\data\\emp.csv'

empiricaldata_sorted.to_csv(fileemp_name, encoding='utf-8', index=False)


# simulate trait data
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


# Data matching
nor_obsZ_delna = (obsZ_delna-np.min(obsZ_delna))/(np.max(obsZ_delna)-np.min(obsZ_delna))


mean_simZ = np.mean(obsZ_sort, axis=0)
nor_mean_simZ = (mean_simZ-np.min(mean_simZ))/(np.max(mean_simZ)-np.min(mean_simZ))


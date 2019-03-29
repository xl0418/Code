import sys
import platform
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from Candipararun import Candi_para

if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam
from multiprocessing import Pool


theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 1000
K=10e8
nu=1/(100*K)
timegap = 1
meantrait = 0.0
population = 100
sigma2 = 0.5
# let's try to find a true simulation:
obs_param = DVParam(gamma=0.01, a=0.5, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500000,
             initpop_sigma = 10.0, break_on_mu=False)
candiparam = np.array([0.0, 0.0, meantrait, 1.0, meantrait, 1.0])
# pop = dvcpp.DVSim(td, obs_param)


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)

params_DR =  np.tile(candiparam,(population,1))
params_DR[:,0]= np.random.uniform(0.0, 1.0, params_DR.shape[0])         # randomize 'gamma'
params_DR[:,1]= np.random.uniform(0.0, 1.0, params_DR.shape[0])         # randomize 'a'
params_DR[:,3]= np.random.uniform(0.0, 5.0, params_DR.shape[0])         # randomize 'm'
params_DR[:, 5] = np.random.uniform(0.0, sigma2, params_DR.shape[0])  # randomize delta

params_nh =  np.tile(candiparam,(population,1))
params_nh[:,0]= np.random.uniform(0.0, 1.0, params_nh.shape[0])     # randomize 'gamma'
params_nh[:,1]= 0
params_nh[:,3]= np.random.uniform(0.0, 5.0, params_nh.shape[0])     # randomize 'm'
params_nh[:, 5] = np.random.uniform(0.0, sigma2, params_nh.shape[0])  # randomize delta

p = Pool(8)
res_DR = p.map( Candi_para, (params_DR[i,:] for i in range(population)) )

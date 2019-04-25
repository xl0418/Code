import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms2/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_py import DVSimLiang
from dvtraitsim_shared import DVTreeData, DVParamLiang


theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 20000
K=10e8
nu=1/(100*K)
# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'


td = DVTreeData(path=files, scalar=scalar)

# gamma, a, K, h, nu, r, theta, V00, V01, Vmax, inittrait, initpop, initpop_sigma, break_on_mu
obs_param = DVParamLiang(gamma=0.0001, a=0.1, K=K,h=1, nu=nu, r=r, theta=theta,V00=0.1,V01=0.1, Vmax=1, inittrait=0, initpop=500,
                    initpop_sigma=10.0, break_on_mu=False)
simresult = DVSimLiang(td,obs_param)

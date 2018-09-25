import os
from time import perf_counter
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import dvtraitsim_py as dvpy


dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/example1/'
td = DVTreeData(path=files, scalar=10000)

param = DVParam(gamma=0.001, a=1, K=10000, nu=0.000001, r=1, theta=0, Vmax=1, inittrait=0, initpop=500, split_stddev=0.2, keep_alive=0.5)
params = np.tile(param, (1000, 1))       # duplicate

print('running', params.shape[0], 'simulations (parallel):')
start = perf_counter()
dispatch = dvcpp.DVSim(td, params)
time = perf_counter() - start
print("valid simulations:", np.where(dispatch['sim_time'] == td.evo_time)[0].size)
print("in", time, "s")

print('running', 1, 'simulation (python):')
start = perf_counter()
R = dvpy.DVSim(td, param)
print("in", perf_counter() - start, "s")

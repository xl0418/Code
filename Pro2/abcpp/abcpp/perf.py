import os
from time import perf_counter
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import dvtraitsim_py as dvpy


dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/example11/'
td = DVTreeData(path=files, scalar=10000)

param = DVParam(gamma=0.001, a=1, K=10000, nu=0.0001, r=1, theta=0, Vmax=1, inittrait=0, initpop=500, split_stddev=0.2)
params = np.tile(param, (100, 1))       # duplicate

print('running ', params.shape[0], ' simulations (cpp):')
start = perf_counter()
dispatch = dvcpp.DVSim(td, params)
time = perf_counter() - start
for x in dispatch['sim_time']:
    if x < td.evo_time: print('inconsistency at time ', x)
print("in", time, "s")

print('running ', params.shape[0], ' simulations (python):')
start = perf_counter()
R = list()
for p in params:
    R.append(dvpy.DVSim(td, p))
print("in", perf_counter() - start, "s")

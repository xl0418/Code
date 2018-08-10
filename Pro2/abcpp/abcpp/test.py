import os
from time import perf_counter
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp


dir_path = os.path.dirname(os.path.realpath(__file__))
files = dir_path + '/../tree_data/example1/'

td = DVTreeData(path=files, scalar=1000)
param = DVParam(gamma=0.001, a=1, K=10000, nu=0.0001, r=1, theta=0, Vmax=1, inittrait=0, initpop=500, split_stddev=0.2)
params = np.tile(param, (1000, 1))       # duplicate

print('running cvtraits_cpp.DVSim on: ', files)
start = perf_counter()
dispatch = dvcpp.DVSim(td, param)
print("in ", perf_counter() - start, "s")

fitnes = np.random.uniform(1, 1000, 1000)
start = perf_counter()
dvcpp.discrete_distribution(fitnes, 1000)
print("in ", perf_counter() - start, "s")



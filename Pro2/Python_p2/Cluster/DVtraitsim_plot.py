from DV_model_sim_along_phy1 import DVtraitsim_tree
import os
import numpy as np


theta = 0  # optimum of natural selection
r = 1  # growth rate
K = 5000  # carrying capacity
delta_pop = .001  # Variance of random walk of population
nu = 0.0001
Vmax = 1
scalor = 1000
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec
count1 = 0
count2 = 0
# trait evolution plot
file = '/home/p274981/Python_p2/DVmodel/'
replicate = 1000
for gamma1 in gamma_vec:
    count1 += 1
    for a in a_vec:
        count2 += 1
        directory = '/data/p274981/Project2/%dgamma%da' % (count1,count2)
        if not os.path.exists(directory):
            os.makedirs(directory)
        for loop in range(0,replicate+1):
            simresult = DVtraitsim_tree(file = file, replicate = loop, gamma1 = gamma1, a = a)
            traitfilestr = directory + "/%dgamma%da%dreplicatetrait.txt" % (count1,count2,loop)
            popfilestr = directory + "/%dgamma%da%dreplicatetraitpop.txt" % (count1,count2,loop)
            np.savetxt(traitfilestr,simresult[0])
            np.savetxt(popfilestr,simresult[1])

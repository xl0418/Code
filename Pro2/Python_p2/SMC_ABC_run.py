import sys, os
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from ABC_MCMC import SMC_ABC, single_trait_sim
import numpy as np
from Trait_sim_in_branches_stat import traitsim
import scipy.stats
from scipy.stats import norm

# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=0)

epsilon = 30
timestep = 8
particlesize = 100
prior = [0.3,0.5,0.2,0.5]
SMC_ABC_result = SMC_ABC(timestep = timestep, particlesize = particlesize, obs = obs, prior = prior, epsilon = epsilon)
np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/SMC_ABC/SMC_ABC_result.txt",SMC_ABC_result)
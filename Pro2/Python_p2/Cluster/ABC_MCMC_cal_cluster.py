import numpy as np
from ABC_SMC_DVmodel1 import SMC_ABC, single_trait_sim,calibrication, MCMC_ABC
from DV_model_sim_along_phy1 import DVtraitsim_tree

# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

# Observation generated
#  load the data for a given tree
file = '/home/p274981/Python_p2/DVmodel/'
simresult = DVtraitsim_tree(file = file, replicate = 1, gamma1 = par_obs[0], a = par_obs[1],scalor = 1000)
evo_time, total_species = simresult[0].shape
evo_time = evo_time-1
trait_RI_dr = simresult[0]
population_RI_dr = simresult[1]
# empirical data for trait and population
trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]
# observation data
obs = np.array([trait_dr_tips,population_tips])


# Calibrication step
cal_size = 40000
# TEST1: Uniform prior distribution example
priorpar = [0.0001,1,0.0001,1]
collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, file = file)
file1 = file + 'cal4w_DVseed13.txt'
np.savetxt(file1,collection)

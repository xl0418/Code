import numpy as np
from ABC_SMC_DVmodel1 import SMC_ABC, single_trait_sim,calibration, MCMC_ABC
from DV_model_sim_along_phy1 import DVtraitsim_tree

# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

# Observation generated
#  load the data for a given tree
file ='/home/p274981/Python_p2/DVmodel/'
simresult = DVtraitsim_tree(file = file,replicate = 13, gamma1 = par_obs[0], a = par_obs[1],scalar = 1000)
evo_time, total_species = simresult[0].shape
evo_time = evo_time-1
trait_RI_dr = simresult[0]
population_RI_dr = simresult[1]
traitvar = simresult[3]
# empirical data for trait and population
trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]
traitvar = traitvar[evo_time,:][~np.isnan(traitvar[evo_time,:])]

# observation data
obs = np.array([trait_dr_tips,population_tips,traitvar])


# Calibrication step 1: rejection process based on trait mean with uniform prior.
calimean_file= '/home/p274981/Python_p2/DVmodel/calimean'
cal_size = 30000
priorpar_trait = [0.0001,1,0.0001,1]
collection_trait = calibration(samplesize = cal_size, priorpar = priorpar_trait,  treefile = file,calidata_file =calimean_file,  calmode = 'uni')
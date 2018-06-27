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
file1 = file + 'cal10w_DV.txt'

collection = np.loadtxt(file1)

threshold = 0.05
num = threshold*cal_size-1
delta = np.sort(collection[:,2])[int(num)]
mn,idx = min( (collection[i,2],i) for i in range(len(collection[:,2])) )
startvalue_par = collection[idx,:2]

filtered_coll = collection[collection[:,2]<=delta]

# ABC_MCMC step
iterations = 100000

posterior = MCMC_ABC(startvalue= startvalue_par, iterations = iterations, delta = delta, obs = obs,sort = 1,
                    priorpar=priorpar, file = file, mode = 'nor')
file2 = file + 'posterior10w_DV.txt'
np.savetxt(file2,posterior)



#Import desired functions
import numpy as np
from ABC_SMC_DVmodel import single_trait_sim,calibration, MCMC_ABC
from DV_model_sim_along_phy import DVtraitsim_tree

# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

# Observation generated
# Load the data for a given tree
# The directory of the tree data
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
# Simulate data
simresult = DVtraitsim_tree(file = file,replicate = 3, gamma1 = par_obs[0], a = par_obs[1],scalar = 1000)
# We only need the data at tips.
evo_time, total_species = simresult[0].shape
evo_time = evo_time-1
trait_RI_dr = simresult[0]
population_RI_dr = simresult[1]
traitvar = simresult[3]
trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]
traitvar = traitvar[evo_time,:][~np.isnan(traitvar[evo_time,:])]
# observation data
obs = np.array([trait_dr_tips,population_tips,traitvar])


# Calibrication step
cal_size = 100
# TEST1: Uniform prior distribution example
priorpar = [0.0001,1,0.0001,1]
collection = calibration(samplesize = cal_size, priorpar = priorpar, obs = obs, file = file)
# file1 = file + 'cal10w_DV.txt'
# np.savetxt(file1,collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/DVmodel/cal4w_DVseed13.txt")


# Processing the calibration data. We would like to chose the 5% closest simulation data to the observation
# to provide the prior information for MCMC step.
threshold = 0.01
num = threshold*cal_size-1
# The minimum of the distance of sorted data
delta = np.sort(collection[:,3])[int(num)]
mn,idx = min( (collection[i,3],i) for i in range(len(collection[:,3])) )
# Start value regarding to the minimum distance.
startvalue_par = collection[idx,:2]
# Filter the calibration data by 5% closest threshold.
filtered_coll = collection[collection[:,2]<=delta]

# ABC_MCMC step
iterations = 100
posterior = MCMC_ABC(startvalue= startvalue_par, iterations = iterations, delta = delta, obs = obs,sort = 1,
                     priorpar=priorpar, file = file, mode = 'nor')
file2 = file + 'posterior10w_DV.txt'
np.savetxt(file2,posterior)
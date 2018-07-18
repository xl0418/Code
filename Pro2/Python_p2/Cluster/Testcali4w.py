import numpy as np
from testABC_sim import DVtraitsim_tree
from testABCfunc import single_trait_sim,calibration, MCMC_ABC

# Observation parameters [gamma,a]
par_obs = np.array([0.01,0.5])

# Observation generated
#  load the data for a given tree
file = '/home/p274981/Python_p2/DVmodel/'
simresult = DVtraitsim_tree(file = file, replicate = 1, gamma1 = par_obs[0], a = par_obs[1],scalar = 1000)
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
calidata_file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\ABCtestcali4w'
collection = calibration(samplesize = cal_size, priorpar = priorpar, treefile = file, calidata_file = calidata_file)

cali1 = np.load(calidata_file+'.npz')
calitrait1 = cali1['calitrait']
calipop1 = cali1['calipop']
calivar1 = cali1['calivar']
calipara1 = cali1['calipar']

coll1 = np.empty(shape=(cal_size, 6))
sort_obstrait = np.sort(obs[0])
sort_obsvar = np.sort(obs[2])

for i in range(0,cal_size):
    meandiff_trait = np.linalg.norm(calitrait1[i] - obs[0])
    meandiff_trait_sort = np.linalg.norm(np.sort(calitrait1[i]) - sort_obstrait)
    meandiff_var = np.linalg.norm(calivar1[i] - obs[2])
    meandiff_var_sort = np.linalg.norm(np.sort(calivar1[i]) - sort_obsvar)
    coll1[i] = np.concatenate((calipara1[i],[meandiff_trait],[meandiff_trait_sort],[meandiff_var],[meandiff_var_sort]))

# Data filtered by trait mean
threshold = 0.01
num = threshold*cal_size-1
# ln = 2: unsorted distance; ln = 3: sorted distance.
ln=3
delta = np.sort(coll1[:,ln])[int(num)]
mn,idx = min( (coll1[i,ln],i) for i in range(len(coll1[:,ln])) )
startvalue_par = coll1[idx,:2]

filtered_coll1 = coll1[coll1[:,ln]<=delta]
priorpar_var = [np.mean(filtered_coll1[:,0]),np.std(filtered_coll1[:,0]),np.mean(filtered_coll1[:,1]),np.std(filtered_coll1[:,1])]
priorpar_var[3] = 1
# ABC_MCMC step
iterations =100000
posterior = MCMC_ABC(startvalue= startvalue_par, iterations = iterations, delta = delta, obs = obs,sort = 1,
                     priorpar=priorpar_var, file = file, mcmcmode = 'nor')
file2 = file + 'testABCMCMC.txt'
np.savetxt(file2,posterior)


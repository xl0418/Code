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

calimean_file= '/home/p274981/Python_p2/DVmodel/calimean'
cal_size = 30000

cali1 = np.load(calimean_file+'.npz')
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
threshold = 0.05
num = threshold*cal_size-1
# ln = 2: unsorted distance; ln = 3: sorted distance.
ln=3
delta = np.sort(coll1[:,ln])[int(num)]
mn,idx = min( (coll1[i,ln],i) for i in range(len(coll1[:,ln])) )
startvalue_par = coll1[idx,:2]

filtered_coll1 = coll1[coll1[:,ln]<=delta]
# Calibrication step 1: rejection process based on trait variance with normal prior.
# Prior information of parameters through the first filter.
priorpar_var = [np.mean(filtered_coll1[:,0]),np.std(filtered_coll1[:,0]),np.mean(filtered_coll1[:,1]),np.std(filtered_coll1[:,1])]
calivar_file= 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\calivar'
collection_var = calibration(samplesize = cal_size, priorpar = priorpar_var, treefile = file,calidata_file =calivar_file, calmode = 'nor')
# cali2 = np.load(calivar_file+'.npz')
# calitrait2 = cali2['calitrait']
# calipop2 = cali2['calipop']
# calivar2 = cali2['calivar']
# calipara2 = cali2['calipar']
#
#
# coll2 = np.empty(shape=(cal_size, 6))
#
# for i in range(0,cal_size):
#     meandiff_trait = np.linalg.norm(calitrait2[i] - obs[0])
#     meandiff_trait_sort = np.linalg.norm(np.sort(calitrait2[i]) - sort_obstrait)
#     meandiff_var = np.linalg.norm(calivar2[i] - obs[2])
#     meandiff_var_sort = np.linalg.norm(np.sort(calivar2[i]) - sort_obsvar)
#     coll2[i] = np.concatenate(
#         (calipara2[i], [meandiff_trait], [meandiff_trait_sort], [meandiff_var], [meandiff_var_sort]))
#
# # Data filtered by trait mean
# threshold = 0.1
# num = threshold*cal_size-1
# # ln = 2: unsorted distance; ln = 3: sorted distance.
# ln=3
# delta = np.sort(coll2[:,ln])[int(num)]
# mn,idx = min( (coll2[i,ln],i) for i in range(len(coll2[:,ln])) )
# startvalue_par = coll2[idx,:2]
#
# filtered_coll2 = coll2[coll2[:,ln]<=delta]

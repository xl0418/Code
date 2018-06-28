from ABC_SMC_DVmodel1 import calibrication
from DV_model_sim_along_phy1 import DVtraitsim_tree
import numpy as np

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

# self study
cal_size = 1000
# Uniform prior distribution example
# priorpar = [0.1,0.1]
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = np.array([0,0.001,0.01,0.1,0.5,1])
priorpar = np.zeros(2)
i= 5
for i in range(0,5):
    for j in range(0,5):
       priorpar[0]=gamma_vec[i]
       priorpar[1] = a_vec[j]
       collection = calibrication(file= file, samplesize = cal_size, priorpar = priorpar, obs = obs,mode='self')
      # str = 'c:/Liang/Googlebox/Research/Project2/python_p2/selfcal%dgamma%da.txt' % (i,j)
       str = '/home/p274981/Python_p2/DVselfcal%dgamma%da.txt' % (i,j)
       np.savetxt(str,collection)

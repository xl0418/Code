import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from ABC_SMC_DVmodel import SMC_ABC, single_trait_sim,calibrication, MCMC_ABC
from DV_model_sim_along_phy import DVtraitsim_tree
import pylab as P
import matplotlib.mlab as mlab
from sklearn.neighbors import KernelDensity
import pymc3 as pm
# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

# Observation generated
#  load the data for a given tree
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
simresult = DVtraitsim_tree(file = file, gamma1 = par_obs[0], a = par_obs[1],scalor = 1000)
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
file1 = file + 'cal10w_DV.txt'
np.savetxt(file1,collection)
collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/DVmodel/DVmodelcalresult4w/cal4w_DV.txt")


#TEST2: Normal prior distribution example
# priorpar = [0.1,0.2,0.1,0.3]
# collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, mode = 'nor')
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/priorresult/calibration2w.txt")

#TEST3: Normal prior distribution with 3 MCMCs
# priorpar = [0.2,0.5,0.1,0.4]
# collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, mode = 'nor')
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/MCMC3/calibration2w_3chains.txt")

collection = filtered_coll

# distance distribution
P.figure()
dis_data = collection[:,[2,3]]
n, bins, patches = P.hist(dis_data, 15, density=1, histtype='bar',
                            color=['crimson', 'burlywood'],
                            label=['distance', 'sorted distance'])
P.legend()
plt.show()
#Estimate prior distribution of parameters
# Generate random samples from a mixture of 2 Gaussians
# with modes at 5 and 10
data = np.array(collection[:,0])
data = data.reshape(-1,1)
# Plot the true distribution
x = np.linspace(0, 1, 100)[:, np.newaxis]
norm_vals = mlab.normpdf(x, 0.1, 0.4)
plt.plot(x, norm_vals)
# Plot the data using a normalized histogram
plt.hist(data, 50, density=True)
# Do kernel density estimation
kd = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(data)
# Plot the estimated densty
kd_vals = np.exp(kd.score_samples(x))
plt.plot(x, kd_vals,'--r')
# Show the plots
plt.show()



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
# posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/001result10w/posterior.txt")
# pm.autocorrplot(posterior)

#
# # Statistic
# posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/Normaldistributionresult/posterior_nor.txt")
# # Distribution plots for parameters
# gamma_samples = posterior[0::,0]
# a_samples = posterior[0::, 1]
# figdis = plt.figure(figsize=(12, 8))
#
# iterations = 20000
# plt.subplot(211)
# plt.title(r"""Distribution of $\gamma$ with %d samples""" % iterations)
#
# plt.hist(gamma_samples, histtype='stepfilled',
#          color = 'darkred', bins=30, alpha=0.8, density=True)
# plt.ylabel('Probability Density')
#
#
# plt.subplot(212)
# plt.title(r"""Distribution of $a$ with %d samples""" % iterations)
# plt.hist(a_samples, histtype='stepfilled',
#          color = 'darkblue', bins=30, alpha=0.8, density=True)
# plt.ylabel('Probability Density')
# plt.show()
#
# figdis.savefig('c:/Liang/Googlebox/Research/Project2/python_p2/posterior.png', dpi=figdis.dpi)
#
# # Trace plots
# figtra = plt.figure(figsize=(12, 6))
#
#
# # Plot alpha trace
# plt.subplot(211)
# plt.title(r'Trace of $\gamma$')
# plt.plot(gamma_samples, color = 'darkred')
# plt.xlabel('Samples'); plt.ylabel('Parameter')
#
# # Plot beta trace
# plt.subplot(212)
# plt.title(r'Trace of $a$')
# plt.plot(a_samples, color='b')
# plt.xlabel('Samples'); plt.ylabel('Parameter')
# plt.tight_layout(h_pad=0.8)
# plt.show()
# figtra.savefig('c:/Liang/Googlebox/Research/Project2/python_p2/posterior_trace.png', dpi=figtra.dpi)

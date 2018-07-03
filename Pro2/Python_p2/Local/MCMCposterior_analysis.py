import numpy as np
import matplotlib.pyplot as plt
import pylab as P
import matplotlib.mlab as mlab
from sklearn.neighbors import KernelDensity
# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

cal_size = 40000
dir = 'c:/Liang/Googlebox/Research/Project2/DVmodel/1stClusterStudy/'
calfile = dir+'cal4w_DVseed13.txt'
collection = np.loadtxt(calfile)

# Raw data analysis
# distance distribution contrasting sorted and unsorted
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
data = np.array(collection[:,1])
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

# Filtered data analysis
# distance distribution contrasting sorted and unsorted
P.figure()
dis_data_filtered = filtered_coll[:,[2,3]]
n1, bins1, patches1 = P.hist(dis_data_filtered, 15, density=1, histtype='bar',
                            color=['crimson', 'burlywood'],
                            label=['distance', 'sorted distance'])
P.legend()
plt.show()

#Estimate prior distribution of parameters
gamma_filtered = np.array(filtered_coll[:,0])
gamma_filtered = gamma_filtered.reshape(-1,1)
a_filtered = np.array(filtered_coll[:,1])
a_filtered = a_filtered.reshape(-1,1)
# Plot the true distribution
x = np.linspace(0, 1, 100)[:, np.newaxis]
norm_vals = mlab.normpdf(x, 0.1, 0.4)
plt.subplot(121)
plt.title('Distribution of $\gamma$')
plt.plot(x, norm_vals)
# Plot the data using a normalized histogram
plt.hist(gamma_filtered, 50, density=True)
# Do kernel density estimation
kd = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(gamma_filtered)
# Plot the estimated densty
kd_vals = np.exp(kd.score_samples(x))
plt.plot(x, kd_vals,'--r')

plt.subplot(122)
plt.title('Distribution of a')
plt.plot(x, norm_vals)
# Plot the data using a normalized histogram
plt.hist(a_filtered, 50, density=True)
# Do kernel density estimation
kd = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(a_filtered)
# Plot the estimated densty
kd_vals = np.exp(kd.score_samples(x))
plt.plot(x, kd_vals,'--r')

# Show the plots
plt.show()



posfile = dir + 'posterior4w_DV.txt'
# ABC_MCMC analysis
iterations = 40000
# Statistic
posterior = np.loadtxt(posfile)

# Distribution plots for parameters
gamma_samples = posterior[5000::,0]
a_samples = posterior[5000::, 1]
figdis = plt.figure(figsize=(12, 8))

plt.subplot(211)
plt.title(r"""Distribution of $\gamma$ with %d samples""" % iterations)

plt.hist(gamma_samples, histtype='stepfilled',
         color = 'darkred', bins=30, alpha=0.8, density=True)
plt.ylabel('Probability Density')


plt.subplot(212)
plt.title(r"""Distribution of $a$ with %d samples""" % iterations)
plt.hist(a_samples, histtype='stepfilled',
         color = 'darkblue', bins=30, alpha=0.8, density=True)
plt.ylabel('Probability Density')
plt.show()

figdir1 = dir + 'posterior.png'
figdis.savefig(figdir1, dpi=figdis.dpi)


# Trace plots
figtra = plt.figure(figsize=(12, 6))
# Plot alpha trace
plt.subplot(211)
plt.title(r'Trace of $\gamma$')
plt.plot(gamma_samples, color = 'darkred')
plt.xlabel('Samples'); plt.ylabel('Parameter')

# Plot beta trace
plt.subplot(212)
plt.title(r'Trace of $a$')
plt.plot(a_samples, color='b')
plt.xlabel('Samples'); plt.ylabel('Parameter')
plt.tight_layout(h_pad=0.8)
plt.show()


figdir2 = dir + 'posterior_trace.png'
figtra.savefig(figdir2, dpi=figtra.dpi)

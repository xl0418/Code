import numpy as np
import matplotlib.pyplot as plt
import pylab as P
import matplotlib.mlab as mlab
from sklearn.neighbors import KernelDensity
# Observation parameters [gamma,a]
par_obs = np.array([0.001,0.1])

cal_size = 20000
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
datagamma = np.array(collection[:,0])
datagamma = datagamma.reshape(-1,1)
dataa = np.array(collection[:,1])
dataa = dataa.reshape(-1,1)
# Plot the true distribution
x = np.linspace(0, 1, 100)[:, np.newaxis]
# norm_vals = mlab.normpdf(x, 0.1, 0.4)
# plt.plot(x, norm_vals)
# Plot the data using a normalized histogram
# Do kernel density estimation
# kd = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(data)
# Plot the estimated densty
# kd_vals = np.exp(kd.score_samples(x))
# plt.plot(x, kd_vals,'--r')
fig_raw = plt.figure(figsize=(12, 6))
# Plot alpha trace
plt.subplot(121)
plt.title('Prior $\gamma$ in raw calibration')
plt.hist(datagamma, 50, density=True)
plt.axvline(par_obs[0],color = 'r',linestyle = '--')
plt.xlabel('Samples')
plt.ylabel('Frequency')

# Plot beta trace
plt.subplot(122)
plt.title('Prior a in raw calibration')
plt.hist(dataa, 50, density=True)
plt.axvline(par_obs[1],color = 'r',linestyle = '--')
plt.xlabel('Samples')
plt.ylabel('')

plt.show()

figraw = dir + 'calraw.png'
fig_raw.savefig(figraw, dpi=fig_raw.dpi)



# Filtered data analysis
sort = 0
if sort == 0:
    ind = 2
else:
    ind = 3
threshold = 0.05
num = threshold*cal_size-1
delta = np.sort(collection[:,ind])[int(num)]
mn,idx = min( (collection[i,ind],i) for i in range(len(collection[:,ind])) )
startvalue_par = collection[idx,:2]

filtered_coll = collection[collection[:,ind]<=delta]


#Estimate prior distribution of parameters
gamma_filtered = np.array(filtered_coll[:,0])
gamma_filtered = gamma_filtered.reshape(-1,1)
a_filtered = np.array(filtered_coll[:,1])
a_filtered = a_filtered.reshape(-1,1)
# Plot the true distribution
fig_filter = plt.figure(figsize=(12, 6))
x = np.linspace(0, 1, 100)[:, np.newaxis]
plt.subplot(121)
plt.title('Prior $\gamma$ in filtered calibration')
plt.hist(gamma_filtered, 50, density=True)
plt.axvline(par_obs[0],color = 'r',linestyle = '--')


plt.subplot(122)
plt.title('Prior a in filtered calibration')
plt.hist(a_filtered, 50, density=True)
plt.axvline(par_obs[1],color = 'r',linestyle = '--')
plt.show()

figfil = dir + 'calfilter.png'
fig_filter.savefig(figfil, dpi=fig_filter.dpi)



# ABC_MCMC analysis

posfile = dir + 'posterior4w_DV.txt'
iterations = 5000
# Statistic
posterior = np.loadtxt(posfile)

# Distribution plots for parameters
gamma_samples = posterior[::10,0]
a_samples = posterior[::10, 1]
figdis = plt.figure(figsize=(12, 8))

plt.subplot(211)
plt.title(r"""Distribution of $\gamma$ with %d samples""" % iterations)

plt.hist(gamma_samples, histtype='stepfilled',
         color = 'darkred', bins=30, alpha=0.8, density=True)
plt.axvline(par_obs[0],color = 'r',linestyle = '--')
plt.ylabel('Probability Density')


plt.subplot(212)
plt.title(r"""Distribution of $a$ with %d samples""" % iterations)
plt.hist(a_samples, histtype='stepfilled',
         color = 'darkblue', bins=30, alpha=0.8, density=True)
plt.ylabel('Probability Density')
plt.axvline(par_obs[1],color = 'r',linestyle = '--')

plt.show()

figdir1 = dir + '4wposterior_distribution.png'
figdis.savefig(figdir1, dpi=figdis.dpi)


# Trace plots
figtra = plt.figure(figsize=(12, 6))
# Plot alpha trace
plt.subplot(211)
plt.title(r'Trace of $\gamma$')
plt.plot(gamma_samples, color = 'darkred')
plt.axhline(par_obs[0],color = 'g',linestyle = '--')
plt.xlabel('Samples'); plt.ylabel('Parameter')

# Plot beta trace
plt.subplot(212)
plt.title(r'Trace of $a$')
plt.plot(a_samples, color='b')
plt.xlabel('Samples'); plt.ylabel('Parameter')
plt.tight_layout(h_pad=0.8)
plt.axhline(par_obs[1],color = 'g',linestyle = '--')

plt.show()


figdir2 = dir + '4wposterior_trace.png'
figtra.savefig(figdir2, dpi=figtra.dpi)

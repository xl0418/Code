import numpy as np
import matplotlib.pyplot as plt

def sample_ztp(lam):

    k = 1
    t = np.exp(-lam) / (1 - np.exp(-lam)) * lam
    s = t
    u = np.random.uniform(0,1)
    while s < u:
        k += 1
        t *= lam / k
        s += t
    return k

sample_ztp(10e9)

sample = [sample_ztp(1) for i in range(100000)]
print(np.mean(sample),np.var(sample))
plt.hist(sample)
sample_poisson = [np.random.poisson(1) for i in range(100000)]
sample_poisson_array = np.asarray(sample_poisson)
sample_poisson_array=sample_poisson_array[np.where(sample_poisson_array>0)]
print(np.mean(sample_poisson_array),np.var(sample_poisson_array))
plt.hist(sample_poisson_array)
print(np.mean(sample_poisson),np.var(sample_poisson))


sample_ztp(np.array([1,2]))

sam_ztp_vec = np.vectorize(sample_ztp)
lamb = [1,2,3,4,5]
sam_ztp_vec(lamb)
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from testABC_sim import DVtraitsim_tree
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'

outestresult = DVtraitsim_tree(file = file, replicate = 13,theta = 0, gamma1 = 0.01, r = 1, a = 0.01,scalar = 1000, K = 100000, nu = 0.0001, Vmax = 1 )

trait_RI_dr = outestresult[0]
population_RI_dr = outestresult[1]
evo_time, total_species = trait_RI_dr.shape
evo_time = evo_time - 1
trait_dr_tips = trait_RI_dr[evo_time, :][~np.isnan(trait_RI_dr[evo_time, :])]
population_tips = population_RI_dr[evo_time, :][~np.isnan(population_RI_dr[evo_time, :])]
trait_RI_dr[np.where(trait_RI_dr == 0)[0], np.where(trait_RI_dr == 0)[1]] = None
population_RI_dr[np.where(population_RI_dr == 0)[0], np.where(population_RI_dr == 0)[1]] = None
num_plots = total_species

# pltstr = directory + "/%dgamma%da%dreplicate.png" % (count1, count2, loop)

x = np.arange(evo_time + 1)
labels = []
fig, ax = plt.subplots(nrows=1, ncols=2)
for i in range(1, num_plots + 1):
    ax[0].plot(x, trait_RI_dr[:, i - 1])
# plt.subplot(122)
for i in range(1, num_plots + 1):
    ax[1].plot(x, population_RI_dr[:,i-1])

# fig.savefig(pltstr)
# plt.close(fig)
# plt.subplot(2, 2, 2)
# for i in range(1, num_plots + 1):
#     plt.plot(x, population_RI_dr[:,i-1])
# plt.subplot(2, 2, 3)
# sns.distplot(trait_dr_tips, hist=False, rug=True)
#
# plt.subplot(2, 2, 4)
# sns.distplot(population_tips, hist=False, rug=True)
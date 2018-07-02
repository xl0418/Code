import matplotlib.pyplot as plt
import os
import numpy as np


theta = 0  # optimum of natural selection
r = 1  # growth rate
K = 5000  # carrying capacity
delta_pop = .001  # Variance of random walk of population
nu = 0.0001
Vmax = 1
scalor = 1000
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec
count1 = 0
count2 = 0
# trait evolution plot
file = '/home/p274981/Python_p2/DVmodel/'
replicate = 1000
for gamma1 in gamma_vec:
    count1 += 1
    for a in a_vec:
        count2 += 1
        directory = '/data/p274981/Project2/%dgamma%da' % (count1,count2)
        if not os.path.exists(directory):
            os.makedirs(directory)
        for loop in range(0,replicate+1):
            traitfilestr = directory + "/%dgamma%da%dreplicatetrait.txt" % (count1,count2,loop)
            popfilestr = directory + "/%dgamma%da%dreplicatetraitpop.txt" % (count1,count2,loop)
            trait_RI_dr=np.loadtxt(traitfilestr)
            population_RI_dr = np.loadtxt(popfilestr)
            evo_time, total_species = trait_RI_dr.shape
            evo_time = evo_time - 1
            trait_dr_tips = trait_RI_dr[evo_time, :][~np.isnan(trait_RI_dr[evo_time, :])]
            population_tips = population_RI_dr[evo_time, :][~np.isnan(population_RI_dr[evo_time, :])]
            trait_RI_dr[np.where(trait_RI_dr == 0)[0], np.where(trait_RI_dr == 0)[1]] = None
            population_RI_dr[np.where(population_RI_dr == 0)[0], np.where(population_RI_dr == 0)[1]] = None
            num_plots = total_species

            pltstr = directory + "/%dgamma%da%dreplicate.png" % (count1, count2, loop)

            x = np.arange(evo_time + 1)
            labels = []
            fig, ax = plt.subplots(nrows=1, ncols=1)
            for i in range(1, num_plots + 1):
                ax.plot(x, trait_RI_dr[:, i - 1])
            fig.savefig(pltstr)
            plt.close(fig)
            # plt.subplot(2, 2, 2)
            # for i in range(1, num_plots + 1):
            #     plt.plot(x, population_RI_dr[:,i-1])
            # plt.subplot(2, 2, 3)
            # sns.distplot(trait_dr_tips, hist=False, rug=True)
            #
            # plt.subplot(2, 2, 4)
            # sns.distplot(population_tips, hist=False, rug=True)
import sys, os
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from DV_model_sim_along_phy import DVtraitsim_tree
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pylab import *
from matplotlib import animation

theta = 0  # optimum of natural selection
r = 1  # growth rate
K = 5000  # carrying capacity
delta_pop = .001  # Variance of random walk of population
nu = 0.0001
Vmax = 1
scalor = 1000
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec
count = 0
# trait evolution plot
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'

for gamma1 in gamma_vec:
    for a in a_vec:
        count += 1
        print(count)
        simresult = DVtraitsim_tree(file = file, gamma1 = gamma1, a = a,scalor = 5000)
        evo_time, total_species = simresult[0].shape
        evo_time = evo_time-1
        trait_RI_dr = simresult[0]
        population_RI_dr = simresult[1]


        trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
        population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]

        trait_RI_dr[np.where(trait_RI_dr == 0)[0],np.where(trait_RI_dr == 0)[1]] = None

        population_RI_dr[np.where(population_RI_dr == 0)[0],np.where(population_RI_dr == 0)[1]] = None
        num_plots = total_species


        x = np.arange(evo_time+1)
        labels = []
        plt.subplot(6, 6, count)
        for i in range(1, num_plots + 1):
            plt.plot(x, trait_RI_dr[:,i-1])

        # plt.subplot(2, 2, 2)
        # for i in range(1, num_plots + 1):
        #     plt.plot(x, population_RI_dr[:,i-1])
        # plt.subplot(2, 2, 3)
        # sns.distplot(trait_dr_tips, hist=False, rug=True)
        #
        # plt.subplot(2, 2, 4)
        # sns.distplot(population_tips, hist=False, rug=True)
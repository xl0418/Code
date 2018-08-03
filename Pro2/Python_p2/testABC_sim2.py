from matplotlib.pylab import *
import numpy as np

def DVtraitsim_tree(file, model, replicate = 0,theta = 0, gamma1 = 0.001, r = 1, a = 0.01,scalar = 1000, K = 100000, nu = 0.0001, Vmax = 1 ):
    valid = True
    if replicate > 0:
        np.random.seed(replicate)  # set random seed
    file = file
    # load tree data
    timelist = np.genfromtxt(file+'timelist.csv', delimiter=',')
    timebranch = np.genfromtxt(file+'timebranch.csv', delimiter=',')
    timeend = np.genfromtxt(file+'timeend.csv', delimiter=',')
    traittable = np.genfromtxt(file+'traittable.csv', delimiter=',')
    ltable = np.genfromtxt(file+'Ltable.csv', delimiter=',')
    # processing data
    ltable = np.delete(ltable, (0), axis=0)
    ltable = np.delete(ltable, (0), axis=1)
    traittable = np.delete(traittable, (0), axis=0)
    traittable = np.delete(traittable, (0), axis=1)
    daughter_index = np.absolute(ltable[:,2])
    daughter_index = [int(x) for x in daughter_index]

    parent_index = np.absolute(ltable[:,1])
    parent_index = [int(x) for x in parent_index]

    evo_timelist = timelist[:,1]
    evo_timelist = np.delete(evo_timelist,0)
    evo_timelist = max(evo_timelist)-evo_timelist
    evo_timelist = evo_timelist * scalar
    evo_timelist = evo_timelist.astype(int)

    timebranch = timebranch[:,1]
    timebranch = np.delete(timebranch,0)
    timebranch = [int(x)-1 for x in timebranch]
    timeend = timeend[:,1]
    timeend = np.delete(timeend,0)
    timeend = [int(x)-1 for x in timeend]

    # evolution time: speciation time
    evo_time = max(evo_timelist)

    speciate_time = evo_timelist[timebranch]
    extinct_time = evo_timelist[timeend]
    extinct_time[np.where(extinct_time == evo_time)[0]] = -1

    extinct_time = np.delete(extinct_time, np.where(extinct_time == evo_time)[0])
    total_species = len(speciate_time)

    # Initialize trait evolution and population evolution matrices
    trait_RI_dr = np.zeros((2, total_species))   # population
    trait_RI_dr[1,] = (model==0)*np.random.normal(loc = gamma1,scale = a,size = total_species) +\
    (model == 1)* np.random.uniform(low=gamma1, high=a, size=total_species)
    population_RI_dr = np.zeros((2, total_species))   # population
    V = np.zeros((2, total_species))          # trait vairance

    return trait_RI_dr, population_RI_dr, valid,V


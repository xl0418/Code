from matplotlib.pylab import *
import numpy as np

def DVtraitsim_tree(file, replicate = 0,theta = 0, gamma1 = 0.001, r = 1, a = 0.01,scalar = 1000, K = 100000, nu = 0.0001, Vmax = 1 ):
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
    trait_RI_dr = np.zeros((evo_time + 1, total_species))    # trait
    population_RI_dr = np.zeros((evo_time + 1, total_species))   # population
    V = np.zeros((evo_time + 1, total_species))          # trait vairance

    #  initialize condition for species trait and population
    trait_RI_dr[0, (0,1)] = 0    #trait for species
    mu_pop, sigma_pop = 500, 10  # mean and standard deviation
    population_RI_dr[0, (0,1)] = np.random.normal(mu_pop, sigma_pop, 2)
    V[0] = (1 / total_species)
    # Existing species matrix
    existing_species = traittable

    for i in range(evo_time):
        num_event = len(np.where(evo_timelist <= i)[0])
        node = num_event - 2
        index_existing_species = np.where(existing_species[node] == 1)[0]

        # trait-population coevolution model
        K_RI_dr = K   # constant carrying capacity

        var_trait = a
        trait_RI_dr[i + 1, index_existing_species] = trait_RI_dr[i, index_existing_species] + gamma1 * (theta -
                                                                    trait_RI_dr[i, index_existing_species]) +\
                                                     np.random.normal(0, var_trait, len(index_existing_species))
        possion_lambda = population_RI_dr[i, index_existing_species] * r * \
                                                          np.exp(1 - np.sum(np.array(population_RI_dr[i, index_existing_species])) / K_RI_dr)

        population_RI_dr[i + 1, index_existing_species] = np.random.poisson(lam = possion_lambda[0],size = (1,len(index_existing_species)))

        # population_RI_dr[i + 1, np.where(population_RI_dr[i + 1] < 1)] = 0
        ext_index_RI_dr = np.where(population_RI_dr[i + 1,index_existing_species] == 0)[0]
        negative_v = np.where(V[i + 1, index_existing_species] < 0)[0]
        if len(ext_index_RI_dr) > 0:
            valid = False
            print('Inconsistent zero population')
            break
        if len(negative_v) > 0:
            valid = False
            print('Negative variance')
            break
        if (i+1) in speciate_time:
            spe_event_index = len(np.where(speciate_time <= (i+1))[0])
            trait_RI_dr[i+1,daughter_index[spe_event_index-1]-1] = trait_RI_dr[i+1,parent_index[spe_event_index-1]-1]
            population_RI_dr[i + 1, daughter_index[spe_event_index-1]-1] =1/2 * population_RI_dr[i + 1,
                                                                                           parent_index[spe_event_index-1]-1]
            population_RI_dr[i + 1, parent_index[spe_event_index-1]-1] = 1 / 2 * population_RI_dr[i + 1,
                                                                                            parent_index[spe_event_index-1]-1]

        if (i + 1) in extinct_time:
            extinct_species = int(np.where(extinct_time == (i+1))[0])
            trait_RI_dr[i+1, extinct_species] = None
            population_RI_dr[i+1, extinct_species] = 0
    row_ext = np.where(population_RI_dr == 0)[0]
    col_ext = np.where(population_RI_dr == 0)[1]
    trait_RI_dr[row_ext,col_ext] = None
    population_RI_dr[row_ext,col_ext] = None
    V[row_ext,col_ext] = None
    return trait_RI_dr, population_RI_dr, valid,V


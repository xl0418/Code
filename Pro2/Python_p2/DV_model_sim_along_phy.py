from matplotlib.pylab import *
import numpy as np


# Natural selection function
def ga(gamma, theta, zi):
    return -gamma * (theta - zi) ** 2


# Dynamic carrying capacity
def Kd(gamma_K, theta, zi, K):
    return max(K * np.exp(-gamma_K * (theta - zi) ** 2), 1)


# Competition function
def beta(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret


# Derivative of the competition function
def sigma(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(
            2 * a * (zi[n1] - np.array(zj)) * np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret


# Derivative of the competition function
def sigmasqr(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(
            4 * a ** 2 * (zi[n1] - np.array(zj)) ** 2 * np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret


# Delta function
def Delta(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(
            -2 * a + 4 * a ** 2 * (zi[n1] - np.array(zj)) ** 2 * np.exp(-a * (zi[n1] - np.array(zj)) ** 2)
            * np.array(nj))
    return zi_ret


def DVtraitsim_tree(file, replicate = 0,theta = 0, gamma1 = 0.001, r = 1, a = 0.01,scalor = 1000, K = 5000, delta_pop = 0.001, nu = 0.0001, Vmax = 1 ):
    if replicate == 1:
        np.random.seed(13)  # set random seed
    file = file
    # load data
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
    evo_timelist = evo_timelist * scalor
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
    trait_RI_dr = np.zeros((evo_time + 1, total_species))
    population_RI_dr = np.zeros((evo_time + 1, total_species))

    #  initialize condition for species trait and population
    trait_RI_dr[0, (0,1)] = 0
    # print(trait_RI_dr[0])
    mu_pop, sigma_pop = 500, 10  # mean and standard deviation
    population_RI_dr[0, (0,1)] = np.random.normal(mu_pop, sigma_pop, 2)
    # print(population_RI_dr[0])
    V = np.zeros((evo_time + 1, total_species))
    V.fill(1 / total_species)
    # Existing species matrix
    existing_species = traittable

    for i in range(evo_time):
        # print(i)
        num_event = len(np.where(evo_timelist <= i)[0])
        node = num_event - 2
        index_existing_species = np.where(existing_species[node] == 1)[0]

        # RI dynamic r model
        K_RI_dr = K
        Beta_RI_dr = beta(a=a, zi=trait_RI_dr[i, index_existing_species], zj=trait_RI_dr[i, index_existing_species],
                          nj=population_RI_dr[i, index_existing_species])
        Sigma_RI_dr = sigma(a=a, zi=trait_RI_dr[i, index_existing_species], zj=trait_RI_dr[i, index_existing_species],
                            nj=population_RI_dr[i, index_existing_species])
        Sigmasqr_RI_dr = sigmasqr(a=a, zi=trait_RI_dr[i, index_existing_species], zj=trait_RI_dr[i, index_existing_species],
                                  nj=population_RI_dr[i, index_existing_species])
        Delta_RI_dr = Delta(a=a, zi=trait_RI_dr[i, index_existing_species], zj=trait_RI_dr[i, index_existing_species],
                            nj=population_RI_dr[i, index_existing_species])
        var_trait = V[i, index_existing_species] / (2 * population_RI_dr[i, index_existing_species])
        trait_RI_dr[i + 1, index_existing_species] = trait_RI_dr[i, index_existing_species] + V[i, index_existing_species] * \
                                                     (2 * gamma1 * (theta -
                                                                    trait_RI_dr[i, index_existing_species]) + 1 / K_RI_dr *
                                                      Sigma_RI_dr) \
                                                     + np.random.normal(0, var_trait, len(index_existing_species))
        population_RI_dr[i + 1, index_existing_species] = population_RI_dr[i, index_existing_species] * r * \
                                                          np.exp(-gamma1 * (
                                                                      theta - trait_RI_dr[i, index_existing_species]) ** 2 +
                                                                 (1 - Beta_RI_dr / K_RI_dr) \
                                                                 + np.random.normal(0, delta_pop,
                                                                                    len(index_existing_species)))
        V[i + 1, index_existing_species] = V[i, index_existing_species] + V[i, index_existing_species] ** 2 * (
                    -2 * gamma1 + 4 * gamma1 ** 2 * (theta - trait_RI_dr[i, index_existing_species]) ** 2 +
                    1 / K_RI_dr *  # V^2 w''/w
                    (Sigma_RI_dr - Sigmasqr_RI_dr) + 4 * gamma1 / K_RI_dr * (
                            theta - trait_RI_dr[i, index_existing_species]) * Sigma_RI_dr +
                    Sigma_RI_dr ** 2 / K_RI_dr ** 2) - \
                                           V[i, index_existing_species] / 2 + 2 * population_RI_dr[
                                               i, index_existing_species] * nu * Vmax / (
                                                   1 + 4 * population_RI_dr[
                                               i, index_existing_species] * nu)  # loss due to sexual selection; gain from mutation

        # population_RI_dr[i + 1, np.where(population_RI_dr[i + 1] < 1)] = 0
        ext_index_RI_dr = np.where(population_RI_dr[i + 1,index_existing_species] == 0)[0]
        if len(ext_index_RI_dr) > 0:
            trait_RI_dr = False
            population_RI_dr = False
            break

        if (i+1) in speciate_time:
            spe_event_index = len(np.where(speciate_time <= (i+1))[0])
            trait_RI_dr[i+1,daughter_index[spe_event_index-1]-1] = trait_RI_dr[i+1,parent_index[spe_event_index-1]-1]
            population_RI_dr[i + 1, daughter_index[spe_event_index-1]-1] =1/2 * population_RI_dr[i + 1,
                                                                                           parent_index[spe_event_index-1]-1]
            population_RI_dr[i + 1, parent_index[spe_event_index-1]-1] = 1 / 2 * population_RI_dr[i + 1,
                                                                                            parent_index[spe_event_index-1]-1]
            V[i+1,daughter_index[spe_event_index-1]-1] = 1/2 * V[i+1,parent_index[spe_event_index-1]-1]
            V[i+1,parent_index[spe_event_index-1]-1] = 1/2 * V[i+1,parent_index[spe_event_index-1]-1]

        if (i + 1) in extinct_time:
            extinct_species = int(np.where(extinct_time == (i+1))[0])
            trait_RI_dr[i+1, extinct_species] = None
            population_RI_dr[i+1, extinct_species] = 0
    trait_RI_dr[np.where(trait_RI_dr == 0)[0],np.where(trait_RI_dr == 0)[1]] = None
    return trait_RI_dr, population_RI_dr


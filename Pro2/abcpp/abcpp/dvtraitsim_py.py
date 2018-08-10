import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam


def competition_functions(a, zi, nj):
    """ competition functions.

    returns beta = Sum_j( exp(-a(zi-zj)^2) * Nj)
            sigma = Sum_j( 2a * (zi-zj) * exp(-a(zi-zj)^2) * Nj)
            sigmaSqr = Sum_j( 4a^2 * (zi-zj)^2 * exp(-a(zi-zj)^2) * Nj)
    """
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t1 = np.exp(-a * T ** 2) * nj
    t2 = (2 * a) * T
    beta = np.sum(t1, axis=1)
    sigma = np.sum(t2 * t1, axis=1)
    sigmasqr = np.sum(t2 ** 2 * t1, axis=1)
    return beta, sigma, sigmasqr


# Sample function within a specific range (0,1)
def PopsplitNormal(mean, sigma):
    while True:
        x = np.random.normal(mean,sigma,1)
        if (x>0 and x<1):
           return x


def DVSim(td, param):
    # parameters from DVParam
    gamma = param[0]
    a = param[1]
    K = param[2]
    nu = param[3]
    r = param[4]
    theta = param[5]
    Vmax = param[6]
    inittrait = param[7]
    initpop = param[8]
    valid = True

    # Initialize trait evolution and population evolution matrices
    trait_RI_dr = np.zeros((td.evo_time + 1, td.total_species))  # trait
    population_RI_dr = np.zeros((td.evo_time + 1, td.total_species))  # population
    V = np.zeros((td.evo_time + 1, td.total_species))  # trait vairance

    #  initialize condition for species trait and population
    trait_RI_dr[0, (0, 1)] = inittrait  # trait for species
    mu_pop, sigma_pop = initpop, 10  # mean and standard deviation
    population_RI_dr[0, (0, 1)] = np.random.normal(mu_pop, sigma_pop, 2)
    V[0] = (1 / td.total_species)
    existing_species = td.traittable
    node = 0;
    next_event = td.events[node];
    idx = np.where(existing_species[node] == 1)[0]    # existing species
    # trait-population coevolution model
    for i in range(td.evo_time):
        # pull current state
        Ni = population_RI_dr[i, idx]
        Vi = V[i, idx]
        zi = trait_RI_dr[i, idx]
        Ki = K
        dtz = theta - zi
        beta, sigma, sigmasqr = competition_functions(a, zi, Ni)

        # update
        var_trait = Vi / (2 * Ni)
        trait_RI_dr[i + 1, idx] = zi + Vi * (2 * gamma * dtz + 1 / Ki * sigma) + np.random.normal(0, var_trait, len(idx))
        possion_lambda = Ni * r * np.exp(-gamma * dtz**2 + (1 - beta / Ki))
        population_RI_dr[i + 1, idx] = np.random.poisson(lam=possion_lambda)
        V[i + 1, idx] = Vi / 2 + 2 * Ni * nu * Vmax / (1 + 4 * Ni * nu) \
                        + Vi ** 2 * (
                            -2 * gamma + 4 * gamma**2 * dtz ** 2 +
                                1 / Ki * (sigma - sigmasqr) + 4 * gamma / Ki *
                                dtz * sigma + sigma ** 2 / Ki**2
                            )
        # sanity check
        if np.any(V[i + 1, idx] < 0):
            valid = False
            print('Inconsistent negative var', i)
            break
        if np.any(population_RI_dr[i + 1, idx] < 1):
            valid = False
            print('Inconsistent zero population', i)
            break
        # events
        if (i + 1) == next_event[0]:
            parent = next_event[1]
            daughter = next_event[2]
            if (daughter == -1):
                # extinction
                extinct_species = next_event[1]
                trait_RI_dr[i + 1, extinct_species] = None
                population_RI_dr[i + 1, extinct_species] = 0
            else:
                # speciation
                trait_RI_dr[i + 1, daughter] = trait_RI_dr[i + 1, parent]
                splitratio = PopsplitNormal(mean=0.5, sigma=0.2)
                tmp = population_RI_dr[i + 1, parent]
                population_RI_dr[i + 1, parent] = splitratio * tmp
                population_RI_dr[i + 1, daughter] = (1 - splitratio) * tmp
                V[i + 1, parent] = 1 / 2 * V[i + 1, parent]
                V[i + 1, daughter] = V[i + 1, parent]
            # advance to next event/node
            node = node + 1
            next_event = td.events[node];
            idx = np.where(existing_species[node] == 1)[0]

    row_ext = np.where(population_RI_dr == 0)[0]
    col_ext = np.where(population_RI_dr == 0)[1]
    trait_RI_dr[row_ext, col_ext] = None
    population_RI_dr[row_ext, col_ext] = None
    V[row_ext, col_ext] = None
    return trait_RI_dr, population_RI_dr, valid, V



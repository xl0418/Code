import numpy as np

def competition_functions_candi(a, zi):
    """ competition functions.

    returns beta = Sum_j( exp(-a(zi-zj)^2) * Nj)
            sigma = Sum_j( 2a * (zi-zj) * exp(-a(zi-zj)^2) * Nj)
            sigmaSqr = Sum_j( 4a^2 * (zi-zj)^2 * exp(-a(zi-zj)^2) * Nj)
    """
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t1 = np.exp(-a * T ** 2)
    t2 = np.sign(T)
    beta = np.sum(t1, axis=1)
    sigma = np.sum(t2 * t1, axis=1)
    return beta, sigma


def Candimodels(td, param):
    # parameters from DVParam
    gamma = param[0]
    a = param[1]
    theta = param[2]
    m = param[3]
    inittrait = param[4]
    var_trait = param[5]

    sim_evo_time = td.sim_evo_time
    events = td.sim_events

    # Initialize trait evolution evolution matrices
    trait_RI_dr = np.zeros((sim_evo_time + 1, td.total_species))  # trait

    #  initialize condition for species trait
    trait_RI_dr[0, (0, 1)] = inittrait  # trait for species
    existing_species = td.traittable
    node = 0
    next_event = events[node]
    idx = np.where(existing_species[node] == 1)[0]    # existing species

    # trait-population coevolution model
    for i in range(sim_evo_time):
        # pull current state
        zi = trait_RI_dr[i, idx]
        dtz = theta - zi
        beta, sigma = competition_functions_candi(a, zi)

        # update
        trait_RI_dr[i + 1, idx] = zi +  gamma * dtz + m * sigma + np.random.normal(0.0, var_trait,size = len(zi))

        # events
        while (i + 1) == next_event[0]:
            daughter = next_event[2]
            if (daughter == -1):
                # extinction
                extinct_species = next_event[1]
                trait_RI_dr[i + 1, extinct_species] = None
            else:
                # speciation
                parent = next_event[1]
                trait_RI_dr[i + 1, daughter] = trait_RI_dr[i + 1, parent]
            # advance to next event/node
            node = node + 1
            next_event = events[node]
            idx = np.where(existing_species[node] == 1)[0]

    return { 'sim_time': i + 1,  'Z': trait_RI_dr }

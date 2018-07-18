from matplotlib.pylab import *
import numpy as np


class DVTreeData:
    def __init__(self, path, scalar):
        self.path = path
        # load tree data
        self.timelist = self._from_txt(path + 'timelist.csv')
        self.timebranch = self._from_txt(path + 'timebranch.csv')
        self.timeend = self._from_txt(path + 'timeend.csv')
        self.traittable = self._from_txt(path + 'traittable.csv')
        self.ltable = self._from_txt(path + 'Ltable.csv')
        # derived data
        self.parent_index = np.absolute(self.ltable[:, 1]).astype(int)
        self.daughter_index = np.absolute(self.ltable[:, 2]).astype(int)
        self.evo_timelist = (scalar * (max(self.timelist[:, 0]) - self.timelist[:, 0])).astype(int)
        self.timebranch = self.timebranch[:, 0].astype(int) - 1
        self.timeend = self.timeend[:, 0].astype(int) - 1

    # returns trimmed table as numpy.ndarray
    def _from_txt(self, file):
        tmp = np.genfromtxt(file, delimiter=',', skip_header=1)
        return np.delete(tmp, (0), axis=1)


# competition functions
# returns beta = Sum_j( exp(-a(zi-zj)^2) * Nj)
#         sigma = Sum_j( 2a * (zi-zj) * exp(-a(zi-zj)^2) * Nj)
#         sigmaSqr = Sum_j( 4a^2 * (zi-zj)^2 * exp(-a(zi-zj)^2) * Nj)
def competition_functions(a, zi, nj):
    T = zi[:, np.newaxis] - zi  # trait-distance matrix (via 'broadcasting')
    t1 = np.exp(-a * T ** 2) * nj
    t2 = (2 * a) * T
    beta = np.sum(t1, axis=1)
    sigma = np.sum(t2 * t1, axis=1)
    sigmasqr = np.sum(t2 ** 2 * t1, axis=1)
    return beta, sigma, sigmasqr


def DVtraitsim_tree(file, replicate=0, theta=0, gamma1=0.001, r=1, a=0.01, scalar=1000, K=100000, nu=0.0001, Vmax=1):
    valid = True
    if replicate > 0:
        np.random.seed(replicate)  # set random seed
    td = DVTreeData(file, scalar)

    # evolution time: speciation time
    evo_time = max(td.evo_timelist)
    speciate_time = td.evo_timelist[td.timebranch]
    extinct_time = td.evo_timelist[td.timeend]
    extinct_time[np.where(extinct_time == evo_time)[0]] = -1
    extinct_time = np.delete(extinct_time, np.where(extinct_time == evo_time)[0])
    total_species = len(speciate_time)

    # Initialize trait evolution and population evolution matrices
    trait_RI_dr = np.zeros((evo_time + 1, total_species))  # trait
    population_RI_dr = np.zeros((evo_time + 1, total_species))  # population
    V = np.zeros((evo_time + 1, total_species))  # trait vairance

    #  initialize condition for species trait and population
    trait_RI_dr[0, (0, 1)] = 0  # trait for species
    mu_pop, sigma_pop = 500, 10  # mean and standard deviation
    population_RI_dr[0, (0, 1)] = np.random.normal(mu_pop, sigma_pop, 2)
    V[0] = (1 / total_species)
    # Existing species matrix
    existing_species = td.traittable

    for i in range(evo_time):
        num_event = len(np.where(td.evo_timelist <= i)[0])
        node = num_event - 2

        # trait-population coevolution model
        # pull current state
        idx = np.where(existing_species[node] == 1)[0]  # index existing species
        Ni = population_RI_dr[i, idx]
        Vi = V[i, idx]
        zi = trait_RI_dr[i, idx]
        Ki = K
        dtz = theta - zi
        beta, sigma, sigmasqr = competition_functions(a=a, zi=zi, nj=Ni)

        # update
        var_trait = Vi / (2 * Ni)
        trait_RI_dr[i + 1, idx] = zi + Vi * (2 * gamma1 * (theta - zi) + 1 / Ki * sigma) + np.random.normal(0,
                                                                                                            var_trait,
                                                                                                            len(idx))
        possion_lambda = Ni * r * np.exp(-gamma1 * (theta - zi) ** 2 + (1 - beta / Ki))
        population_RI_dr[i + 1, idx] = np.random.poisson(lam=possion_lambda, size=(1, len(idx)))
        V[i + 1, idx] = Vi / 2 + 2 * Ni * nu * Vmax / (1 + 4 * Ni * nu) \
                        + Vi ** 2 * (
                            -2 * gamma1 + 4 * gamma1 ** 2 * dtz ** 2 +
                             1 / Ki * (sigma - sigmasqr) + 4 * gamma1 / Ki *
                             dtz * sigma + sigma ** 2 / Ki ** 2
                            )

        # population_RI_dr[i + 1, np.where(population_RI_dr[i + 1] < 1)] = 0
        ext_index_RI_dr = np.where(population_RI_dr[i + 1, idx] == 0)[0]
        negative_v = np.where(V[i + 1, idx] < 0)[0]
        if len(ext_index_RI_dr) > 0:
            valid = False
            print('Inconsistent zero population')
            break
        if len(negative_v) > 0:
            valid = False
            print('Negative variance')
            break
        if (i + 1) in speciate_time:
            spe_event_index = len(np.where(speciate_time <= (i + 1))[0])
            trait_RI_dr[i + 1, td.daughter_index[spe_event_index - 1] - 1] = trait_RI_dr[
                i + 1, td.parent_index[spe_event_index - 1] - 1]
            population_RI_dr[i + 1, td.daughter_index[spe_event_index - 1] - 1] = 1 / 2 * population_RI_dr[i + 1,
                                                                                                           td.parent_index[
                                                                                                               spe_event_index - 1] - 1]
            population_RI_dr[i + 1, td.parent_index[spe_event_index - 1] - 1] = 1 / 2 * population_RI_dr[i + 1,
                                                                                                         td.parent_index[
                                                                                                             spe_event_index - 1] - 1]
            V[i + 1, td.daughter_index[spe_event_index - 1] - 1] = 1 / 2 * V[
                i + 1, td.parent_index[spe_event_index - 1] - 1]
            V[i + 1, td.parent_index[spe_event_index - 1] - 1] = 1 / 2 * V[
                i + 1, td.parent_index[spe_event_index - 1] - 1]

        if (i + 1) in extinct_time:
            extinct_species = int(np.where(extinct_time == (i + 1))[0])
            trait_RI_dr[i + 1, extinct_species] = None
            population_RI_dr[i + 1, extinct_species] = 0
    row_ext = np.where(population_RI_dr == 0)[0]
    col_ext = np.where(population_RI_dr == 0)[1]
    trait_RI_dr[row_ext, col_ext] = None
    population_RI_dr[row_ext, col_ext] = None
    V[row_ext, col_ext] = None
    return trait_RI_dr, population_RI_dr, valid, V


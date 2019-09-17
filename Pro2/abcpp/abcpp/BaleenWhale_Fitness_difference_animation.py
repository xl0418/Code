import sys, os
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
from Dvtraitsim_TVP import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from HeatMap import heatmap,annotate_heatmap
theta = 1300  # optimum of natural selection
gamma = 1.006e-06  # intensity of natural selection
r = 1  # growth rate
a = 2.346e-04  # intensity of competition
K = 10e5  # carrying capacity
kscale=1000
delta_pop = .001  # Variance of random walk of population
nu=1.319e-04
Vmax = 265
scalar = 20000
ani_gap = 200


def competition_functions_Liang(a, zi, nj):
	""" competition functions, Liang's model.

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




#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)

# parameter settings
obs_param = DVParamLiang(gamma=gamma, a=a, K=K,h=1, nu=nu, r=1, theta=theta,V00=.1,V01=.1, Vmax=Vmax, inittrait=theta, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)

for rep in range(100):
    simresult = DVSimTVP(td,obs_param)
    if simresult['sim_time'] == td.sim_evo_time:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)




if simresult['sim_time'] == td.sim_evo_time:
    evo_time, total_species = simresult['N'].shape
    evo_time = evo_time-1
    trait_RI_dr = simresult['Z']
    population_RI_dr = simresult['N']
    V = simresult['V']

    trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
    population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]

    # trait_RI_dr[np.where(trait_RI_dr == 0)[0],np.where(trait_RI_dr == 0)[1]] = np.nan
    population_RI_dr = population_RI_dr.astype(float)
    population_RI_dr[np.where(population_RI_dr == 0)[0],np.where(population_RI_dr == 0)[1]] = np.nan
    num_plots = total_species



    x = np.arange(evo_time+1)
    labels = []
    plt.subplot(2, 2, 1)
    for i in range(1, num_plots + 1):
        plt.plot(x, trait_RI_dr[:,i-1])

    plt.subplot(2, 2, 2)
    for i in range(1, num_plots + 1):
        plt.plot(x, population_RI_dr[:,i-1])
    plt.subplot(2, 2, 3)
    sns.distplot(trait_dr_tips, hist=False, rug=True)

    plt.subplot(2, 2, 4)
    sns.distplot(population_tips, hist=False, rug=True)
    plt.show()
    print(trait_RI_dr[evo_time,:])

    # Animating trait evolution along a tree
    timelist = td.timelist # np.genfromtxt(files + 'timelist.csv', delimiter=',')
    timebranch = td.timebranch # np.genfromtxt(files + 'timebranch.csv', delimiter=',')
    timeend = td.timeend # np.genfromtxt(files + 'timeend.csv', delimiter=',')
    traittable = td.traittable #np.genfromtxt(files + 'traittable.csv', delimiter=',')
    ltable = td.ltable #np.genfromtxt(files + 'Ltable.csv', delimiter=',')
    # processing data
    ltable = np.delete(ltable, (0), axis=0)
    ltable = np.delete(ltable, (0), axis=1)
    traittable = np.delete(traittable, (0), axis=0)
    traittable = np.delete(traittable, (0), axis=1)
    daughter_index = np.absolute(ltable[:, 2])
    daughter_index = [int(x) for x in daughter_index]

    parent_index = np.absolute(ltable[:, 1])
    parent_index = [int(x) for x in parent_index]

    evo_timelist = timelist[:, 0]
    evo_timelist = max(evo_timelist) - evo_timelist
    evo_timelist = evo_timelist * scalar
    evo_timelist = evo_timelist.astype(int)

    timebranch = timebranch
    # timebranch = np.delete(timebranch, 0)
    timebranch = [int(x) for x in timebranch]
    timeend = timeend
    # timeend = np.delete(timeend, 0)
    timeend = [int(x) for x in timeend]

    # evolution time: speciation time
    evo_time = max(evo_timelist)//ani_gap

    speciate_time = evo_timelist[timebranch]//ani_gap
    extinct_time = evo_timelist[timeend]//ani_gap
    extinct_time[np.where(extinct_time == evo_time)[0]] = evo_time+10

    # extinct_time = np.delete(extinct_time, np.where(extinct_time == evo_time)[0])
    total_species = len(speciate_time)
    x = np.array(range(evo_time+1))

    existing_species = td.traittable
    idx = np.where(existing_species[0] == 1)[0]

    trait_RI_dr = simresult['Z'][::ani_gap]
    population_RI_dr = simresult['N'][::ani_gap]
    V=simresult['V'][::ani_gap]
    fitness_species = np.zeros(V.shape)
    fitness_species_deltav = np.zeros(V.shape)
    fitness_difference_list = list()
    for i in range(trait_RI_dr.shape[0]):
        fitness_di = np.zeros((15,15))
        ind = np.where(np.isnan(V[i,])==False)[0]
        Ni = population_RI_dr[i,ind]
        Vi = V[i,ind]
        zi = trait_RI_dr[i,ind]
        zi_deltav = trait_RI_dr[i,ind]+np.sqrt(Vi)
        Ki = K
        dtz = theta - zi_deltav
        beta, sigma, sigmasqr = competition_functions_Liang(a, zi, Ni)
        beta_delta, sigma_delta, sigmasqr_delta = competition_functions_Liang(a, zi_deltav, Ni)

        fitness_species[i,ind] = r * np.exp(-gamma * dtz ** 2 + (1 - beta / Ki))
        fitness_species_deltav[i,ind] = r * np.exp(-gamma * dtz ** 2 + (1 - beta_delta / Ki))
        fitslice = fitness_species[i,ind]
        fitness_di[:len(ind),:len(ind)] = fitslice[:, np.newaxis] - fitslice
        np.fill_diagonal(fitness_di,fitness_species[i,] - fitness_species_deltav[i,])
        fitness_difference_list.append(abs(fitness_di))

    spe_label = ['Spe %i' % i for i in range(1,16)]

    fig, ax = plt.subplots()

    im = heatmap(fitness_difference_list[0], spe_label, spe_label, ax=ax,
                       cmap="YlGn")
    texts = annotate_heatmap(im, valfmt="",threshold=0.5)
    im.set_clim(0,0.5)
    # Create colorbar
    cbarlabel = "Difference in fitness between species"
    fig.colorbar(im,ax=ax)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im = heatmap(fitness_difference_list[0], spe_label, spe_label, ax=ax,cmap="YlGn")

    def init():
        im.set_data(np.zeros((15, 15)))


    def animate(i):
        i = (i + 1) % (len(x) + 1)
        time_text.set_text(time_template % (i))
        data = fitness_difference_list[i]
        im = heatmap(data, spe_label, spe_label, ax=ax,
                           cmap="YlGn")
        texts = annotate_heatmap(im, valfmt="")
        return im,texts

    time_template = 'Time = %d G'    # prints running simulation time
    time_text = ax.text(1.02, 1.05, '', transform=ax.transAxes)
    anim = animation.FuncAnimation(fig, animate, interval= 10,
                                   init_func=init, frames=evo_time, repeat=False,blit=False)
    # Create colorbar
    cbarlabel = "Difference in fitness between species"
    cbar = fig.colorbar(im,ax=ax)

    im.set_clim(0, 0.5)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    plt.show()

    # #
    # Writer = animation.writers['ffmpeg']
    # writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)
    # ani.save('C:\\Liang\\Googlebox\\Research\\Project2\\BaleenWhales\\BaleenWhale2w.mp4', writer=writer)
else:
    print('Junk simulation! Please try again or elevate K.')


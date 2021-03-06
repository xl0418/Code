import sys, os
sys.path.append('C:/Liang/abcpp_ms6/abcpp')
from Dvtraitsim_TVP import DVSimTVP
from Dvtraitsim_TVM import DVSimTVM
from dvtraitsim_shared import DVTreeData, DVParamLiang
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from matplotlib import animation
import pandas as pd
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from HeatMap import heatmap,annotate_heatmap
theta = 1300  # optimum of natural selection
gamma = 1.006e-06  # intensity of natural selection
r = 1  # growth rate
a = 2.346e-04  # intensity of competition
K = 10e5  # carrying capacity
kscale=1000
nu=1.1e-03
Vmax = 70
scalar = 20000
ani_gap = 10


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
    # labels = []
    # plt.subplot(2, 2, 1)
    # for i in range(1, num_plots + 1):
    #     plt.plot(x, trait_RI_dr[:,i-1])
    #
    # plt.subplot(2, 2, 2)
    # for i in range(1, num_plots + 1):
    #     plt.plot(x, population_RI_dr[:,i-1])
    # plt.subplot(2, 2, 3)
    # sns.distplot(trait_dr_tips, hist=False, rug=True)
    #
    # plt.subplot(2, 2, 4)
    # sns.distplot(population_tips, hist=False, rug=True)
    # plt.show()
    # print(trait_RI_dr[evo_time,:])

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


    ext_times_RI = []
    pop_nan20 = population_RI_dr.astype(float)
    pop_sum = np.sum(pop_nan20,axis=1)
    # evo_time, total_species = simresult[0].shape
    # evo_time = evo_time-1
    ext_spec_index_RI = np.where(population_RI_dr[evo_time,] == 0)[0]
    # ext_spec_index_RI_rs = np.where(population_RI_rs[evo_time,] == 0)[0]
    RI_traits = []
    RI_sizes = []

    for i in np.arange(total_species):
        RI_trait = trait_RI_dr[x,i]
        RI_size = population_RI_dr[x,i]
        RI_traits.append(RI_trait)
        RI_sizes.append(RI_size)


    for j in ext_spec_index_RI:
        ext_time = np.where(population_RI_dr[:,j] == 0)
        ext_times_RI.append(ext_time[0][0])

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
        fitness_di[:len(ind),:len(ind)] = abs(fitslice[:, np.newaxis] - fitslice)
        temp_mat = fitness_di[:len(ind),:len(ind)]
        temp_mat[np.where(temp_mat==0)[0],np.where(temp_mat==0)[1]]=2
        min_inter_fitness = sort(temp_mat)[:,0]
        # if(len(np.where(min_inter_fitness==0)[0])>0 and i>0):
        #     print(i)
        #     break
        np.fill_diagonal(fitness_di,abs(fitness_species[i,] - fitness_species_deltav[i,]))
        fitness_di[ind,:] = fitness_di[ind,:]/min_inter_fitness[:,None]
        fitness_difference_list.append(fitness_di)

    diff_list = pd.DataFrame([])
    for generation in range(len(fitness_difference_list)):
        diff_list_temp = pd.concat([pd.DataFrame(fitness_difference_list[generation][:, spei]) for spei in range(15)],
                                 ignore_index=True)
        diff_list = pd.concat([diff_list,diff_list_temp],ignore_index=True)

    generation_time = np.repeat([i for i in range(evo_time+1)],15*15)
    spe_label = ['Spe %i' % i for i in range(1,16)]
    spe_label_col = np.tile(np.repeat(spe_label,15),evo_time+1)
    spe_label_row = np.tile(spe_label,15*(evo_time+1))

    diff_list.insert(0,'Spe1',spe_label_col)
    diff_list.insert(1,'Spe2',spe_label_row)
    diff_list.insert(0,'Generation',generation_time)

    diff_list.rename(columns={0:'Diff'}, inplace=True)
    test_file = dir_path+'test_diff_full_ratio.csv'
    diff_list.to_csv(test_file)
    # fig, ax = plt.subplots()
    #
    # im = heatmap(fitness_difference_list[2000], spe_label, spe_label, ax=ax,
    #                    cmap="YlGn")
    # texts = annotate_heatmap(im, valfmt="",threshold=0.05)
    # im.set_clim(0,0.05)
    # # Create colorbar
    # cbarlabel = "Difference in fitness between species"
    # fig.colorbar(im,ax=ax)
    # fig.tight_layout()
    # plt.show()

    # fig, ax = plt.subplots()

    fig = figure(num=0, figsize=(10, 10))  # , dpi = 100)
    fig.suptitle("Trait Evolution", fontsize=12)
    ax01 = subplot2grid((2, 1), (0, 0))
    ax02 = subplot2grid((2, 1), (1, 0))

    # Subplot ax01 for trait evolution
    i = theta
    ax01.axhline(y=i, color='k', linestyle='--', alpha=0.7)
    ax01.set_title('Animation of the trait-population coevolution model')

    trait_max = np.nanmax(trait_RI_dr)
    trait_min = np.nanmin(trait_RI_dr)

    # set y-limits
    ax01.set_ylim(trait_min - 10, trait_max + 10)
    ax01.set_xlim(0, evo_time)

    # Turn on grids
    ax01.grid(True)
    ax01.set_xlabel("")
    ax01.set_ylabel("Trait Value")
    ax01.set_xlabel("Generation")
    popu_RI_spec_texts = []

    text_y1 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)

    for i in np.arange(total_species):
        popu_RI_spec_text1 = ax01.text(1, text_y1[i], '', transform=ax01.transAxes)
        popu_RI_spec_texts.append(popu_RI_spec_text1)

    RI_lines = []
    RI_scatters = []
    pop_line = []
    for i in np.arange(total_species):
        RI_line, = ax01.plot([], [], 'b-')
        RI_scatter = ax01.scatter([], [], s=0, c='r', alpha=0.3)
        RI_lines.append(RI_line)
        RI_scatters.append(RI_scatter)

    # Subplot ax02 for fitness difference heatmap
    im = heatmap(fitness_difference_list[0], spe_label, spe_label, cmap="YlGn")
    texts = annotate_heatmap(im, valfmt="")
    # Create colorbar
    cbarlabel = "Difference in fitness between species"
    cbar = fig.colorbar(im,ax=ax02,ticks=[0,0.5,1])

    im.set_clim(0, 1)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    cbar.ax.set_yticklabels(['Min','Medium','Max'])
    def init():
        im.set_data(np.zeros((15, 15)))


    def animate(i):
        i = (i + 1) % (len(x) + 1)
        time_text.set_text(time_template % (i))
        data = fitness_difference_list[i]
        im = heatmap(data, spe_label, spe_label, ax=ax02,
                           cmap="YlGn")
        texts = annotate_heatmap(im, valfmt="")

        RI_datas = []

        for j in np.arange(total_species):
            RI_data = np.hstack((x[i], RI_traits[j][i]))
            RI_datas.append(RI_data)
            RI_lines[j].set_data(x[0:i], RI_traits[j][0:i])
            RI_scatters[j].set_offsets(RI_datas[j])
            RI_scatters[j].set_sizes([RI_sizes[j][i] / kscale])

            # Animating labels
            # popu_RI_spec_texts[j].set_text('POS %d = %.1f' % (j+1, population_RI_dr[i,j]))
            if (i < speciate_time[j]):
                RI_lines[j].set_data([], [])

            elif (i >= speciate_time[j] and i < extinct_time[j]):
                RI_lines[j].set_data(x[speciate_time[j]:i], RI_traits[j][speciate_time[j]:i])

            else:
                RI_lines[j].set_data(x[speciate_time[j]:extinct_time[j]],
                                     RI_traits[j][speciate_time[j]:extinct_time[j]])

            if (j in ext_spec_index_RI):
                end_time_RI = ext_times_RI[np.where(j == ext_spec_index_RI)[0][0]]
                if (i >= end_time_RI):
                    RI_lines[j].set_data(x[0:end_time_RI], RI_traits[j][0:end_time_RI])
        return im,texts,RI_lines, RI_scatters

    time_template = 'Time = %d G'    # prints running simulation time
    time_text = ax02.text(1.02, 1.05, '', transform=ax02.transAxes)
    anim = animation.FuncAnimation(fig, animate, interval= 10,
                                   init_func=init, frames=evo_time, repeat=False,blit=False)


    plt.show()

    # #
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('C:\\Liang\\Googlebox\\Research\\Project2\\BaleenWhales\\BaleenWhaleFitnessDiffTT'
               '.mp4',writer=writer)
else:
    print('Junk simulation! Please try again or elevate K.')


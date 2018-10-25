import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_py import DVSim
from dvtraitsim_shared import DVTreeData, DVParam
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
theta = 0  # optimum of natural selection
gamma1 = 0.001  # intensity of natural selection
r = 1  # growth rate
a = 0.1  # intensity of competition
K = 10e8  # carrying capacity
kscale=100000
delta_pop = .001  # Variance of random walk of population
nu=1/(100*K)
Vmax = 1
scalar = 1000
no_tree = 1


# trait evolution plot
tree = 'tree'+'%d' % no_tree
example = 'example'+'%d' % no_tree
if platform.system()=='Windows':
    dir_path = 'c:/Liang/Googlebox/Research/Project2'
    files = dir_path + '/treesim_newexp/'+example+'/'
    td = DVTreeData(path=files, scalar=scalar)
elif platform.system()=='Darwin':
    file = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/tree_data/'+example+'/'

# parameter settings
obs_param = DVParam(gamma=gamma1, a=a, K=K, nu=nu, r=r, theta=theta, Vmax=1, inittrait=0, initpop=500,
                    initpop_sigma=10.0, break_on_mu=False)

for rep in range(100):
    simresult = DVSim(td,obs_param)
    if simresult['sim_time'] == td.sim_evo_time:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)


if simresult['sim_time'] == td.sim_evo_time:
    evo_time, total_species = simresult['N'].shape
    evo_time = evo_time-1
    trait_RI_dr = simresult['Z']
    population_RI_dr = simresult['N']


    trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
    population_tips = population_RI_dr[evo_time,:][~np.isnan(population_RI_dr[evo_time,:])]

    # trait_RI_dr[np.where(trait_RI_dr == 0)[0],np.where(trait_RI_dr == 0)[1]] = np.nan

    # population_RI_dr[np.where(population_RI_dr == 0)[0],np.where(population_RI_dr == 0)[1]] = np.nan
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
    timelist = np.genfromtxt(files + 'timelist.csv', delimiter=',')
    timebranch = np.genfromtxt(files + 'timebranch.csv', delimiter=',')
    timeend = np.genfromtxt(files + 'timeend.csv', delimiter=',')
    traittable = np.genfromtxt(files + 'traittable.csv', delimiter=',')
    ltable = np.genfromtxt(files + 'Ltable.csv', delimiter=',')
    # processing data
    ltable = np.delete(ltable, (0), axis=0)
    ltable = np.delete(ltable, (0), axis=1)
    traittable = np.delete(traittable, (0), axis=0)
    traittable = np.delete(traittable, (0), axis=1)
    daughter_index = np.absolute(ltable[:, 2])
    daughter_index = [int(x) for x in daughter_index]

    parent_index = np.absolute(ltable[:, 1])
    parent_index = [int(x) for x in parent_index]

    evo_timelist = timelist[:, 1]
    evo_timelist = np.delete(evo_timelist, 0)
    evo_timelist = max(evo_timelist) - evo_timelist
    evo_timelist = evo_timelist * scalar
    evo_timelist = evo_timelist.astype(int)

    timebranch = timebranch[:, 1]
    timebranch = np.delete(timebranch, 0)
    timebranch = [int(x) - 1 for x in timebranch]
    timeend = timeend[:, 1]
    timeend = np.delete(timeend, 0)
    timeend = [int(x) - 1 for x in timeend]

    # evolution time: speciation time
    evo_time = max(evo_timelist)

    speciate_time = evo_timelist[timebranch]
    extinct_time = evo_timelist[timeend]
    extinct_time[np.where(extinct_time == evo_time)[0]] = evo_time+10

    # extinct_time = np.delete(extinct_time, np.where(extinct_time == evo_time)[0])
    total_species = len(speciate_time)



    trait_RI_dr = simresult['Z']
    population_RI_dr = simresult['N']
    ext_times_RI = []
    pop_nan20 = np.nan_to_num(population_RI_dr)
    pop_sum = np.sum(pop_nan20,axis=1)
    # evo_time, total_species = simresult[0].shape
    # evo_time = evo_time-1
    ext_spec_index_RI = np.where(population_RI_dr[evo_time,] == 0)[0]
    # ext_spec_index_RI_rs = np.where(population_RI_rs[evo_time,] == 0)[0]


    for j in ext_spec_index_RI:
        ext_time = np.where(population_RI_dr[:,j] == 0)
        ext_times_RI.append(ext_time[0][0])


    x = np.array(range(evo_time))
    RI_traits = []
    RI_sizes = []

    for i in np.arange(total_species):
        RI_trait = trait_RI_dr[x,i]
        RI_size = population_RI_dr[x,i]
        RI_traits.append(RI_trait)
        RI_sizes.append(RI_size)

    f0 = plt.figure(num = 0, figsize = (12, 8))#, dpi = 100)
    f0.suptitle("Trait Evolution", fontsize=12)
    ax01 = plt.subplot2grid((1, 1), (0, 0))


    i = theta
    ax01.axhline(y = i, color='k', linestyle='--',alpha=0.7)
    ax01.set_title('Animation of the trait-population coevolution model')

    trait_max = np.nanmax(trait_RI_dr)
    trait_min = np.nanmin(trait_RI_dr)

    # set y-limits
    ax01.set_ylim(trait_min-10,trait_max+10)
    ax01.set_xlim(0, evo_time)

    # Turn on grids
    ax01.grid(True)
    ax01.set_xlabel("")
    ax01.set_ylabel("Trait Value")
    ax01.set_xlabel("Generation")
    ax01.set_ylabel("Trait Value")
    ax02 = plt.twinx()
    ax02.set_ylabel("Total abundance")
    popu_RI_spec_texts= []


    text_y1 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)

    for i in np.arange(total_species):
        popu_RI_spec_text1 = ax01.text(1,text_y1[i],'', transform = ax01.transAxes)
        popu_RI_spec_texts.append(popu_RI_spec_text1)

    RI_lines = []
    RI_scatters = []
    pop_line = []
    for i in  np.arange(total_species):
        RI_line, = ax01.plot([],[], 'b-')
        RI_scatter = ax01.scatter([], [], s=0, c='r', alpha=0.3)
        RI_lines.append(RI_line)
        RI_scatters.append(RI_scatter)
    # pop_line, = ax02.plot([],[],'k--')


    def animate(i):
        i = (i+1)%(len(x)+1)
        time_text.set_text(time_template % (i))
        RI_datas = []
        # pop_line.set_data(x[0:i], pop_sum[0:i])
        ax02.plot(x[0:i:20], pop_sum[0:i:20],'k-')
        for j in np.arange(total_species):
            RI_data = np.hstack((x[i], RI_traits[j][i]))
            RI_datas.append(RI_data)
            RI_lines[j].set_data(x[0:i], RI_traits[j][0:i])
            RI_scatters[j].set_offsets(RI_datas[j])
            RI_scatters[j].set_sizes([RI_sizes[j][i]/kscale])


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


        return  RI_lines, RI_scatters,pop_line


    ##
    time_template = 'Time = %d G'    # prints running simulation time
    time_text = ax01.text(0.05, 1, '', transform=ax01.transAxes)
    #
    ani = animation.FuncAnimation(f0, animate, interval= 1, frames= evo_time, repeat=False, blit=False) #, init_func=init)
    plt.show()
    #
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=2, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('C:\\Liang\\DVmodeltree1_popsum1q.mp4', writer=writer)
else:
    print('Junk simulation! Please try again or elevate K.')
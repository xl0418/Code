import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/Code/Pro2/Python_p2')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from DV_model_sim_along_phy import DVtraitsim_tree
import seaborn as sns
from matplotlib.pylab import *
from matplotlib import animation

theta = 0  # optimum of natural selection
gamma1 = 0.00  # intensity of natural selection
r = 1  # growth rate
a = 0.00  # intensity of competition
K = 1000000000  # carrying capacity
kscale=1000
delta_pop = .001  # Variance of random walk of population
nu=1/(100*K)
Vmax = 1
scalar = 1000



# trait evolution plot
if platform.system()=='Windows':
    file = 'C:\\Liang\\Code\\Pro2\\abcpp\\tree_data\\example3\\'
elif platform.system()=='Darwin':
    file = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/tree_data/example3/'

for rep in range(100):
    simresult = DVtraitsim_tree(file = file, gamma1 = gamma1, a = a,K=K, scalar=scalar,replicate=1)
    if simresult[2]:
        break
    else:
        print('%d simulations are all junks! Try more!' % rep)


if simresult[2]:
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
    timelist = np.genfromtxt(file + 'timelist.csv', delimiter=',')
    timebranch = np.genfromtxt(file + 'timebranch.csv', delimiter=',')
    timeend = np.genfromtxt(file + 'timeend.csv', delimiter=',')
    traittable = np.genfromtxt(file + 'traittable.csv', delimiter=',')
    ltable = np.genfromtxt(file + 'Ltable.csv', delimiter=',')
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



    trait_RI_dr = simresult[0]
    population_RI_dr = simresult[1]
    ext_times_RI = []

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

    f0 = figure(num = 0, figsize = (12, 8))#, dpi = 100)
    f0.suptitle("Trait Evolution", fontsize=12)
    ax01 = subplot2grid((1, 1), (0, 0))


    i = theta
    ax01.axhline(y = i, color='k', linestyle='--',alpha=0.7)
    ax01.set_title('RI-sd')

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

    popu_RI_spec_texts= []


    text_y1 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)

    for i in np.arange(total_species):
        # popu_rs_spec_text1 = ax02.text(1,text_y2[i], '', transform = ax02.transAxes)
        popu_RI_spec_text1 = ax01.text(1,text_y1[i],'', transform = ax01.transAxes)
        # popu_rs_spec_texts.append(popu_rs_spec_text1)
        popu_RI_spec_texts.append(popu_RI_spec_text1)

    RI_lines = []
    RI_scatters = []

    for i in  np.arange(total_species):
        RI_line, = ax01.plot([],[], 'b-')
        # rs_line, = ax02.plot([],[], 'b-')

        RI_scatter = ax01.scatter([], [], s=0, c='r', alpha=0.3)
        # rs_scatter = ax02.scatter([], [], s=0, c='r', alpha=0.3)

        RI_lines.append(RI_line)
        RI_scatters.append(RI_scatter)
        # rs_lines.append(rs_line)
        # rs_scatters.append(rs_scatter)

    def animate(i):
        i = (i+1)%(len(x)+1)
        time_text.set_text(time_template % (i))
        RI_datas = []
        for j in np.arange(total_species):
            # rs_data = np.hstack((x[i], RI_traits_rs[j][i]))
            RI_data = np.hstack((x[i], RI_traits[j][i]))
            # rs_datas.append(rs_data)
            RI_datas.append(RI_data)

        # for j in np.arange(total_species):
        #     rs_lines[j].set_data(x[0:i], RI_traits_rs[j][0:i])
            RI_lines[j].set_data(x[0:i], RI_traits[j][0:i])
            # rs_scatters[j].set_offsets(rs_datas[j])
            RI_scatters[j].set_offsets(RI_datas[j])
            # rs_scatters[j].set_sizes([RI_sizes_rs[j][i]])
            RI_scatters[j].set_sizes([RI_sizes[j][i]/kscale])
            # BH_lines[j].set_label("spec %d" % j)
            # Extinct species being labeled by dashed lines
            # if (population_BH[i, j] == 0):
            #     BH_lines[j].set_dashes([2, 2, 2, 2])
            #     BH_lines[j].set_color("red")
            # if (population_RI[i, j] == 0):
            #     RI_lines[j].set_dashes([2, 2, 2, 2])
            #     RI_lines[j].set_color("red")

            # Animating labels
            popu_RI_spec_texts[j].set_text('POS %d = %.1f' % (j+1, population_RI_dr[i,j]))
            # popu_rs_spec_texts[j].set_text('POS %d = %.1f' % (j+1, population_RI_rs[i,j]))
            if (i < speciate_time[j]):
                # rs_lines[2].set_data([], [])
                RI_lines[j].set_data([], [])
                # rs_lines[1].set_data(x[0:i], RI_traits_rs[1][0:i])
                # RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
            elif (i >= speciate_time[j] and i < extinct_time[j]):
                # rs_lines[2].set_data(x[1000:i], RI_traits_rs[2][1000:i])
                RI_lines[j].set_data(x[speciate_time[j]:i], RI_traits[j][speciate_time[j]:i])
                # rs_lines[1].set_data(x[0:i], RI_traits_rs[1][0:i])
                # RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
            else:
                # rs_lines[2].set_data(x[1000:i], RI_traits_rs[2][1000:i])
                RI_lines[j].set_data(x[speciate_time[j]:extinct_time[j]],
                                        RI_traits[j][speciate_time[j]:extinct_time[j]])
                # rs_lines[1].set_data(x[0:2000], RI_traits_rs[1][0:2000])
                # RI_lines[1].set_data(x[0:2000], RI_traits[1][0:2000])

            if (j in ext_spec_index_RI):
                end_time_RI = ext_times_RI[np.where(j == ext_spec_index_RI)[0][0]]
                if (i >= end_time_RI):
                    RI_lines[j].set_data(x[0:end_time_RI], RI_traits[j][0:end_time_RI])

        # RI_lines[1].set_color("green")
        # RI_lines[2].set_color("green")

        return  RI_lines, RI_scatters


    ##
    time_template = 'Time = %d G'    # prints running simulation time
    time_text = ax01.text(0.05, 1, '', transform=ax01.transAxes)
    #
    ani = animation.FuncAnimation(f0, animate, interval= 1, frames= evo_time, repeat=False, blit=False) #, init_func=init)
    plt.show()
    #
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('C:\\Liang\\DVmodeltree4_n.mp4', writer=writer)
else:
    print('Junk simulation! Please try again or elevate K.')
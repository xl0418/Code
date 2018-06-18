import matplotlib.animation as animation
from matplotlib.pylab import *
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(12)
theta =  0   # optimum of natural selection
gamma1 = 0.001 # intensity of natural selection
r = 2  # growth rate
a = 0.1 # intensity of competition
K = 5000  # carrying capacity
nu = 0.0001
Vmax = 1

speciate_time = 1000
extinction_time = 2000
evo_time = 3000
total_species = 10

# Natural selection function
def ga(gamma, theta, zi):
    return -gamma * (theta - zi) ** 2

# Dynamic carrying capacity
def Kd(gamma_K, theta, zi, K):
    return max(K * np.exp(-gamma_K * (theta - zi) ** 2),1)


# Competition function
def beta(a, zi, zj, nj):
    zi_ret = np.ndarray((1,len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0,n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Derivative of the competition function
def sigma(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(2 * a * (zi[n1]-np.array(zj)) * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Derivative of the competition function
def sigmasqr(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(4 * a**2 * (zi[n1]-np.array(zj))**2 * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Delta function
def Delta(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(-2 * a +4* a **2  * (zi[n1]-np.array(zj))**2 * np.exp( -a * (zi[n1] - np.array(zj)) ** 2)
                               * np.array(nj))
    return zi_ret


# Variance of random walk of trait evolution
delta_trait = 0.1
# build trait evolution and population evolution matrices
trait_RI_dr = np.zeros((evo_time + 1, total_species))
population_RI_dr = np.zeros((evo_time + 1, total_species))

V = np.zeros((evo_time + 1, total_species))
V.fill(1 / total_species)

trait_RI_rs = np.zeros((evo_time + 1, total_species))
population_RI_rs = np.zeros((evo_time + 1, total_species))

V_rs = np.zeros((evo_time + 1, total_species))
V_rs.fill(1 / total_species)


existing_species = np.zeros(shape = (3,total_species))  # np.matrix([[1,1,0], [1,1,1],[1,0,1]])
existing_species.fill(1)
existing_species[0,2]=existing_species[2,1] = 0

#  initialize condition for species trait and population
mu_trait, sigma_trait = 0, 10  # mean and standard deviation
trait_RI_dr[0] = 0   #existing_species[0] * np.random.normal(mu_trait, sigma_trait, total_species)
trait_RI_rs[0] = 0
print(trait_RI_dr[0])
mu_pop, sigma_pop = 50, 10  # mean and standard deviation
population_RI_dr[0] = existing_species[0] * np.random.normal(mu_pop, sigma_pop, total_species)
population_RI_rs[0] = population_RI_dr[0]
print(population_RI_dr[0])

# vectorize function ga
# ga_vector = np.vectorize(ga)
Kd_vector = np.vectorize(Kd)


for i in range(evo_time):
    delta_pop = 0.01
    # parameter controls the type of competition: 1, competition; -1 cooperation.
    m = 1
    if(i < speciate_time):
        node = 0
    elif(i < extinction_time and i>= speciate_time):
        node = 1
    else:
        node = 2
    # num_species = node + 2
    index_existing_species = np.where(existing_species[node] == 1)[0]
    K_RI_dr = K
    Beta_RI_dr = beta(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species],
                      nj=population_RI_dr[i,index_existing_species])
    Sigma_RI_dr = sigma(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species], 
                        nj=population_RI_dr[i,index_existing_species])
    Sigmasqr_RI_dr = sigmasqr(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species],
                              nj=population_RI_dr[i,index_existing_species])
    Delta_RI_dr = Delta(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species],
                        nj=population_RI_dr[i,index_existing_species])
    var_trait = V[i,index_existing_species] / (2 * population_RI_dr[i,index_existing_species])
    trait_RI_dr[i + 1,index_existing_species] = trait_RI_dr[i,index_existing_species] + V[i,index_existing_species] *\
                                                (2 * gamma1 * (theta -
                                                    trait_RI_dr[i,index_existing_species]) + 1 / K_RI_dr *
                                                  Sigma_RI_dr) \
                         + np.random.normal(0, var_trait, len(index_existing_species))
    population_RI_dr[i + 1,index_existing_species] = population_RI_dr[i,index_existing_species] * r * \
                              np.exp(-gamma1 * (theta - trait_RI_dr[i,index_existing_species]) ** 2 + 
                                     (1 - Beta_RI_dr / K_RI_dr) \
                                     + np.random.normal(0, delta_pop, len(index_existing_species)))
    V[i + 1,index_existing_species] = V[i,index_existing_species]+ V[i,index_existing_species] ** 2 * (-2 * gamma1 + 4 * gamma1 ** 2 * (theta - trait_RI_dr[i,index_existing_species]) ** 2 +
                            1 / K_RI_dr *  # V^2 w''/w
                            (Sigma_RI_dr - Sigmasqr_RI_dr) + 4 * gamma1 / K_RI_dr * (
                                        theta - trait_RI_dr[i,index_existing_species]) * Sigma_RI_dr +
                            Sigma_RI_dr ** 2 / K_RI_dr ** 2) - \
               V[i,index_existing_species] / 2 + 2 * population_RI_dr[i,index_existing_species] * nu * Vmax / (
                           1 + 4 * population_RI_dr[i,index_existing_species] * nu)  # loss due to sexual selection; gain from mutation

    #
    # K_RI_rs = K
    # Beta_RI_rs = beta(a=a, zi=trait_RI_rs[i, index_existing_species], zj=trait_RI_rs[i, index_existing_species],
    #                   nj=population_RI_rs[i, index_existing_species])
    # Sigma_RI_rs = sigma(a=a, zi=trait_RI_rs[i, index_existing_species], zj=trait_RI_rs[i, index_existing_species],
    #                     nj=population_RI_rs[i, index_existing_species])
    # Sigmasqr_RI_rs = sigmasqr(a=a, zi=trait_RI_rs[i, index_existing_species], zj=trait_RI_rs[i, index_existing_species],
    #                           nj=population_RI_rs[i, index_existing_species])
    # Delta_RI_rs = Delta(a=a, zi=trait_RI_rs[i, index_existing_species], zj=trait_RI_rs[i, index_existing_species],
    #                     nj=population_RI_rs[i, index_existing_species])
    # var_trait_rs = V_rs[i, index_existing_species] / (2 * population_RI_rs[i, index_existing_species])
    # trait_RI_rs[i + 1, index_existing_species] = trait_RI_rs[i, index_existing_species] + V[i, index_existing_species] * \
    #                                              (2 * gamma1 * (theta -
    #                                                             trait_RI_rs[i, index_existing_species]) + 1 / K_RI_rs *
    #                                               Sigma_RI_rs) \
    #                                              + np.random.normal(0, var_trait, len(index_existing_species))
    # population_RI_rs[i + 1, index_existing_species] = population_RI_rs[i, index_existing_species] * r * \
    #                                                   np.exp(-gamma1 * (
    #                                                               theta - trait_RI_rs[i, index_existing_species]) ** 2 +
    #                                                          (1 - Beta_RI_rs / K_RI_rs) \
    #                                                          + np.random.normal(0, delta_pop,
    #                                                                             len(index_existing_species)))
    # V_rs[i + 1, index_existing_species] = V_rs[i, index_existing_species] + V_rs[i, index_existing_species] ** 2 * (
    #             -2 * gamma1 + 4 * gamma1 ** 2 * (theta - trait_RI_rs[i, index_existing_species]) ** 2 +
    #             1 / K_RI_rs *  # V^2 w''/w
    #             (Sigma_RI_rs - Sigmasqr_RI_rs) + 4 * gamma1 / K_RI_rs * (
    #                     theta - trait_RI_rs[i, index_existing_species]) * Sigma_RI_rs +
    #             Sigma_RI_rs ** 2 / K_RI_rs** 2) - \
    #                                    V_rs[i, index_existing_species] / 2 + 2 * population_RI_rs[
    #                                        i, index_existing_species] * nu * Vmax / (
    #                                            1 + 4 * population_RI_rs[
    #                                        i, index_existing_species] * nu)  # loss due to sexual selection; gain from mutation

    if (i+1) == speciate_time:
         trait_RI_dr[i+1,2] = trait_RI_dr[i+1,1] #+ np.random.normal(0, 0.01, 1)
         population_RI_dr[i + 1, 2] =1/2 * population_RI_dr[i + 1, 1]
         population_RI_dr[i + 1, 1] = 1 / 2 * population_RI_dr[i + 1, 1]
         V[i+1,2] = 1/2 * V[i+1,1]
         V[i+1,1] = 1/2 * V[i+1,1]
         #
         # trait_RI_rs[i + 1, 2] = trait_RI_rs[i + 1, 1]# + np.random.normal(0, 0.01, 1)
         # population_RI_rs[i + 1, 2] = 1 / 2 * population_RI_rs[i + 1, 1]
         # population_RI_rs[i + 1, 1] = 1 / 2 * population_RI_rs[i + 1, 1]
         # V_rs[i + 1, 2] = 1 / 2 * V_rs[i + 1, 1]
         # V_rs[i + 1, 1] = 1 / 2 * V_rs[i + 1, 1]

    if (i + 1) == extinction_time:
        trait_RI_dr[i+1, 1] = None
        population_RI_dr[i+1, 1] = 0
        # trait_RI_rs[i + 1, 1] = None
        # population_RI_rs[i + 1, 1] = 0




ext_times_RI = []
# ext_times_RI_rs = []

ext_spec_index_RI = np.where(population_RI_dr[evo_time,] == 0)[0]
# ext_spec_index_RI_rs = np.where(population_RI_rs[evo_time,] == 0)[0]


for j in ext_spec_index_RI:
    ext_time = np.where(population_RI_dr[:,j] == 0)
    ext_times_RI.append(ext_time[0][0])

#
# for j in ext_spec_index_RI_rs:
#     ext_time_rs = np.where(population_RI_rs[:,j] == 0)
#     ext_times_RI_rs.append(ext_time_rs[0][0])

x = np.array(range(evo_time))
RI_traits = []
RI_sizes = []
# RI_traits_rs = []
# RI_sizes_rs = []
for i in np.arange(total_species):
    RI_trait = trait_RI_dr[x,i]
    RI_size = population_RI_dr[x,i]
    RI_traits.append(RI_trait)
    RI_sizes.append(RI_size)

    # RI_trait_rs = trait_RI_rs[x,i]
    # RI_size_rs = population_RI_rs[x,i]
    # RI_traits_rs.append(RI_trait_rs)
    # RI_sizes_rs.append(RI_size_rs)


f0 = figure(num = 0, figsize = (12, 8))#, dpi = 100)
f0.suptitle("Trait Evolution", fontsize=12)
ax01 = subplot2grid((1, 1), (0, 0))
# ax02 = subplot2grid((2, 1), (1, 0))


i = 0
ax01.axhline(y = i, color='k', linestyle='--',alpha=0.7)
# ax02.axhline(y = i, color='k', linestyle='--',alpha=0.7)




ax01.set_title('RI-sd')
# ax02.set_title('RI-rs')


trait_max = max(np.nanmax(trait_RI_dr),np.nanmax(trait_RI_rs))
trait_min = min(np.nanmin(trait_RI_dr),np.nanmin(trait_RI_rs))


# set y-limits
ax01.set_ylim(trait_min-10,trait_max+10)
# ax02.set_ylim(trait_min-10,trait_max+10)

# ax03.set_ylim(-0,5)
# ax04.set_ylim(-10,10)

# sex x-limits
ax01.set_xlim(0, evo_time)
# ax02.set_xlim(0, evo_time)

# ax03.set_xlim(0,5.0)
# ax04.set_xlim(0,5.0)

# Turn on grids
ax01.grid(True)
ax02.grid(True)


#
ax01.set_xlabel("")
ax01.set_ylabel("Trait Value")
ax01.set_xlabel("Generation")
ax01.set_ylabel("Trait Value")

popu_RI_spec_texts= []

popu_rs_spec_texts= []


text_y1 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)
text_y2 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)

for i in np.arange(total_species):
    # popu_rs_spec_text1 = ax02.text(1,text_y2[i], '', transform = ax02.transAxes)
    popu_RI_spec_text1 = ax01.text(1,text_y1[i],'', transform = ax01.transAxes)
    # popu_rs_spec_texts.append(popu_rs_spec_text1)
    popu_RI_spec_texts.append(popu_RI_spec_text1)

RI_lines = []
RI_scatters = []
rs_lines = []
rs_scatters = []


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
    rs_datas = []
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
        RI_scatters[j].set_sizes([RI_sizes[j][i]])
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

        # if (1 not in ext_spec_index_BH):
        if (i < 1000):
            # rs_lines[2].set_data([], [])
            RI_lines[2].set_data([], [])
            # rs_lines[1].set_data(x[0:i], RI_traits_rs[1][0:i])
            RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
        elif (i >= 1000 and i < 2000):
            # rs_lines[2].set_data(x[1000:i], RI_traits_rs[2][1000:i])
            RI_lines[2].set_data(x[1000:i], RI_traits[2][1000:i])
            # rs_lines[1].set_data(x[0:i], RI_traits_rs[1][0:i])
            RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
        else:
            # rs_lines[2].set_data(x[1000:i], RI_traits_rs[2][1000:i])
            RI_lines[2].set_data(x[1000:i], RI_traits[2][1000:i])
            # rs_lines[1].set_data(x[0:2000], RI_traits_rs[1][0:2000])
            RI_lines[1].set_data(x[0:2000], RI_traits[1][0:2000])

        # if (j in ext_spec_index_RI_rs):
        #     end_time_BH = ext_times_RI_rs[np.where(j == ext_spec_index_RI_rs)[0][0]]
        #     if (i >= end_time_BH):
        #         rs_lines[j].set_data(x[0:end_time_BH], RI_traits_rs[j][0:end_time_BH])

        if (j in ext_spec_index_RI):
            end_time_RI = ext_times_RI[np.where(j == ext_spec_index_RI)[0][0]]
            if (i >= end_time_RI):
                RI_lines[j].set_data(x[0:end_time_RI], RI_traits[j][0:end_time_RI])

    # rs_lines[1].set_color("green")
    # rs_lines[2].set_color("green")
    RI_lines[1].set_color("green")
    RI_lines[2].set_color("green")

    return  RI_lines, RI_scatters


##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 1, '', transform=ax01.transAxes)
#
ani = animation.FuncAnimation(f0, animate, interval= 1, frames= 3000, repeat=False, blit=False) #, init_func=init)
plt.show()
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save('/Users/dudupig/Google 云端硬盘/Python/Project2/S+C.mp4', writer=writer)
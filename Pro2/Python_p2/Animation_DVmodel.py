import matplotlib.animation as animation
from matplotlib.pylab import *
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(12)
theta =  0   # optimum of natural selection
gamma1 = 0.01 # intensity of natural selection
r = 1  # growth rate
a = 0.05 # intensity of competition
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
trait_RI_dk_dk = np.zeros((evo_time + 1, total_species))
population_RI_dk_dk = np.zeros((evo_time + 1, total_species))
V = np.zeros((evo_time + 1, total_species))
V.fill(1 / total_species)

#  initialize condition for species trait and population
mu_trait, sigma_trait = 0, 10  # mean and standard deviation
trait_RI_dr[0, (0,1)] = np.random.normal(mu_trait, sigma_trait, 2)
trait_RI_dk_dk[0] = trait_RI_dr[0]
print(trait_RI_dr[0])
mu_pop, sigma_pop = 50, 10  # mean and standard deviation
population_RI_dr[0, (0,1)] = np.random.normal(mu_pop, sigma_pop, 2)
population_RI_dk_dk[0] = population_RI_dr[0]
print(population_RI_dr[0])

# vectorize function ga
# ga_vector = np.vectorize(ga)
Kd_vector = np.vectorize(Kd)
# Existing species matrix
existing_species = np.matrix([[1,1,0], [1,1,1],[1,0,1]])
speciation_event = np.array([1])


for i in range(evo_time):
    delta_pop = 0.1
    # parameter controls the type of competition: 1, competition; -1 cooperation.
    m = 1
    if(i < speciate_time):
        node = 0
    elif(i < extinction_time and i>= speciate_time):
        node = 1
    else:
        node = 2
    # num_species = node + 2
    index_existing_species = np.where(existing_species[node] == 1)[1]
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

  
    if (i+1) == speciate_time:
         trait_RI_dr[i+1,2] = trait_RI_dr[i+1,1] + np.random.normal(0, 0.01, 1)
         population_RI_dr[i + 1, 2] =1/2 * population_RI_dr[i + 1, 1]
         population_RI_dr[i + 1, 1] = 1 / 2 * population_RI_dr[i + 1, 1]
         V[i+1,2] = 1/2 * V[i+1,1]
         V[i+1,1] = 1/2 * V[i+1,1]

         # trait_RI_dk[i + 1, 2] = trait_RI_dk[i + 1, 1] + np.random.normal(0, 0.01, 1)
         # population_RI_dk[i + 1, 2] = 1 / 2 * population_RI_dk[i + 1, 1]
         # population_RI_dk[i + 1, 1] = 1 / 2 * population_RI_dk[i + 1, 1]
    if (i + 1) == extinction_time:
        trait_RI_dr[i+1, 1] = None
        population_RI_dr[i+1, 1] = 0
        # trait_RI_dk[i + 1, 1] = None
        # population_RI_dk[i + 1, 1] = 0





x = np.array(range(evo_time))
BH_trait1 = trait_RI_dr[x,0]
BH_trait2 = trait_RI_dr[x,1]
BH_trait3 = trait_RI_dr[x,2]
#
# RI_trait1 = trait_RI_dk[x,0]
# RI_trait2 = trait_RI_dk[x,1]
# RI_trait3 = trait_RI_dk[x,2]

BH_size1 = population_RI_dr[x,0]
BH_size2 = population_RI_dr[x,1]
BH_size3 = population_RI_dr[x,2]
#
# RI_size1 = population_RI_dk_dk[x,0]
# RI_size2 = population_RI_dk[x,1]
# RI_size3 = population_RI_dk[x,2]
#


f0 = figure(num = 0, figsize = (12, 8))#, dpi = 100)
f0.suptitle("Trait Evolution", fontsize=12)
ax01 = subplot2grid((1, 1), (0, 0))
# ax02 = subplot2grid((2, 1), (1, 0))
# ax03 = subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)
# ax04 = ax03.twinx()

ax01.set_title('RI: dynamical r, constant K')
# ax02.set_title('RI: dynamic r, constant K')


# set y-limits
ax01.set_ylim(-40,40)
# ax02.set_ylim(-40,40)
# ax03.set_ylim(-0,5)
# ax04.set_ylim(-10,10)

# sex x-limits
ax01.set_xlim(0, evo_time)
# ax02.set_xlim(0, evo_time)
# ax03.set_xlim(0,5.0)
# ax04.set_xlim(0,5.0)

# Turn on grids
ax01.grid(True)
# ax02.grid(True)
# ax03.grid(True)
#
ax01.set_xlabel("")
ax01.set_ylabel("Trait Value")
# ax02.set_xlabel("Generation")
# ax02.set_ylabel("Trait Value")

popu_BH_spec1_text = ax01.text(0.8,0.9, '', transform = ax01.transAxes)
popu_BH_spec2_text = ax01.text(0.8,0.85, '', transform = ax01.transAxes)
popu_BH_spec3_text = ax01.text(0.8,0.8, '', transform = ax01.transAxes)

#
# popu_RI_spec1_text = ax02.text(0.8,0.4,'', transform = ax02.transAxes)
# popu_RI_spec2_text = ax02.text(0.8,0.35, '', transform = ax02.transAxes)
# popu_RI_spec3_text = ax02.text(0.8,0.3, '' ,transform = ax02.transAxes)



# ax03.set_xlabel("t")
# ax03.set_ylabel("py")
# ax04.set_ylabel("vy")
# ax = plt.axes(xlim = (0,2001), ylim = (-20,20))
# ax2 = plt.axes(xlim = (0,2001), ylim = (-20,20))

BH_line1, = ax01.plot([],[], 'b-')
BH_line2, = ax01.plot([],[], 'r-')
BH_line3, = ax01.plot([],[], 'g-')
#
# RI_line1, = ax02.plot([],[], 'b-')
# RI_line2, = ax02.plot([],[], 'r-')
# RI_line3, = ax02.plot([],[], 'g-')

BH_scatter1 = ax01.scatter([],[], s = 0, c = 'b', alpha = 0.5)
BH_scatter2 = ax01.scatter([],[], s = 0, c = 'r', alpha = 0.5)
BH_scatter3 = ax01.scatter([],[], s = 0, c = 'g', alpha = 0.5)
#
#
# RI_scatter1 = ax02.scatter([],[], s = 0, c = 'b', alpha = 0.5)
# RI_scatter2 = ax02.scatter([],[], s = 0, c = 'r', alpha = 0.5)
# RI_scatter3 = ax02.scatter([],[], s = 0, c = 'g', alpha = 0.5)
# #
# def init():
#     BH_line1.set_data([], [])
#     BH_line2.set_data([], [])
#     BH_line3.set_data([], [])
#     RI_line1.set_data([], [])
#     RI_line2.set_data([], [])
#     RI_line3.set_data([], [])
#     BH_scatter1.set_offsets([])
#     BH_scatter2.set_offsets([])
#     BH_scatter3.set_offsets([])
#     RI_scatter1.set_offsets([])
#     RI_scatter2.set_offsets([])
#     RI_scatter3.set_offsets([])
#     return BH_line1, BH_line2, BH_line3, RI_line1, RI_line2,RI_line3,BH_scatter1,BH_scatter2,BH_scatter3,RI_scatter1 \
# , RI_scatter2, RI_scatter3

def animate(i):
    i = (i+1)%(len(x)+1)
    time_text.set_text(time_template % (i))
    BH_data1 = np.hstack((x[i], BH_trait1[i]))
    BH_data2 = np.hstack((x[i], BH_trait2[i]))
    BH_data3 = np.hstack((x[i], BH_trait3[i]))
    # RI_data1 = np.hstack((x[i], RI_trait1[i]))
    # RI_data2 = np.hstack((x[i], RI_trait2[i]))
    # RI_data3 = np.hstack((x[i], RI_trait3[i]))
    BH_line1.set_data(x[0:i], BH_trait1[0:i])

    # RI_line1.set_data(x[0:i], RI_trait1[0:i])

    BH_scatter1.set_offsets(BH_data1)
    BH_scatter2.set_offsets(BH_data2)
    # RI_scatter1.set_offsets(RI_data1)
    # RI_scatter2.set_offsets(RI_data2)
    BH_scatter3.set_offsets(BH_data3)
    # RI_scatter3.set_offsets(RI_data3)
    BH_scatter1.set_sizes([BH_size1[i]])
    BH_scatter2.set_sizes([BH_size2[i]])
    # RI_scatter1.set_sizes([RI_size1[i]])
    # RI_scatter2.set_sizes([RI_size2[i]])
    BH_scatter3.set_sizes([BH_size3[i]])
    # RI_scatter3.set_sizes([RI_size3[i]])

    # Animating labels
    popu_BH_spec1_text.set_text('Pop of Spec 1 = %.1f' % population_RI_dr[i,0])
    popu_BH_spec2_text.set_text('Pop of Spec 2 = %.1f' % population_RI_dr[i,1])
    popu_BH_spec3_text.set_text('Pop of Spec 3 = %.1f' % population_RI_dr[i,2])
    #
    # popu_RI_spec1_text.set_text('Pop of Spec 1 = %.1f' % population_RI_dk[i,0])
    # popu_RI_spec2_text.set_text('Pop of Spec 2 = %.1f' % population_RI_dk[i,1])
    # popu_RI_spec3_text.set_text('Pop of Spec 3 = %.1f' % population_RI_dk[i,2])

    if (i < speciate_time):
        BH_line3.set_data([], [])
        # RI_line3.set_data([], [])
        BH_line2.set_data(x[0:i], BH_trait2[0:i])
        # RI_line2.set_data(x[0:i], RI_trait2[0:i])
    elif (i>= speciate_time and i< extinction_time):
        BH_line3.set_data(x[speciate_time:i], BH_trait3[speciate_time:i])
        # RI_line3.set_data(x[speciate_time:i], RI_trait3[speciate_time:i])
        BH_line2.set_data(x[0:i], BH_trait2[0:i])
        # RI_line2.set_data(x[0:i], RI_trait2[0:i])
    else:
        BH_line3.set_data(x[speciate_time:i], BH_trait3[speciate_time:i])
        # RI_line3.set_data(x[speciate_time:i], RI_trait3[speciate_time:i])
        BH_line2.set_data(x[0:extinction_time], BH_trait2[0:extinction_time])
        # RI_line2.set_data(x[0:extinction_time], RI_trait2[0:extinction_time])


    return BH_line1, BH_line2, BH_line3,BH_scatter1,BH_scatter2,BH_scatter3,#RI_line1, RI_line2,RI_line3,RI_scatter1 \
#, RI_scatter2, RI_scatter3

##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 0.9, '', transform=ax01.transAxes)



ani = animation.FuncAnimation(f0, animate, interval= 1, frames= 6000, repeat=False, blit=False) #, init_func=init)
plt.show()
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save('/Users/dudupig/Google 云端硬盘/Python/Project2/S+C.mp4', writer=writer)
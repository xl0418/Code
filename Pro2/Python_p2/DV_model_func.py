import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


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



# Trait simulation function under the new coevolution model with changing variance
def DVtraitsim_branch(num_time, num_species, num_iteration, gamma1, a, r, nu,Vmax, theta,K ,replicate = 1):
    j = 0   # initialize the iteration number
    delta_pop = 0.001
    num_vec = np.arange(1,(num_iteration+1),1) # iteration vector
    nrow = len(num_vec)    # Row number of the trait evolution history matrix
    stat_rate_trait_RI_dr = np.empty((nrow,num_species))   # trait evolution matrix under BH

    stat_rate_popu_RI_dr = np.empty((nrow,num_species))  # population evolution matrix under BH
    # vectorize ga function
    ga_vector = np.vectorize(ga)

    # loop for preset iteration
    for loop in num_vec:
        if replicate == 1:
            np.random.seed(loop)  # set random seed from 1 to max of iteration
        # print(j)
        # Matrices for trait and populations under BH and RI
        trait_RI_dr = np.zeros((num_time+1, num_species))
        population_RI_dr = np.zeros((num_time+1, num_species))
        V = np.zeros((num_time+1, num_species))

        # initialize the input for trait values and populations
        # mu_trait, sigma_trait = mean_trait, dev_trait  # mean and standard deviation
        # trait_RI_dr[0] = np.random.normal(mu_trait, sigma_trait, num_species)
        trait_RI_dr[0] = np.zeros( num_species)

        # mu_pop, sigma_pop = mean_pop, dev_pop  # mean and standard deviation
        # population_RI_dr[0] = np.random.normal(mu_pop, sigma_pop, num_species)
        pop_ini = np.empty(num_species)
        pop_ini.fill(100)
        population_RI_dr[0] = pop_ini

        V.fill(1/num_species)

        # trait evolution simulation
        for i in range(num_time):
            # print(i)
            # RI dynamic r model
            K_RI_dr = K
            Beta_RI_dr = beta(a=a, zi=trait_RI_dr[i], zj=trait_RI_dr[i], nj=population_RI_dr[i])
            Sigma_RI_dr = sigma(a=a, zi=trait_RI_dr[i], zj=trait_RI_dr[i], nj=population_RI_dr[i])
            Sigmasqr_RI_dr = sigmasqr(a=a, zi=trait_RI_dr[i], zj=trait_RI_dr[i], nj=population_RI_dr[i])
            Delta_RI_dr = Delta(a=a, zi=trait_RI_dr[i], zj=trait_RI_dr[i], nj=population_RI_dr[i])
            var_trait = V[i]/(2*population_RI_dr[i])
            trait_RI_dr[i + 1] = trait_RI_dr[i] + V[i]*(2*gamma1*(theta - trait_RI_dr[i])+1/K_RI_dr *
                                                        Sigma_RI_dr) \
                                +  np.random.normal(0, var_trait, num_species)
            population_RI_dr[i + 1] = population_RI_dr[i] * r * \
                                      np.exp(-gamma1*(theta - trait_RI_dr[i])**2 + (1 - Beta_RI_dr / K_RI_dr) \
                                      + np.random.normal(0, delta_pop, num_species))
            V[i+1] = V[i]+ V[i]**2 * (-2* gamma1 + 4* gamma1**2 * (theta - trait_RI_dr[i])**2+ 1/K_RI_dr *      #V^2 w''/w
                                (Sigma_RI_dr-Sigmasqr_RI_dr) +4*gamma1/K_RI_dr *(theta - trait_RI_dr[i])*Sigma_RI_dr +
                                Sigma_RI_dr**2/K_RI_dr**2) - \
             V[i]/2 + 2 * population_RI_dr[i] * nu * Vmax/(1+4*population_RI_dr[i]*nu)  #loss due to sexual selection; gain from mutation

            # population_RI_dr[i + 1, np.where(population_RI_dr[i + 1] < 1)] = 0
            ext_index_RI_dr = np.where(population_RI_dr[i +1] == 0)[0]
            if len(ext_index_RI_dr) > 0:
                trait_RI_dr[i+1][ext_index_RI_dr] = 0


        # Diversity statistics
        stat_rate_trait_RI_dr[j,:] =trait_RI_dr[num_time,:]
        stat_rate_popu_RI_dr[j, :] = population_RI_dr[num_time, :]
        j += 1

    return stat_rate_trait_RI_dr,stat_rate_popu_RI_dr,





# boxplots and trait distributions
def drawplot(traitdata):
    # read simulated data
    stat_rate_trait_RI_dr = traitdata[0]
    stat_rate_trait_RI_dk = traitdata[1]
    stat_rate_popu_RI_dr = traitdata[2]
    stat_rate_popu_RI_dk = traitdata[3]
    num_species = traitdata[0].shape[1]
    # find out extinct species and remove the responding trait values
    ext_index_RI_dr = np.where(stat_rate_popu_RI_dr == 0)
    ext_index_RI_dk = np.where(stat_rate_popu_RI_dk == 0)
    statplot_trait_RI_dr = stat_rate_trait_RI_dr
    statplot_trait_RI_dr[ext_index_RI_dr[0],ext_index_RI_dr[1]] = np.nan
    statplot_trait_RI_dr_sorted = np.sort(statplot_trait_RI_dr)
    statplot_trait_RI_dk = stat_rate_trait_RI_dk
    statplot_trait_RI_dk[ext_index_RI_dk[0],ext_index_RI_dk[1]] = np.nan
    statplot_trait_RI_dk_sorted = np.sort(statplot_trait_RI_dk)

    # filter the missing data
    mask_RI_dr = ~np.isnan(statplot_trait_RI_dr_sorted)
    filtered_data_RI_dr = [d[m] for d, m in zip(statplot_trait_RI_dr_sorted.T, mask_RI_dr.T)]
    mask_RI_dk = ~np.isnan(statplot_trait_RI_dk_sorted)
    filtered_data_RI_dk = [d[m] for d, m in zip(statplot_trait_RI_dk_sorted.T, mask_RI_dk.T)]
    # convert the list to an array
    merge_RI_dr = np.concatenate( filtered_data_RI_dr, axis=0 )
    merge_RI_dk = np.concatenate( filtered_data_RI_dk, axis=0 )
    # create a fig
    fig = plt.figure(1, figsize=(12, 9))
    # determine the limits for axises
    min_RI_dr = np.amin(merge_RI_dr)-5
    max_RI_dr = np.amax(merge_RI_dr)+5
    min_RI_dk = np.amin(merge_RI_dk)-5
    max_RI_dk = np.amax(merge_RI_dk)+5
    global_min = -50         # min(min_RI_dr, min_RI_dk)
    global_max = 50          # max(max_RI_dr, max_RI_dk)
    # multiple plots arrangement
    gs = gridspec.GridSpec(1, 4)

    # Create an axes instance
    ax1 = fig.add_subplot(gs[0,:-2])
    # boxplot BH trait data
    bh = ax1.boxplot(filtered_data_RI_dr, 0 , "", patch_artist=True)
    ## change outline color, fill color and linewidth of the boxes
    for box in bh['boxes']:
        # change outline color
        box.set( color='DarkBlue', linewidth=0.5)
        # change fill color
        box.set( facecolor = '#95d0fc' ,  alpha=0.5)

    ## change color and linewidth of the whiskers
    for whisker in bh['whiskers']:
        whisker.set(color='#95d0fc', linewidth=0.5)

    ## change color and linewidth of the caps
    for cap in bh['caps']:
        cap.set(color='#95d0fc', linewidth=0.5)

    ## change color and linewidth of the medians
    for median in bh['medians']:
        median.set(color='#95d0fc', linewidth=0.5)

    ## change the style of fliers and their fill
    for flier in bh['fliers']:
        flier.set(marker='o', color='#95d0fc', alpha=0.5)

    # boxplot RI trait data
    ri = ax1.boxplot(filtered_data_RI_dk, 0 , "", patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    for box in ri['boxes']:
        # change outline color
        box.set( color='#ff000d', linewidth=0.5)
        # change fill color
        box.set( facecolor = '#fc5a50' ,  alpha=0.5)

    ## change color and linewidth of the whiskers
    for whisker in ri['whiskers']:
        whisker.set(color='#fc5a50', linewidth=0.5)

    ## change color and linewidth of the caps
    for cap in ri['caps']:
        cap.set(color='#fc5a50', linewidth=0.5)

    ## change color and linewidth of the medians
    for median in ri['medians']:
        median.set(color='#fc5a50', linewidth=0.5)

    ## change the style of fliers and their fill
    for flier in ri['fliers']:
        flier.set(marker='o', color='#fc5a50', alpha=0.5)

    # add legends
    ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['DR', 'DK'], loc='upper left')
    # add title
    ax1.set_title('Trait distribution')
    # add x label and y label
    ax1.set_xlabel('Species ordered by trait values')
    ax1.set_ylabel('Trait values')
    ax1.set_ylim(global_min,global_max)

    # Major ticks, minor ticks for axises
    # major_ticks_x = np.arange(0, num_species+1, 20)
    # minor_ticks_x = np.arange(0, num_species+1, 1)
    major_ticks_x = np.arange(0, num_species + 1)
    minor_ticks_x = np.arange(0, num_species + 1)
    major_ticks_y = np.arange(global_min, global_max, 10)
    minor_ticks_y = np.arange(global_min, global_max, 1)

    ax1.set_xticks(major_ticks_x)
    ax1.set_xticklabels(major_ticks_x)
    ax1.set_xticks(minor_ticks_x, minor=True)
    ax1.set_yticks(major_ticks_y)
    ax1.set_yticks(minor_ticks_y, minor=True)

    # Or if you want different settings for the grids:
    ax1.grid(which='minor', alpha=0.2)
    ax1.grid(which='major', alpha=0.5)

    # trait distribution plot for BH
    ax2 = fig.add_subplot(gs[0, 2])
    # ax2.spines["top"].set_visible(False)
    # ax2.spines["right"].set_visible(False)
    # ax2.spines["left"].set_visible(False)
    # ax2.get_xaxis().set_ticks([])
    # ax2.get_yaxis().set_ticks([])
    ax2.set_title('DR model')
    # add x label and y label
    ax2.set_xlabel('Count')
    # ax1.set_ylabel('Trait values')
    # major_ticks_x2 = np.arange(0, 41, 10)
    #
    # ax2.set_xticks(major_ticks_x2)
    # ax2.set_xticklabels(major_ticks_x2)
    ax2.set_ylim(global_min,global_max)
    # ax2.axes.get_xaxis().set_ticklabels([])
    ax2.axes.get_yaxis().set_ticklabels([])
    ax2.hist(merge_RI_dr, bins = 100, color="#95d0fc", alpha=0.5, orientation='horizontal')

    # trait distribution plot for RI
    ax3 = fig.add_subplot(gs[0, 3])
    # ax3.spines["top"].set_visible(False)
    # ax3.spines["right"].set_visible(False)
    # ax3.spines["left"].set_visible(False)
    ax3.set_title('DK model')
    # add x label and y label
    ax3.set_xlabel('Count')
    # major_ticks_x3 = np.arange(0, 41, 10)
    #
    # ax3.set_xticks(major_ticks_x3)
    # ax3.set_xticklabels(major_ticks_x3)
    ax3.set_ylim(global_min,global_max)
    # ax3.axes.get_xaxis().set_ticklabels([])
    ax3.axes.get_yaxis().set_ticklabels([])
    ax3.hist(merge_RI_dk, bins = 100, color="#fc5a50",alpha=0.5, orientation='horizontal')

    plt.show()
    return fig


# trait 1 v.s. trait 2 dotplot
def dotplot(traitdata):
    # read simulated data
    stat_rate_trait_RI_dr = traitdata[0]
    stat_rate_trait_RI_dk = traitdata[1]
    stat_rate_popu_RI_dr = traitdata[2]
    stat_rate_popu_RI_dk = traitdata[3]

    # find out extinct species and remove the responding trait values
    ext_index_RI_dr = np.where(stat_rate_popu_RI_dr == 0)
    ext_index_RI_dk = np.where(stat_rate_popu_RI_dk == 0)
    statplot_trait_RI_dr = stat_rate_trait_RI_dr
    statplot_trait_RI_dr[ext_index_RI_dr[0],ext_index_RI_dr[1]] = np.nan
    # statplot_trait_RI_dr_sorted = np.sort(statplot_trait_RI_dr)
    statplot_trait_RI_dk = stat_rate_trait_RI_dk
    statplot_trait_RI_dk[ext_index_RI_dk[0],ext_index_RI_dk[1]] = np.nan
    # statplot_trait_RI_dk_sorted = np.sort(statplot_trait_RI_dk)
    # filter the missing data
    mask_RI_dr = ~np.isnan(statplot_trait_RI_dr)
    filtered_data_RI_dr = [d[m] for d, m in zip(statplot_trait_RI_dr.T, mask_RI_dr.T)]
    mask_RI_dk = ~np.isnan(statplot_trait_RI_dk)
    filtered_data_RI_dk = [d[m] for d, m in zip(statplot_trait_RI_dk.T, mask_RI_dk.T)]
    # convert the list to an array
    merge_RI_dr = np.concatenate(filtered_data_RI_dr, axis=0)
    merge_RI_dk = np.concatenate(filtered_data_RI_dk, axis=0)

    statplot_popu_RI_dr = stat_rate_popu_RI_dr
    statplot_popu_RI_dr[ext_index_RI_dr[0], ext_index_RI_dr[1]] = np.nan
    # statplot_popu_RI_dr_sorted = np.sort(statplot_popu_RI_dr)
    statplot_popu_RI_dk = stat_rate_popu_RI_dk
    statplot_popu_RI_dk[ext_index_RI_dk[0], ext_index_RI_dk[1]] = np.nan
    # statplot_popu_RI_sorted = np.sort(statplot_popu_RI_dk)

    # filter the missing data
    mask_popu_RI_dr = ~np.isnan(statplot_popu_RI_dr)
    filtered_popu_data_RI_dr = [d[m] for d, m in zip(statplot_popu_RI_dr.T, mask_popu_RI_dr.T)]
    mask_popu_RI_dk = ~np.isnan(statplot_popu_RI_dk)
    filtered_popu_data_RI_dk = [d[m] for d, m in zip(statplot_popu_RI_dk.T, mask_popu_RI_dk.T)]
    # convert the list to an array
    merge_popu_RI_dr = np.concatenate(filtered_popu_data_RI_dr, axis=0)
    merge_popu_RI_dk = np.concatenate(filtered_popu_data_RI_dk, axis=0)


    # create a fig
    fig = plt.figure(1, figsize=(12, 9))
    # determine the limits for axises
    min_RI_dr_t = np.amin(merge_RI_dr)-5
    max_RI_dr_t = np.amax(merge_RI_dr)+5
    min_RI_dk_t = np.amin(merge_RI_dk)-5
    max_RI_dk_t = np.amax(merge_RI_dk)+5
    global_min_t =  min(min_RI_dr_t, min_RI_dk_t)
    global_max_t =  max(max_RI_dr_t, max_RI_dk_t)

    min_RI_dr_p = np.amin(merge_popu_RI_dr) - 5
    max_RI_dr_p = np.amax(merge_popu_RI_dr) + 5
    min_RI_dk_p = np.amin(merge_popu_RI_dk) - 5
    max_RI_dk_p = np.amax(merge_popu_RI_dk) + 5
    global_min_p =  min(min_RI_dr_p,min_RI_dk_p)
    global_max_p =  max(max_RI_dr_p,max_RI_dk_p)
    # multiple plots arrangement
    gs = gridspec.GridSpec(1, 2)

    # Create an axes instance
    ax1 = fig.add_subplot(gs[0,0])
    # boxplot BH trait data
    bh = ax1.scatter(merge_popu_RI_dr, merge_RI_dr, s=10, c = "#95d0fc", marker='o', alpha = 0.5)

    # add legends
    # ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['BH', 'RI'], loc='upper left')
    # add title
    ax1.set_title('Trait v.s. Population under RI-dr')
    # add x label and y label
    ax1.set_xlabel('Population size')
    ax1.set_ylabel('Trait values')
    ax1.set_ylim(global_min_t,global_max_t)

    # Create an axes instance
    ax2 = fig.add_subplot(gs[0, 1])
    # boxplot BH trait data
    ri = ax2.scatter(merge_popu_RI_dk, merge_RI_dk, s=10, c="#fc5a50", marker='o', alpha=0.5)

    # add legends
    # ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['BH', 'RI'], loc='upper left')
    # add title
    ax2.set_title('Trait v.s. Population under RI-dk')
    # add x label and y label
    ax2.set_xlabel('Population size')
    ax2.set_ylabel('')
    ax2.set_ylim(global_min_t, global_max_t)

    plt.show()
    return fig




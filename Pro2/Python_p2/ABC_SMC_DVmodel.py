import sys, os
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
import numpy as np
from DV_model_sim_along_phy import DVtraitsim_tree
import scipy.stats
from scipy.stats import norm
import timeit

def single_trait_sim(par,file,replicate):
    do = 0
    while(do == 0):
        sim = DVtraitsim_tree(file = file, replicate = replicate,K = 10000000, gamma1 = par[0],a = par[1])
        if sim[2]:
            do = 1
        else:
            print('Retry')
            do = 0
    trait_RI_dr = sim[0]
    population_RI_dr = sim[1]
    traitVar = sim[3]
    evo_time, total_species = sim[0].shape
    evo_time = evo_time - 1
    # valid = sim[2]
    # empirical data for trait and population
    trait_dr_tips = trait_RI_dr[evo_time, :][~np.isnan(trait_RI_dr[evo_time, :])]
    population_tips = population_RI_dr[evo_time, :][~np.isnan(population_RI_dr[evo_time, :])]
    traitVar_tips = traitVar[evo_time, :][~np.isnan(traitVar[evo_time, :])]
    simtips = np.array([trait_dr_tips,population_tips,traitVar_tips])
    return simtips

def PosNormal(mean, sigma):
    x = np.random.normal(mean,sigma,1)
    return(x if x>=0 else PosNormal(mean,sigma))
#
# def calibration(samplesize, priorpar, treefile,calidata_file, calmode = 'uni'):
#     collection = np.zeros(shape=(samplesize,2))
#     cali_traitdata = ([])
#     cali_popdata = ([])
#     cali_vardata = ([])
#
#     if calmode == 'uni':
#         uniform_gamma = np.random.uniform(priorpar[0],priorpar[1],samplesize)
#         uniform_a = np.random.uniform(priorpar[2],priorpar[3],samplesize)
#         do = True
#     elif calmode == 'nor':
#         uniform_gamma = np.zeros(samplesize)
#         uniform_a = np.zeros(samplesize)
#         for i in range(samplesize):
#             uniform_gamma[i] = PosNormal(priorpar[0],priorpar[1])
#             uniform_a[i] = PosNormal(priorpar[2],priorpar[3])
#         do = True
#     elif calmode == 'self':
#         uniform_gamma = np.repeat(priorpar[0],samplesize)
#         uniform_a = np.repeat(priorpar[1],samplesize)
#         do = True
#     else:
#         print('Please indicate one mode!')
#         uniform_a = 0
#         uniform_gamma = 0
#         do = False
#
#     if do == True:
#         for i in range(samplesize):
#             print(i)
#             par_cal = np.zeros(2)
#             par_cal[0] = uniform_gamma[i]
#             par_cal[1] = uniform_a[i]
#             sample_cal =  single_trait_sim(par = par_cal,file = treefile,replicate=0)
#             collection[i] = par_cal  #,[diff],[diff_sort]
#             cali_traitdata.append(sample_cal[0])
#             cali_popdata.append(sample_cal[1])
#             cali_vardata.append(sample_cal[2])
#         cali_traitdataarray = np.asarray(cali_traitdata)
#         cali_popdataarray = np.asarray(cali_popdata)
#         cali_vardataarray = np.asarray(cali_vardata)
#         calipar=collection[:,:2]
#         np.savez(calidata_file,calipar = calipar, calitrait=cali_traitdataarray, calipop =cali_popdataarray, calivar=cali_vardataarray)
#         return collection


def calibration(samplesize, priorpar, treefile,calidata_file):
    collection = np.zeros(shape=(samplesize,2))
    cali_traitdata = ([])
    cali_popdata = ([])
    cali_vardata = ([])
    par_picked = ([])
    for i in range(samplesize):
        do = 0
        while(do==0):
            uniform_gamma = np.random.uniform(priorpar[0], priorpar[1], 1)
            uniform_a = np.random.uniform(priorpar[2], priorpar[3], 1)
            print(i)
            par_cal = np.zeros(2)
            par_cal[0] = uniform_gamma
            par_cal[1] = uniform_a
            par_picked.append(par_cal)
            sample_cal =  DVtraitsim_tree(file = treefile, replicate = 0,K = 100000, gamma1 = uniform_gamma,a = uniform_a)
            if sample_cal[2]:
                do = 1
            else:
                print('Retry')
                do = 0
        trait_RI_dr = sample_cal[0]
        population_RI_dr = sample_cal[1]
        traitVar = sample_cal[3]
        evo_time, total_species = sample_cal[0].shape
        evo_time = evo_time - 1
        trait_dr_tips = trait_RI_dr[evo_time, :][~np.isnan(trait_RI_dr[evo_time, :])]
        population_tips = population_RI_dr[evo_time, :][~np.isnan(population_RI_dr[evo_time, :])]
        traitVar_tips = traitVar[evo_time, :][~np.isnan(traitVar[evo_time, :])]
        collection[i] = np.array(par_cal)
        cali_traitdata.append(trait_dr_tips)
        cali_popdata.append(population_tips)
        cali_vardata.append(traitVar_tips)
    cali_traitdataarray = np.asarray(cali_traitdata)
    cali_popdataarray = np.asarray(cali_popdata)
    cali_vardataarray = np.asarray(cali_vardata)
    calipar=collection[:,:2]
    par_pickedarray = np.asarray(par_picked)
    np.savez(calidata_file+'picked',picked = par_pickedarray)
    np.savez(calidata_file,calipar = calipar, calitrait=cali_traitdataarray, calipop =cali_popdataarray, calivar=cali_vardataarray)
    return collection



def ABC_acceptance(par,delta,obs,sort, file,abcmode='mean'):
    sample = single_trait_sim(par = par,file = file,replicate=0)
    if sort == 0:
        samplearray = np.array([sample[0],sample[2]])
        obsarray = np.array([obs[0],obs[2]])
        if abcmode == 'mean':
            diff = np.linalg.norm(sample[0] - obs[0])
        elif abcmode == 'variance':
            diff = np.linalg.norm(sample[2] - obs[2])
        elif abcmode == 'both':
            diff = np.linalg.norm(samplearray - obsarray)
        else:
            print('Please indicate mode')
            diff = np.inf
        if diff < delta:
            return True
        else:
            return False
    else:
        samplearray = np.array([sample[0], sample[2]])
        obsarray = np.array([obs[0], obs[2]])
        samplearray_sort = samplearray[:,samplearray[0,:].argsort()]
        obsarray_sort = obsarray[:,obsarray[0,:].argsort()]
        if abcmode == 'mean':
            diff_sort = np.linalg.norm(np.sort(sample[0]) - np.sort(obs[0]))
        elif abcmode == 'variance':
            diff_sort = np.linalg.norm(np.sort(sample[2]) - np.sort(obs[2]))
        elif abcmode == 'both':
            diff_sort = np.linalg.norm(samplearray_sort - obsarray_sort)
        else:
            print('Please indicate mode')
            diff_sort = np.inf
        if diff_sort < delta:
            return True
        else:
            return False




def MCMC_ABC(startvalue, iterations,delta,obs,sort,priorpar, file, mcmcmode = 'uni',abcmode='mean'):
    tic = timeit.default_timer()
    MCMC = np.zeros(shape=(iterations+1,2))
    MCMC[0,] = startvalue
    par_jump = np.empty(2)
    if mcmcmode == 'uni':
        for i in range(iterations):
            par_jump[0] = np.random.uniform(priorpar[0],priorpar[1])
            par_jump[1] = np.random.uniform(priorpar[2],priorpar[3])

            if (ABC_acceptance(par_jump,delta = delta, obs = obs,sort = sort, file = file,abcmode=abcmode)):
                MCMC[i+1,] = par_jump
                print("MCMC : %d Accepted" % (i+1))

            else:
                MCMC[i + 1,] = MCMC[i ,]
                print("MCMC : %d Rejected" % (i+1))
    elif mcmcmode == 'nor':
        for i in range(iterations):
            par_jump[0] = abs(np.random.normal(loc=MCMC[i, 0], scale=0.01))
            par_jump[1] = abs(np.random.normal(loc=MCMC[i, 1], scale=0.01))

            pro = np.random.uniform(0,1,1)[0]
            pro_gamma1 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(par_jump[0])
            pro_gamma2 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(MCMC[i ,0])
            pro_a1 = scipy.stats.norm(priorpar[2], priorpar[3]).pdf(par_jump[1])
            pro_a2 = scipy.stats.norm(priorpar[2], priorpar[3]).pdf(MCMC[i ,1])

            pro_ratio = (pro_gamma1*pro_a1)/(pro_gamma2*pro_a2)
            accept_criterion = np.min([1,pro_ratio])
            if ABC_acceptance(par = par_jump, delta=delta, obs=obs, sort=sort, file = file,abcmode=abcmode
                              ) and (pro <= accept_criterion):
                MCMC[i + 1,] = par_jump
                print("MCMC : %d Accepted" % (i + 1))
            else:
                MCMC[i + 1,] = MCMC[i,]
                print("MCMC : %d Rejected" % (i + 1))
    toc = timeit.default_timer()
    elapse = toc - tic
    timetext = 'Elapsed time: %.2f' % elapse
    print(timetext)
    return MCMC


# SMC ABC
def SMC_ABC(timestep, particlesize, obs, epsilon, prior, file, sort = 0):
    tic = timeit.default_timer()
    d = np.zeros(shape = (timestep, particlesize))  #distance matrix of simulations and obs
    gamma = np.zeros(shape = (timestep, particlesize))  # gamma jumps
    a =  np.zeros(shape = (timestep, particlesize))     # a  jumps
    total_simulation = np.zeros(timestep)
    # prior information [mean_gamma,var_gamma,mean_a,var_a]
    gamma_prior_mean =  prior[0]
    gamma_prior_var = prior[1]
    a_prior_mean = prior[2]
    a_prior_var = prior[3]
    # Initial thredhold
    epsilon = epsilon
    # Weight vectors for gamma and a
    weight_gamma = np.zeros(shape = (timestep, particlesize))
    weight_gamma.fill(1/particlesize)
    weight_a = np.zeros(shape = (timestep, particlesize))
    weight_a.fill(1/particlesize)
    for t in range(timestep):
        sim_count = 0
        str = 'Time step : %d;' % t
        #Initial round
        if t == 0:
            for i in range(particlesize):
                str_p = str + ' Particle size : %d' % i
                print(str_p)
                d[t,i] = epsilon + 1
                while d[t,i] > epsilon:
                    sim_count += 1
                    # draw parameters from prior information
                    propose_gamma = abs(np.random.normal(gamma_prior_mean,gamma_prior_var))
                    propose_a = abs(np.random.normal(a_prior_mean,a_prior_var))
                    par = [propose_gamma,propose_a]
                    # simulate under the parameters
                    sample = single_trait_sim(par = par, file = file,replicate=0)
                    samplearray = np.array([sample[0], sample[2]])
                    obsarray = np.array([obs[0], obs[2]])
                    # calculate the distance between simulation and obs
                    if sort == 0:
                        diff = np.linalg.norm(samplearray - obsarray)
                    else:
                        samplearray_sort = samplearray[:, samplearray[0, :].argsort()]
                        obsarray_sort = obsarray[:, obsarray[0, :].argsort()]
                        diff = np.linalg.norm(samplearray_sort - obsarray_sort)
                    d[t,i] = diff
                # record the accepted values
                gamma[t, i] = propose_gamma
                a[t,i] = propose_a
        else:
            # shrink the threshold by 75% for each time step
            epsilon = np.append(epsilon, np.percentile(d[t-1,],75))
            # calculate weighted variance of the parameters at previous time step
            gamma_pre_mean = np.sum(gamma[t-1,] * weight_gamma[t-1,])
            gamma_pre_var = np.sum(( gamma[t-1,] - gamma_pre_mean)**2 * weight_gamma[t-1,])
            a_pre_mean = np.sum(a[t - 1,] * weight_a[t - 1,])
            a_pre_var = np.sum((a[t - 1,] - a_pre_mean) ** 2 * weight_a[t - 1,])
            for i in range(particlesize):
                str_p = str + ' Particle size : %d' % i
                print(str_p)
                d[t, i] = epsilon[t] + 1
                while d[t,i] > epsilon[t]:
                    sim_count += 1
                    # sample the parameters by the weight
                    sample_gamma_index = np.random.choice(particlesize,1, p = weight_gamma[t-1,])
                    sample_a_index = np.random.choice(particlesize,1, p = weight_a[t-1,])
                    # mean of the sample for gamma
                    propose_gamma0 = gamma[t-1,sample_gamma_index-1]
                    # draw new gamma with mean and variance
                    propose_gamma = abs(np.random.normal(propose_gamma0,np.sqrt(2*gamma_pre_var)))
                    # mean of the sample for a
                    propose_a0 = a[t-1,sample_a_index-1]
                    # draw new a with mean and variance
                    propose_a = abs(np.random.normal(propose_a0,np.sqrt(2* a_pre_var)))
                    par = [propose_gamma,propose_a]
                    sample = single_trait_sim(par, file = file,replicate=0)
                    samplearray = np.array([sample[0], sample[2]])
                    obsarray = np.array([obs[0], obs[2]])
                    # calculate the distance between simulation and obs
                    if sort == 0:
                        diff = np.linalg.norm(samplearray - obsarray)
                    else:
                        samplearray_sort = samplearray[:, samplearray[0, :].argsort()]
                        obsarray_sort = obsarray[:, obsarray[0, :].argsort()]
                        diff = np.linalg.norm(samplearray_sort - obsarray_sort)
                    d[t, i] = diff
                gamma[t, i] = propose_gamma
                a[t,i] = propose_a
                # compute new weights for gamma and a
                weight_gamma_denominator = np.sum(weight_gamma[t-1,]* norm.pdf(propose_gamma,gamma[t-1,] ,
                                                                               np.sqrt(2*gamma_pre_var)))
                weight_gamma_numerator = norm.pdf(propose_gamma,gamma_prior_mean,gamma_prior_var)
                weight_gamma[t,i] = weight_gamma_numerator / weight_gamma_denominator

                weight_a_denominator = np.sum(weight_a[t - 1,] * norm.pdf(propose_a, a[t - 1,],
                                                                                  np.sqrt(2 * a_pre_var)))
                weight_a_numerator = norm.pdf(propose_a, a_prior_mean, a_prior_var)
                weight_a[t, i] = weight_a_numerator / weight_a_denominator
        # normalize the weights
        total_simulation[t] = sim_count
        weight_gamma[t,] = weight_gamma[t,]/sum(weight_gamma[t,])
        weight_a[t,] = weight_a[t,]/sum(weight_a[t,])
    # create the list for output
    SMC_ABC = {'gamma': gamma, 'a': a, 'weight_gamma':weight_gamma,'weight_a':weight_a,'error':epsilon,'diff':d,
               'tot_sim':total_simulation}
    toc = timeit.default_timer()
    elapse = toc - tic
    timetext = 'Elapsed time: %.2f' % elapse
    print(timetext)
    return SMC_ABC



# SMC ABC for model selection
def SMC_ABC_MS(timestep, particlesize, obs, epsilon, prior, file, sort = 0):
    tic = timeit.default_timer()
    d = np.zeros(shape = (timestep, particlesize))  #distance matrix of simulations and obs
    model = np.zeros(shape = (timestep, particlesize))
    gamma = np.zeros(shape = (timestep, particlesize))  # gamma jumps
    a =  np.zeros(shape = (timestep, particlesize))     # a  jumps
    total_simulation = np.zeros(timestep)
    # prior information [mean_gamma,var_gamma,mean_a,var_a]
    gamma_prior_mean =  prior[0]
    gamma_prior_var = prior[1]
    a_prior_mean = prior[2]
    a_prior_var = prior[3]
    # Initialize thredhold
    epsilon = epsilon
    # Weight vectors for gamma and a
    weight_gamma = np.zeros(shape = (timestep, particlesize))
    weight_gamma.fill(1/particlesize)
    weight_a = np.zeros(shape = (timestep, particlesize))
    weight_a.fill(1/particlesize)
    for t in range(timestep):
        sim_count = 0
        str = 'Time step : %d;' % t
        #Initial round
        if t == 0:
            for i in range(particlesize):
                str_p = str + ' Particle size : %d' % i
                print(str_p)
                d[t,i] = epsilon + 1
                while d[t,i] > epsilon:
                    sim_count += 1
                    # Sample model and parameters
                    propose_model = np.random.randint(4)
                    propose_gamma = abs(np.random.normal(gamma_prior_mean, gamma_prior_var))
                    propose_a = abs(np.random.normal(a_prior_mean, a_prior_var))
                    if propose_model == 0:  # BM model
                        # draw parameters from prior information
                        propose_gamma = 0
                        propose_a = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 1:  # Competition model
                        # draw parameters from prior information
                        propose_gamma = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 2: # OU model / Natural selection model
                        # draw parameters from prior information
                        propose_a = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 3: # Natural selection & competition model
                        # draw parameters from prior information
                        par = [propose_gamma,propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par = par, file = file,replicate=0)

                    samplearray = np.array([sample[0], sample[2]])
                    obsarray = np.array([obs[0], obs[2]])
                    # calculate the distance between simulation and obs
                    if sort == 0:
                        diff = np.linalg.norm(samplearray - obsarray)
                    else:
                        samplearray_sort = samplearray[:, samplearray[0, :].argsort()]
                        obsarray_sort = obsarray[:, obsarray[0, :].argsort()]
                        diff = np.linalg.norm(samplearray_sort - obsarray_sort)
                    d[t, i] = diff
                # record the accepted values
                gamma[t, i] = propose_gamma
                a[t,i] = propose_a
                model[t,i] = propose_model
        else:
            # shrink the threshold by 75% for each time step
            epsilon = np.append(epsilon, np.percentile(d[t-1,],75))
            # calculate weighted variance of the parameters at previous time step
            gamma_pre_mean = np.sum(gamma[t-1,] * weight_gamma[t-1,])
            gamma_pre_var = np.sum(( gamma[t-1,] - gamma_pre_mean)**2 * weight_gamma[t-1,])
            a_pre_mean = np.sum(a[t - 1,] * weight_a[t - 1,])
            a_pre_var = np.sum((a[t - 1,] - a_pre_mean) ** 2 * weight_a[t - 1,])
            for i in range(particlesize):
                str_p = str + ' Particle size : %d' % i
                print(str_p)
                d[t, i] = epsilon[t] + 1
                while d[t,i] > epsilon[t]:
                    sim_count += 1
                    # Sample model
                    propose_model = np.random.randint(4)
                    # sample the parameters by the weight
                    sample_gamma_index = np.random.choice(particlesize,1, p = weight_gamma[t-1,])
                    sample_a_index = np.random.choice(particlesize,1, p = weight_a[t-1,])
                    # mean of the sample for gamma
                    propose_gamma0 = gamma[t-1,sample_gamma_index-1]
                    # draw new gamma with mean and variance
                    propose_gamma = abs(np.random.normal(propose_gamma0,np.sqrt(2*gamma_pre_var)))
                    # mean of the sample for a
                    propose_a0 = a[t-1,sample_a_index-1]
                    # draw new a with mean and variance
                    propose_a = abs(np.random.normal(propose_a0,np.sqrt(2* a_pre_var)))
                    if propose_model == 0:  # BM model
                        # draw parameters from prior information
                        propose_gamma = 0
                        propose_a = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 1:  # Competition model
                        # draw parameters from prior information
                        propose_gamma = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 2:  # OU model / Natural selection model
                        # draw parameters from prior information
                        propose_a = 0
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    elif propose_model == 3:  # Natural selection & competition model
                        # draw parameters from prior information
                        par = [propose_gamma, propose_a]
                        # simulate under the parameters
                        sample = single_trait_sim(par=par, file=file,replicate=0)
                    samplearray = np.array([sample[0], sample[2]])
                    obsarray = np.array([obs[0], obs[2]])
                    # calculate the distance between simulation and obs
                    if sort == 0:
                        diff = np.linalg.norm(samplearray - obsarray)
                    else:
                        samplearray_sort = samplearray[:, samplearray[0, :].argsort()]
                        obsarray_sort = obsarray[:, obsarray[0, :].argsort()]
                        diff = np.linalg.norm(samplearray_sort - obsarray_sort)
                    d[t, i] = diff
                gamma[t, i] = propose_gamma
                a[t,i] = propose_a
                model[t,i] = propose_model
                # compute new weights for gamma and a
                weight_gamma_denominator = np.sum(weight_gamma[t-1,]* norm.pdf(propose_gamma,gamma[t-1,] ,
                                                                               np.sqrt(2*gamma_pre_var)))
                weight_gamma_numerator = norm.pdf(propose_gamma,gamma_prior_mean,gamma_prior_var)
                weight_gamma[t,i] = weight_gamma_numerator / weight_gamma_denominator

                weight_a_denominator = np.sum(weight_a[t - 1,] * norm.pdf(propose_a, a[t - 1,],
                                                                                  np.sqrt(2 * a_pre_var)))
                weight_a_numerator = norm.pdf(propose_a, a_prior_mean, a_prior_var)
                weight_a[t, i] = weight_a_numerator / weight_a_denominator
        # normalize the weights
        total_simulation[t] = sim_count
        weight_gamma[t,] = weight_gamma[t,]/sum(weight_gamma[t,])
        weight_a[t,] = weight_a[t,]/sum(weight_a[t,])
    # create the list for output
    SMC_ABC_model = {'gamma': gamma, 'a': a, 'model': model, 'weight_gamma':weight_gamma,'weight_a':weight_a,'error':epsilon,'diff':d,
               'tot_sim':total_simulation}
    toc = timeit.default_timer()
    elapse = toc - tic
    timetext = 'Elapsed time: %.2f' % elapse
    print(timetext)
    return SMC_ABC_model
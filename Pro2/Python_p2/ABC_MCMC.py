import sys, os
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
import numpy as np
from Trait_sim_in_branches_stat import traitsim, drawplot, dotplot
import scipy.stats
from scipy.stats import norm

def single_trait_sim(par):
    sim = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par[0],gamma_K2=par[0],a = par[1],r = 1,theta = 0,K = 5000
                   , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate = 0)
    return sim

def PosNormal(mean, sigma):
    x = np.random.normal(mean,sigma,1)
    return(x if x>=0 else PosNormal(mean,sigma))

def calibrication(samplesize, priorpar, obs, mode = 'uni'):
    collection = np.zeros(shape=(samplesize,4))
    if mode == 'uni':
        uniform_gamma = np.random.uniform(priorpar[0],priorpar[1],samplesize)
        uniform_a = np.random.uniform(priorpar[2],priorpar[3],samplesize)
        do = True
    elif mode == 'nor':
        uniform_gamma = np.zeros(samplesize)
        uniform_a = np.zeros(samplesize)
        for i in range(samplesize):
            uniform_gamma[i] = PosNormal(priorpar[0],priorpar[1])
            uniform_a[i] = PosNormal(priorpar[2],priorpar[3])
        do = True
    elif mode == 'self':
        uniform_gamma = np.repeat(priorpar[0],samplesize)
        uniform_a = np.repeat(priorpar[1],samplesize)
        do = True
    else:
        print('Please indicate one mode!')
        uniform_a = 0
        uniform_gamma = 0
        do = False

    if do == True:
        for i in range(samplesize):
            print(i)
            par_cal = np.zeros(2)
            par_cal[0] = uniform_gamma[i]
            par_cal[1] = uniform_a[i]
            sample_cal =  single_trait_sim(par_cal)
            diff =  np.linalg.norm(sample_cal[0]-obs[0])
            diff_sort = np.linalg.norm(np.sort(sample_cal[0])-np.sort(obs[0]))
            collection[i] = np.concatenate((par_cal,[diff],[diff_sort]))
        return collection




def ABC_acceptance(par,delta,obs,sort):
    sample = single_trait_sim(par)
    if sort == 0:
        diff = np.linalg.norm(sample[0] - obs[0])
        if diff<delta:
            return True
        else:
            return False
    else:
        diff_sort = np.linalg.norm(np.sort(sample[0]) - np.sort(obs[0]))
        if diff_sort < delta:
            return True
        else:
            return False

def MCMC_ABC(startvalue, iterations,delta,obs,sort,priorpar,mode = 'uni'):
    MCMC = np.zeros(shape=(iterations+1,2))
    MCMC[0,] = startvalue
    par_jump = np.empty(2)
    if mode == 'uni':
        for i in range(iterations):
            par_jump[0] = abs(np.random.normal(loc=MCMC[i,0], scale= 0.01 ))
            par_jump[1] = abs(np.random.normal(loc=MCMC[i,1], scale= 0.01 ))

            if (ABC_acceptance(par_jump,delta = delta, obs = obs,sort = sort)):
                MCMC[i+1,] = par_jump
                print("MCMC : %d Accepted" % (i+1))

            else:
                MCMC[i + 1,] = MCMC[i ,]
                print("MCMC : %d Rejected" % (i+1))
    elif mode == 'nor':
        for i in range(iterations):
            par_jump[0] = abs(np.random.normal(loc=MCMC[i, 0], scale=0.01))
            par_jump[1] = abs(np.random.normal(loc=MCMC[i, 1], scale=0.01))

            pro = np.random.uniform(0,1,1)[0]
            pro_gamma1 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(par_jump[0])
            pro_gamma2 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(MCMC[i ,0])
            pro_a1 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(par_jump[1])
            pro_a2 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(MCMC[i ,1])

            pro_ratio = (pro_gamma1*pro_a1)/(pro_gamma2*pro_a2)
            accept_creterion = np.min([1,pro_ratio])
            if ABC_acceptance(par = par_jump, delta=delta, obs=obs, sort=sort) and (pro <= accept_creterion):
                MCMC[i + 1,] = par_jump
                print("MCMC : %d Accepted" % (i + 1))
            else:
                MCMC[i + 1,] = MCMC[i,]
                print("MCMC : %d Rejected" % (i + 1))

    return MCMC


# SMC ABC
def SMC_ABC(timestep, particlesize, obs, epsilon, prior,sort = 0):
    d = np.zeros(shape = (timestep, particlesize))  #distance matrix of simulations and obs
    gamma = np.zeros(shape = (timestep, particlesize))  # gamma jumps
    a =  np.zeros(shape = (timestep, particlesize))     # a  jumps
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
        str = 'Time step : %d;' % t
        #Initial round
        if t == 0:
            for i in range(particlesize):
                str_p = str + ' Particle size : %d' % i
                print(str_p)
                d[t,i] = epsilon + 1
                while d[t,i] > epsilon:
                    # draw parameters from prior information
                    propose_gamma = abs(np.random.normal(gamma_prior_mean,gamma_prior_var))
                    propose_a = abs(np.random.normal(a_prior_mean,a_prior_var))
                    par = [propose_gamma,propose_a]
                    # simulate under the parameters
                    sample = single_trait_sim(par)
                    # calculate the distance between simulation and obs
                    if sort == 0:
                        diff = np.linalg.norm(sample[0] - obs[0])
                    else:
                        diff = np.linalg.norm(np.sort(sample[0]) - np.sort(obs[0]))
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
                    sample = single_trait_sim(par)
                    if sort == 0:
                        diff = np.linalg.norm(sample[0] - obs[0])
                    else:
                        diff = np.linalg.norm(np.sort(sample[0]) - np.sort(obs[0]))
                    d[t,i] = diff
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
        weight_gamma[t,] = weight_gamma[t,]/sum(weight_gamma[t,])
        weight_a[t,] = weight_a[t,]/sum(weight_a[t,])
    # create the list for output
    SMC_ABC = {'gamma': gamma, 'a': a, 'weight_gamma':weight_gamma,'weight_a':weight_a,'error':epsilon,'diff':d}
    return SMC_ABC
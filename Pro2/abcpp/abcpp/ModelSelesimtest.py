import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp

sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from PhyloDiff_model_sim import Candimodels
from tp_update import tp_update
from bm_update import bm_update
from ou_update import ou_update
from dr_update import dr_update
from nh_update import nh_update

#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err

# prep data
# dir_path = os.path.dirname(os.path.realpath(__file__))
# files = dir_path + '/../tree_data/example1/'
dir_path = 'C:\\Liang\\Googlebox\\Research\\Project2\\'
files = dir_path + 'treesim_newexp/example1/'
savedir = dir_path+'modelsele/'
file2 = savedir + 'example1/modelsele.npy'

td = DVTreeData(path=files, scalar=1000)
print(td.total_species)
K=10e8
nu=1/(100*K)

sigma2 = 0.02  # Brownian Motion variance
meantrait = 0.0
del_mute = 'on'
allowmodeldie = 'off'
# the generating params for models

generating = 'TP'
# TP
gamma = 0.01
a = .1
# OU
gene_gamma_ou = 0.01
# DR
gene_gamma_dr = 0.01
gene_a_dr=.1
gene_m_dr=1
# NH
gene_gamma_nh = 0.01
gene_m_nh =1


modelnum = 5
population = 1000
total_population = population * modelnum
generations = 20
# let's try to find a true simulation:
# The generating paras
gene_obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                    initpop_sigma=10.0, break_on_mu=False)
gene_candiparam = np.array([0.0, 0.0, meantrait, 1.0, meantrait, 1.0])

# When the generating model is TP...
if generating == 'TP':
    for r in range(10000):
        obs = dvcpp.DVSim(td, gene_obs_param)
        if obs['sim_time'] == td.sim_evo_time:
            break
    if obs['sim_time'] < td.sim_evo_time:
        print('hopeless, does not compute.')
        sys.exit(-1)
    s = np.argsort(obs['Z'])
    obsN = obs['N'][s]
    obsZ = obs['Z'][s]
    obsV = obs['V'][s]
    obsN = np.nan_to_num(obsN)
    obsZ = np.nan_to_num(obsZ)
    obsV = np.nan_to_num(obsV)

# When the generating model is the candidates...
elif generating == 'BM':
    gene_candiparam[ 3] = 0.0

    gene_candiparam[ 5] = sigma2  # randomize delta
    gene_modelBM = Candimodels(td, gene_candiparam)
    obsZ = gene_modelBM['Z'][gene_modelBM['Z'].shape[0] - 1, :]

elif generating == 'OU':
    gene_candiparam[ 0] = gene_gamma_ou  # randomize 'gamma'
    gene_candiparam[ 3] = 0.0
    gene_candiparam[ 5] = sigma2  # randomize delta
    gene_modelOU = Candimodels(td, gene_candiparam)
    obsZ = gene_modelOU['Z'][gene_modelOU['Z'].shape[0] - 1, :]

elif generating == 'DR':
    gene_candiparam[ 0] = gene_gamma_dr  # randomize 'gamma'
    gene_candiparam[ 1] = gene_a_dr  # randomize 'a'
    gene_candiparam[ 3] = gene_m_dr  # randomize 'm'
    gene_candiparam[ 5] = sigma2  # randomize delta
    gene_modelDR = Candimodels(td, gene_candiparam)
    obsZ = gene_modelDR['Z'][gene_modelDR['Z'].shape[0] - 1, :]

elif generating == 'NH':
    gene_candiparam[ 0] = gene_gamma_nh  # randomize 'gamma'
    gene_candiparam[ 1] = 0
    gene_candiparam[ 3] = gene_m_nh # randomize 'm'
    gene_candiparam[ 5] = sigma2  # randomize delta
    gene_modelNH = Candimodels(td, gene_candiparam)
    obsZ = gene_modelNH['Z'][gene_modelNH['Z'].shape[0] - 1, :]

else:
    print('Please specify the generating model...')

print('Try to find the generating model: %s ...' % generating)
# let's try to find a true simulation:
obs_param = DVParam(gamma=0.01, a=0.5, K=K, nu=nu, r=1, theta=meantrait, Vmax=1, inittrait=meantrait, initpop=500000,
             initpop_sigma = 10.0, break_on_mu=False)
candiparam = np.array([0.0, 0.0, meantrait, 1.0, meantrait, 1.0])
# pop = dvcpp.DVSim(td, obs_param)


params_TP = np.tile(obs_param, (population, 1))  # duplicate
params_TP[:, 0] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'gamma'
params_TP[:, 1] = np.random.uniform(0.0, 1.0, params_TP.shape[0])  # randomize 'a'
params_TP[:, 3] = np.random.uniform(0.0, 1e-10, params_TP.shape[0])  # randomize 'nu'

params_BM = np.tile(candiparam,(population,1))
params_BM[:,3] = 0.0
params_BM[:, 5] = np.random.uniform(0.0, sigma2, params_BM.shape[0])  # randomize delta


params_OU = np.tile(candiparam,(population,1))
params_OU[:,0] = np.random.uniform(0.0, 1.0, params_OU.shape[0])   # randomize 'gamma'
params_OU[:,3]= 0.0
params_OU[:, 5] = np.random.uniform(0.0, sigma2, params_OU.shape[0])  # randomize delta


params_DR =  np.tile(candiparam,(population,1))
params_DR[:,0]= np.random.uniform(0.0, 1.0, params_DR.shape[0])         # randomize 'gamma'
params_DR[:,1]= np.random.uniform(0.0, 1.0, params_DR.shape[0])         # randomize 'a'
params_DR[:,3]= np.random.uniform(0.0, 5.0, params_DR.shape[0])         # randomize 'm'
params_DR[:, 5] = np.random.uniform(0.0, sigma2, params_DR.shape[0])  # randomize delta

params_nh =  np.tile(candiparam,(population,1))
params_nh[:,0]= np.random.uniform(0.0, 1.0, params_nh.shape[0])     # randomize 'gamma'
params_nh[:,1]= 0
params_nh[:,3]= np.random.uniform(0.0, 5.0, params_nh.shape[0])     # randomize 'm'
params_nh[:, 5] = np.random.uniform(0.0, sigma2, params_nh.shape[0])  # randomize delta

# model choice
model_index = np.array(range(modelnum))
model_params = np.repeat(model_index,repeats = population)
model_data = np.zeros(shape=(generations+1, total_population))
model_data[0,:] = model_params
propose_model = model_params

# store parameters used
gamma_data = np.zeros(shape=(generations, population))
a_data = np.zeros(shape=(generations, population))
nu_data = np.zeros(shape=(generations, population))

fitness= np.zeros(shape=(generations, total_population))
# Initialize the weights.
weight_model = np.zeros(total_population)
weight_model.fill(1 / total_population)

# weights for paras of TP
weight_gamma_TP = np.zeros(population)
weight_gamma_TP.fill(1 / population)
weight_a_TP = np.zeros(population)
weight_a_TP.fill(1 / population)
weight_nu_TP = np.zeros(population)
weight_nu_TP.fill(1 / population)

# weights for paras of BM
weight_del_BM = np.zeros(population)
weight_del_BM.fill(1 / population)


# weights for paras of OU
weight_gamma_OU = np.zeros(population)
weight_gamma_OU.fill(1 / population)
weight_del_OU = np.zeros(population)
weight_del_OU.fill(1 / population)


# weights for paras of dr
weight_gamma_dr = np.zeros(population)
weight_gamma_dr.fill(1 / population)
weight_a_dr = np.zeros(population)
weight_a_dr.fill(1 / population)
weight_m_dr = np.zeros(population)
weight_m_dr.fill(1 / population)
weight_del_dr = np.zeros(population)
weight_del_dr.fill(1 / population)

# weights for paras of nh
weight_gamma_nh = np.zeros(population)
weight_gamma_nh.fill(1 / population)
weight_m_nh = np.zeros(population)
weight_m_nh.fill(1 / population)
weight_del_nh = np.zeros(population)
weight_del_nh.fill(1 / population)


for g in range(generations):
    # model 0
    TP_sample_length = len(np.where(model_data[g,:]==0)[0])
    BM_sample_length = len(np.where(model_data[g,:]==1)[0])
    OU_sample_length = len(np.where(model_data[g,:]==2)[0])
    DR_sample_length = len(np.where(model_data[g,:]==3)[0])
    NH_sample_length = len(np.where(model_data[g,:]==4)[0])
    Z = np.zeros((1, td.total_species))

    if TP_sample_length>0:
        simmodelTP = dvcpp.DVSim(td, params_TP)
        valid_TP = np.where(simmodelTP['sim_time'] == td.sim_evo_time)[0]
        Z_modelTP = simmodelTP['Z'][valid_TP]
        i, j = argsort2D(Z_modelTP)
        Z_modelTP = Z_modelTP[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTP = np.nan_to_num(Z_modelTP)
        Z = np.vstack([Z,Z_modelTP])
            # V = np.nan_to_num(V)
            #GOF: Goodness of fit
    else:
        valid_TP = []

    if BM_sample_length>0:
        valid_BM = []
        print('BM start')
        iter_num_BM = 0
        # model 1
        for param_BM in params_BM:
            simmodelBM = Candimodels(td, param_BM)
            if simmodelBM['sim_time']==td.sim_evo_time:
                valid_BM.append(iter_num_BM)
                Z_modelBM = simmodelBM['Z'][simmodelBM['Z'].shape[0]-1,:]
                Z = np.vstack([Z, sorted(Z_modelBM)])
            iter_num_BM += 1
    else:
        valid_BM = []

    # model 2
    if OU_sample_length>0:
        valid_OU = []
        iter_num_OU = 0
        print('OU start')
        for param_OU in params_OU:
            simmodelOU = Candimodels(td, param_OU)
            if simmodelOU['sim_time']==td.sim_evo_time:
                valid_OU.append(iter_num_OU)
                Z_modelOU = simmodelOU['Z'][simmodelOU['Z'].shape[0]-1,:]
                Z = np.vstack([Z, sorted(Z_modelOU)])
            iter_num_OU += 1
    else:
        valid_OU = []

    if DR_sample_length>0:
        valid_DR = []
        iter_num_DR = 0
        print('DR start')

        # model 3
        for param_drury in params_DR:
            simmodeldr= Candimodels(td, param_drury)
            if simmodeldr['sim_time']==td.sim_evo_time:
                valid_DR.append(iter_num_DR)
                Z_modeldr = simmodeldr['Z'][simmodeldr['Z'].shape[0]-1,:]
                Z = np.vstack([Z, sorted(Z_modeldr)])
            iter_num_DR += 1
    else:
        valid_DR = []

    if NH_sample_length>0:
        Z_NH = np.zeros((1, td.total_species))
        valid_NH = []
        iter_num_NH = 0
        print('NH start')
        for param_nh in params_nh:
            simmodelnh= Candimodels(td, param_nh)
            if simmodelnh['sim_time']==td.sim_evo_time:
                valid_NH.append(iter_num_NH)
                Z_modelnh = simmodelnh['Z'][simmodelnh['Z'].shape[0]-1,:]
                Z = np.vstack([Z, sorted(Z_modelnh)])
            iter_num_NH += 1
    else:
        valid_NH = []

    Z = Z[1:,]

    valid = np.concatenate([valid_TP,np.array(valid_BM)+TP_sample_length,
                            np.array(valid_OU)+TP_sample_length+BM_sample_length,
                            np.array(valid_DR)+TP_sample_length+BM_sample_length+OU_sample_length,
                            np.array(valid_NH)+TP_sample_length+BM_sample_length+OU_sample_length+DR_sample_length
                            ]).astype(int)

    eudis = normalized_norm(Z, obsZ)
    fitness[g,valid] += 1.0 - eudis


    # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g,:])[-int(total_population// 4)]  # best 25%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]

    modelTPperc = len(np.where(propose_model[fit_index]==0)[0])/len(fit_index)
    modelBMperc = len(np.where(propose_model[fit_index]==1)[0])/len(fit_index)
    modelOUperc = len(np.where(propose_model[fit_index]==2)[0])/len(fit_index)
    modeldrperc = len(np.where(propose_model[fit_index]==3)[0])/len(fit_index)
    modelnhperc = len(np.where(propose_model[fit_index]==4)[0])/len(fit_index)

    print('Iteration = %d 25th Model TP: %.1f%% ; Model BM: %.1f%% ; Model OU: %.1f%% ; Model DR: %.1f%% ; Model NH: %.1f%%...'
          % (g,modelTPperc*100 , modelBMperc*100,modelOUperc*100 ,modeldrperc*100,modelnhperc*100))
    print('Average fitness: %f' % np.mean(fitness[g,fit_index]))
    # reevaluate the weight of the best fitted  models
    weight_model_bestfitted = weight_model[fit_index]*fitness[g,fit_index]/sum(weight_model[fit_index]*fitness[g,fit_index])
    # sample new models from the fitness of previous best fitted models
    sample_model_index = sorted(np.random.choice(fit_index, total_population, p=weight_model_bestfitted))

    if allowmodeldie == 'on':
        propose_model = propose_model[sample_model_index]
        previous_bestfitted_model = propose_model[fit_index]

    else:
        propose_model = model_params
        q5_TP = np.argsort(fitness[g, :population-1])[-int(population // 4)]  # best 25%
        q5_BM = np.argsort(fitness[g, population:2*population-1])[-int(population // 4)]+population  # best 25%
        q5_OU = np.argsort(fitness[g, 2*population:3*population-1])[-int(population // 4)]+2*population  # best 25%
        q5_DR = np.argsort(fitness[g,3*population:4*population-1])[-int(population // 4)]+3*population   # best 25%
        q5_NH = np.argsort(fitness[g, 4*population:])[-int(population // 4)]+4*population  # best 25%

        fit_index_TP = np.where(fitness[g, :population-1] > fitness[g, q5_TP])[0]
        fit_index_BM = np.where(fitness[g, population:2*population-1] > fitness[g, q5_BM])[0]+population
        fit_index_OU = np.where(fitness[g, 2*population:3*population-1] > fitness[g, q5_OU])[0]+2*population
        fit_index_DR = np.where(fitness[g, 3*population:4*population-1] > fitness[g, q5_DR])[0]+3*population
        fit_index_NH = np.where(fitness[g, 4*population:] > fitness[g, q5_NH])[0]+4*population
        previous_bestfitted_model = propose_model[np.concatenate([fit_index_TP,fit_index_BM,
                                                    fit_index_OU,fit_index_DR,
                                                    fit_index_NH])]
        paraTP_index = q5_TP
        paraBM_index = q5_BM-population
        paraOU_index = q5_OU-2*population
        paraDR_index = q5_DR-3*population
        paraNH_index = q5_NH-4*population
        chosengamma_TP,chosena_TP,chosennu_TP = np.mean(params_TP[paraTP_index,0]),np.mean(params_TP[paraTP_index,1]),\
                                                    np.mean(params_TP[paraTP_index, 3])

        chosengamma_OU = np.mean(params_OU[paraOU_index,0])

        chosengamma_DR,chosena_DR,chosenm_DR = np.mean(params_DR[paraDR_index,0]),np.mean(params_DR[paraDR_index,1]),\
                                                    np.mean(params_DR[paraDR_index, 3])

        chosengamma_NH,chosenm_NH = np.mean(params_nh[paraNH_index,0]),\
                                                    np.mean(params_nh[paraNH_index, 3])

        print('Mean estimates: TP gamma: %f ; a: %f ; nu: %f' % ( chosengamma_TP,chosena_TP,chosennu_TP))
        print('Mean estimates: BM gamma: %f ; a: %f ; del: %f' % ( 0.0,0.0,sigma2))
        print('Mean estimates: OU gamma: %f ; a: %f ; del: %f' % ( chosengamma_OU,0.0,sigma2))
        print('Mean estimates: DR gamma: %f ; a: %f ; m: %f' % ( chosengamma_DR,chosena_DR,chosenm_DR))
        print('Mean estimates: TP gamma: %f ; a: %f ; m: %f' % ( chosengamma_NH, 0.0,chosenm_NH))



    model_data[g+1,:] = propose_model

    if len(np.where(propose_model==0)[0])>0:

        # update TP paras and weights
        weight_gamma_TP,weight_a_TP,weight_nu_TP,propose_gamma_TP,propose_a_TP,propose_nu_TP=\
            tp_update(previous_bestfitted_model, propose_model, params_TP, weight_gamma_TP, weight_a_TP, weight_nu_TP)
        modelTP = np.where(propose_model==0)
        params_TP = np.tile(obs_param, (len(modelTP[0]), 1))
        params_TP[:, 0] = propose_gamma_TP
        params_TP[:, 1] = propose_a_TP
        params_TP[:, 3] = propose_nu_TP

    if len(np.where(propose_model==1)[0])>0:

        # update BM paras and weights
        weight_del_BM, propose_del_BM = \
            bm_update(previous_bestfitted_model, propose_model, params_BM, weight_del_BM)
        modelBM = np.where(propose_model == 1)
        params_BM = np.tile(candiparam, (len(modelBM[0]), 1))
        params_BM[:,3] = 0.0
        if del_mute=='on':
            params_BM[:, 5]=sigma2
            weight_del_BM.fill(1/len(weight_del_BM))
        else:
            params_BM[:, 5] = propose_del_BM

    if len(np.where(propose_model==2)[0])>0:

        # update OU paras and weights
        weight_gamma_OU, weight_del_OU, propose_gamma_OU,propose_del_OU = \
            ou_update(previous_bestfitted_model, propose_model, params_OU, weight_gamma_OU,weight_del_OU)
        modelOU = np.where(propose_model == 2)
        params_OU = np.tile(candiparam, (len(modelOU[0]), 1))
        params_OU[:, 0] = propose_gamma_OU
        params_OU[:, 3] = 0.0

        if del_mute=='on':
            params_OU[:, 5]=sigma2
            weight_del_OU.fill(1/len(weight_del_OU))

        else:
            params_OU[:, 5] = propose_del_OU
    #
    if len(np.where(propose_model==3)[0])>0:

        # update DR paras and weights
        weight_gamma_dr, weight_a_dr, weight_m_dr,weight_del_dr, propose_gamma_dr, propose_a_dr, propose_m_dr,propose_del_dr = \
            dr_update(previous_bestfitted_model, propose_model, params_DR, weight_gamma_dr, weight_a_dr,weight_m_dr, weight_del_dr)
        modelDR = np.where(propose_model == 3)
        params_DR = np.tile(candiparam, (len(modelDR[0]), 1))
        params_DR[:, 0] = propose_gamma_dr
        params_DR[:, 1] = propose_a_dr
        params_DR[:, 3] = propose_m_dr
        if del_mute=='on':
            params_DR[:, 5]=sigma2
            weight_m_dr.fill(1/len(weight_m_dr))

        else:
            params_DR[:, 5] = propose_del_dr
#
    if len(np.where(propose_model==4)[0])>0:

        # update nh paras and weights
        weight_gamma_nh,  weight_m_nh,weight_del_nh, propose_gamma_nh, propose_m_nh,propose_del_nh = \
            nh_update(previous_bestfitted_model, propose_model, params_nh, weight_gamma_nh,weight_m_nh, weight_del_nh)
        modelnh = np.where(propose_model == 4)
        params_nh = np.tile(candiparam, (len(modelnh[0]), 1))
        params_nh[:, 0] = propose_gamma_nh
        params_nh[:, 3] = propose_m_nh
        if del_mute=='on':
            params_nh[:, 5]=sigma2
            weight_del_nh.fill(1/len(weight_del_nh))

        else:
            params_nh[:, 5] = propose_del_nh
#

para_data = {'model_data': model_data,'fitness': fitness}
# file='C:/Liang/Code/Pro2/abcpp/abcpp/smcndata/'
# # file = '/home/p274981/abcpp/abcpp/'
# file2 = file + 'paraT11.npy'
# np.save(file2,para_data)
#
# smc = np.load(file2).item()
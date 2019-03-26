import numpy as np
from scipy.stats import norm


def bm_update(previous_bestfitted_model, propose_model, params_BM, weight_del_BM):
    # Room in BM model to sample new paras
    previous_bestfitted_index_BM = np.where(previous_bestfitted_model == 1)[0]-\
                                   len(np.where(previous_bestfitted_model == 0)[0])
    previous_del_BM = params_BM[previous_bestfitted_index_BM, 5]

    weight_del_BM = weight_del_BM[previous_bestfitted_index_BM] / sum(weight_del_BM[previous_bestfitted_index_BM])

    del_pre_mean_BM = np.sum(previous_del_BM * weight_del_BM)
    del_pre_var_BM = abs(np.max((np.sum((previous_del_BM - del_pre_mean_BM) ** 2 * weight_del_BM),
                             0.0001*np.max(previous_del_BM))))

    # sample parameters by the weights computed in last loop.
    population_BM = len(np.where(propose_model == 1)[0])
    sample_del_index_BM = np.random.choice(previous_bestfitted_index_BM, population_BM, p=weight_del_BM)

    # mean of the sample for gamma
    propose_del0_BM = params_BM[sample_del_index_BM, 5]
    # draw new gamma with mean and variance
    propose_del_BM = abs(np.random.normal(propose_del0_BM, np.sqrt(2 * del_pre_var_BM)))

    extend_weight_del_BM = weight_del_BM[sample_del_index_BM]

    # compute new weights for gamma and a
    weight_del_denominator_BM = np.sum(extend_weight_del_BM * norm.pdf(propose_del_BM, propose_del0_BM,
                                                                           np.sqrt(2 * del_pre_var_BM)))
    weight_del_numerator_BM = norm.pdf(propose_del_BM, del_pre_mean_BM, np.sqrt(2 * del_pre_var_BM))
    weight_del_BM = weight_del_numerator_BM / weight_del_denominator_BM

    weight_del_BM = weight_del_BM / sum(weight_del_BM)


    return weight_del_BM, propose_del_BM
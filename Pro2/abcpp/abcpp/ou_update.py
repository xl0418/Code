import numpy as np
from scipy.stats import norm


def ou_update(previous_bestfitted_model, propose_model, params_OU, weight_gamma_OU,weight_del_OU):
    # Room in TP model to sample new paras
    previous_bestfitted_index_OU = np.where(previous_bestfitted_model == 2)[0]- \
                                   len(np.where(previous_bestfitted_model == 0)[0])-\
                                    len(np.where(previous_bestfitted_model == 1)[0])
    previous_gamma_OU = params_OU[previous_bestfitted_index_OU, 0]
    previous_del_OU = params_OU[previous_bestfitted_index_OU, 5]

    weight_gamma_OU = weight_gamma_OU[previous_bestfitted_index_OU] / sum(weight_gamma_OU[previous_bestfitted_index_OU])
    weight_del_OU = weight_del_OU[previous_bestfitted_index_OU] / sum(weight_del_OU[previous_bestfitted_index_OU])

    gamma_pre_mean_OU = np.sum(previous_gamma_OU * weight_gamma_OU)
    gamma_pre_var_OU = np.sum((previous_gamma_OU - gamma_pre_mean_OU) ** 2 * weight_gamma_OU)
    del_pre_mean_OU = np.sum(previous_del_OU * weight_del_OU)
    del_pre_var_OU = np.sum((previous_del_OU - del_pre_mean_OU) ** 2 * weight_del_OU)

    # sample parameters by the weights computed in last loop.
    population_OU = len(np.where(propose_model == 2)[0])
    sample_gamma_index_OU = np.random.choice(previous_bestfitted_index_OU, population_OU, p=weight_gamma_OU)
    sample_del_index_OU = np.random.choice(previous_bestfitted_index_OU, population_OU, p=weight_del_OU)

    # mean of the sample for gamma
    propose_gamma0_OU = params_OU[sample_gamma_index_OU, 0]
    # draw new gamma with mean and variance
    propose_gamma_OU = abs(np.random.normal(propose_gamma0_OU, np.sqrt(2 * gamma_pre_var_OU)))
    # mean of the sample for a
    propose_del0_OU = params_OU[sample_del_index_OU, 5]
    # draw new a with mean and variance
    propose_del_OU = abs(np.random.normal(propose_del0_OU, np.sqrt(2 * del_pre_var_OU)))


    extend_weight_gamma_OU = weight_gamma_OU[sample_gamma_index_OU]
    extend_weight_del_OU = weight_del_OU[sample_del_index_OU]

    # compute new weights for gamma and a
    weight_gamma_denominator_OU = np.sum(extend_weight_gamma_OU * norm.pdf(propose_gamma_OU, propose_gamma0_OU,
                                                                           np.sqrt(2 * gamma_pre_var_OU)))
    weight_gamma_numerator_OU = norm.pdf(propose_gamma_OU, gamma_pre_mean_OU, np.sqrt(2 * gamma_pre_var_OU))
    weight_gamma_OU = weight_gamma_numerator_OU / weight_gamma_denominator_OU

    weight_del_denominator_OU = np.sum(extend_weight_del_OU * norm.pdf(propose_del_OU, propose_del0_OU,
                                                                       np.sqrt(2 * del_pre_var_OU)))
    weight_del_numerator_OU = norm.pdf(propose_del_OU, del_pre_mean_OU, np.sqrt(2 * del_pre_var_OU))
    weight_del_OU = weight_del_numerator_OU / weight_del_denominator_OU
    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma_OU = weight_gamma_OU / sum(weight_gamma_OU)
    weight_del_OU = weight_del_OU / sum(weight_del_OU)

    return weight_gamma_OU, weight_del_OU, propose_gamma_OU,propose_del_OU
import numpy as np
from scipy.stats import norm
def tp_update(previous_bestfitted_index_TP,propose_model,params_TP,weight_gamma_TP,weight_a_TP,
              weight_nu_TP):
    # Room in TP model to sample new paras
    previous_gamma_TP = params_TP[previous_bestfitted_index_TP, 0]
    previous_a_TP = params_TP[previous_bestfitted_index_TP, 1]
    previous_nu_TP = params_TP[previous_bestfitted_index_TP, 4]

    weight_gamma_TP = weight_gamma_TP[previous_bestfitted_index_TP] / sum(weight_gamma_TP[previous_bestfitted_index_TP])
    weight_a_TP = weight_a_TP[previous_bestfitted_index_TP] / sum(weight_a_TP[previous_bestfitted_index_TP])
    weight_nu_TP = weight_nu_TP[previous_bestfitted_index_TP] / sum(weight_nu_TP[previous_bestfitted_index_TP])

    gamma_pre_mean_TP = np.sum(previous_gamma_TP * weight_gamma_TP)
    gamma_pre_var_TP = np.sum((previous_gamma_TP - gamma_pre_mean_TP) ** 2 * weight_gamma_TP)

    a_pre_mean_TP = np.sum(previous_a_TP * weight_a_TP)
    a_pre_var_TP = np.sum((previous_a_TP - a_pre_mean_TP) ** 2 * weight_a_TP)
    nu_pre_mean_TP = np.sum(previous_nu_TP * weight_nu_TP)
    nu_pre_var_TP = np.sum((previous_nu_TP - nu_pre_mean_TP) ** 2 * weight_nu_TP)

    # sample parameters by the weights computed in last loop.
    population_TP = len(np.where(propose_model == 0)[0])
    sample_gamma_index_TP = np.random.choice(previous_bestfitted_index_TP, population_TP, p=weight_gamma_TP)
    sample_a_index_TP = np.random.choice(previous_bestfitted_index_TP, population_TP, p=weight_a_TP)
    sample_nu_index_TP = np.random.choice(previous_bestfitted_index_TP, population_TP, p=weight_nu_TP)

    # mean of the sample for gamma
    propose_gamma0_TP = params_TP[sample_gamma_index_TP, 0]
    # draw new gamma with mean and variance
    propose_gamma_TP = abs(np.random.normal(propose_gamma0_TP, np.sqrt(2 * gamma_pre_var_TP)))
    # mean of the sample for a
    propose_a0_TP = params_TP[sample_a_index_TP, 1]
    # draw new a with mean and variance
    propose_a_TP = abs(np.random.normal(propose_a0_TP, np.sqrt(2 * a_pre_var_TP)))
    # mean of the sample for nu
    propose_nu0_TP = params_TP[sample_nu_index_TP, 4]
    # draw new nu with mean and variance
    propose_nu_TP = abs(np.random.normal(propose_nu0_TP, np.sqrt(2 * nu_pre_var_TP)))

    extend_weight_gamma_TP = weight_gamma_TP[previous_bestfitted_index_TP.searchsorted(sample_gamma_index_TP)]
    extend_weight_a_TP = weight_a_TP[previous_bestfitted_index_TP.searchsorted(sample_a_index_TP)]
    extend_weight_nu_TP = weight_nu_TP[previous_bestfitted_index_TP.searchsorted(sample_nu_index_TP)]

    # compute new weights for gamma and a
    weight_gamma_denominator_TP = np.sum(extend_weight_gamma_TP * norm.pdf(propose_gamma_TP, propose_gamma0_TP,
                                                                           np.sqrt(2 * gamma_pre_var_TP)))
    weight_gamma_numerator_TP = norm.pdf(propose_gamma_TP, gamma_pre_mean_TP, np.sqrt(2 * gamma_pre_var_TP))
    weight_gamma_TP = weight_gamma_numerator_TP / weight_gamma_denominator_TP

    weight_a_denominator_TP = np.sum(extend_weight_a_TP * norm.pdf(propose_a_TP, propose_a0_TP,
                                                                   np.sqrt(2 * a_pre_var_TP)))
    weight_a_numerator_TP = norm.pdf(propose_a_TP, a_pre_mean_TP, np.sqrt(2 * a_pre_var_TP))
    weight_a_TP = weight_a_numerator_TP / weight_a_denominator_TP

    weight_nu_denominator_TP = np.sum(extend_weight_nu_TP * norm.pdf(propose_nu_TP, propose_nu0_TP,
                                                                     np.sqrt(2 * nu_pre_var_TP)))
    weight_nu_numerator_TP = norm.pdf(propose_nu_TP, nu_pre_mean_TP, np.sqrt(2 *nu_pre_var_TP))
    weight_nu_TP = weight_nu_numerator_TP / weight_nu_denominator_TP
    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma_TP = weight_gamma_TP / sum(weight_gamma_TP)
    weight_a_TP = weight_a_TP / sum(weight_a_TP)
    weight_nu_TP = weight_nu_TP / sum(weight_nu_TP)

    return weight_gamma_TP,weight_a_TP,weight_nu_TP,propose_gamma_TP,propose_a_TP,propose_nu_TP
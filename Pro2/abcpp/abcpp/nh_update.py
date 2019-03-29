import numpy as np
from scipy.stats import norm


def nh_update(previous_bestfitted_index_NH, propose_model, params_nh, weight_gamma_nh,weight_m_nh,
              weight_del_nh):
    # Room in TP model to sample new paras
    previous_bestfitted_index_nh = previous_bestfitted_index_NH
    previous_gamma_nh = params_nh[previous_bestfitted_index_nh, 0]
    previous_m_nh = params_nh[previous_bestfitted_index_nh, 3]
    previous_del_nh = params_nh[previous_bestfitted_index_nh, 5]

    weight_gamma_nh = weight_gamma_nh[previous_bestfitted_index_nh] / sum(weight_gamma_nh[previous_bestfitted_index_nh])
    weight_m_nh = weight_m_nh[previous_bestfitted_index_nh] / sum(weight_m_nh[previous_bestfitted_index_nh])
    weight_del_nh = weight_del_nh[previous_bestfitted_index_nh] / sum(weight_del_nh[previous_bestfitted_index_nh])

    gamma_pre_mean_nh = np.sum(previous_gamma_nh * weight_gamma_nh)
    gamma_pre_var_nh = abs(np.max((np.sum((previous_gamma_nh - gamma_pre_mean_nh) ** 2 * weight_gamma_nh),
                               0.01*np.max(previous_gamma_nh))))
    m_pre_mean_nh = np.sum(previous_m_nh * weight_m_nh)
    m_pre_var_nh = abs(np.max((np.sum((previous_m_nh - m_pre_mean_nh) ** 2 * weight_m_nh),
                           0.01*np.max(previous_m_nh))))
    del_pre_mean_nh = np.sum(previous_del_nh * weight_del_nh)
    del_pre_var_nh = abs(np.max((np.sum((previous_del_nh - del_pre_mean_nh) ** 2 * weight_del_nh),
                             0.0001*np.max(previous_del_nh))))

    # sample parameters by the weights computed in last loop.
    population_nh = len(np.where(propose_model == 2)[0])
    sample_gamma_index_nh = np.random.choice(previous_bestfitted_index_nh, population_nh, p=weight_gamma_nh)
    sample_m_index_nh = np.random.choice(previous_bestfitted_index_nh, population_nh, p=weight_m_nh)
    sample_del_index_nh = np.random.choice(previous_bestfitted_index_nh, population_nh, p=weight_del_nh)

    # mean of the sample for gamma
    propose_gamma0_nh = params_nh[sample_gamma_index_nh, 0]
    # draw new gamma with mean and variance
    propose_gamma_nh = abs(np.random.normal(propose_gamma0_nh, np.sqrt(2 * gamma_pre_var_nh)))
    # mean of the sample for nu
    propose_m0_nh = params_nh[sample_m_index_nh, 3]
    # draw new nu with mean and variance
    propose_m_nh = abs(np.random.normal(propose_m0_nh, np.sqrt(2 * m_pre_var_nh)))
    # mean of the sample for del
    propose_del0_nh = params_nh[sample_del_index_nh, 5]
    # draw new nu with mean and variance
    propose_del_nh = abs(np.random.normal(propose_del0_nh, np.sqrt(2 * del_pre_var_nh)))


    extend_weight_gamma_nh = weight_gamma_nh[previous_bestfitted_index_NH.searchsorted(sample_gamma_index_nh)]
    extend_weight_m_nh = weight_m_nh[previous_bestfitted_index_NH.searchsorted(sample_m_index_nh)]
    extend_weight_del_nh = weight_del_nh[previous_bestfitted_index_NH.searchsorted(sample_del_index_nh)]


    # compute new weights for gamma and a
    weight_gamma_denominator_nh = np.sum(extend_weight_gamma_nh * norm.pdf(propose_gamma_nh, propose_gamma0_nh,
                                                                           np.sqrt(2 * gamma_pre_var_nh)))
    weight_gamma_numerator_nh = norm.pdf(propose_gamma_nh, gamma_pre_mean_nh, np.sqrt(2 * gamma_pre_var_nh))
    weight_gamma_nh = weight_gamma_numerator_nh / weight_gamma_denominator_nh

    weight_m_denominator_nh = np.sum(extend_weight_m_nh * norm.pdf(propose_m_nh, propose_m0_nh,
                                                                     np.sqrt(2 * m_pre_var_nh)))
    weight_m_numerator_nh = norm.pdf(propose_m_nh, m_pre_mean_nh, np.sqrt(2 * m_pre_var_nh))
    weight_m_nh = weight_m_numerator_nh / weight_m_denominator_nh
    weight_del_denominator_nh = np.sum(extend_weight_del_nh * norm.pdf(propose_del_nh, propose_del0_nh,
                                                                     np.sqrt(2 * del_pre_var_nh)))
    weight_del_numerator_nh = norm.pdf(propose_del_nh, del_pre_mean_nh, np.sqrt(2 * del_pre_var_nh))
    weight_del_nh = weight_del_numerator_nh / weight_del_denominator_nh

    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma_nh = weight_gamma_nh / sum(weight_gamma_nh)
    weight_m_nh= weight_m_nh / sum(weight_m_nh)
    weight_del_nh= weight_del_nh / sum(weight_del_nh)


    return weight_gamma_nh, weight_m_nh,weight_del_nh, propose_gamma_nh, propose_m_nh,propose_del_nh
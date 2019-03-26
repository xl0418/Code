import numpy as np
from scipy.stats import norm


def dr_update(previous_bestfitted_model, propose_model, params_DR, weight_gamma_dr, weight_a_dr,weight_m_dr, weight_del_dr):
    # Room in TP model to sample new paras
    previous_bestfitted_index_dr = np.where(previous_bestfitted_model == 3)[0]-\
                            len(np.where(previous_bestfitted_model == 0)[0]) - \
                            len(np.where(previous_bestfitted_model == 1)[0]) - \
                                   len(np.where(previous_bestfitted_model == 2)[0])
    previous_gamma_dr = params_DR[previous_bestfitted_index_dr, 0]
    previous_a_dr = params_DR[previous_bestfitted_index_dr, 1]
    previous_m_dr = params_DR[previous_bestfitted_index_dr, 3]
    previous_del_dr = params_DR[previous_bestfitted_index_dr, 5]

    weight_gamma_dr = weight_gamma_dr[previous_bestfitted_index_dr] / sum(weight_gamma_dr[previous_bestfitted_index_dr])
    weight_a_dr = weight_a_dr[previous_bestfitted_index_dr] / sum(weight_a_dr[previous_bestfitted_index_dr])
    weight_m_dr = weight_m_dr[previous_bestfitted_index_dr] / sum(weight_m_dr[previous_bestfitted_index_dr])
    weight_del_dr = weight_del_dr[previous_bestfitted_index_dr] / sum(weight_del_dr[previous_bestfitted_index_dr])

    gamma_pre_mean_dr = np.sum(previous_gamma_dr * weight_gamma_dr)
    gamma_pre_var_dr = abs(np.max((np.sum((previous_gamma_dr - gamma_pre_mean_dr) ** 2 * weight_gamma_dr),
                               0.01*np.max(previous_gamma_dr))))
    a_pre_mean_dr = np.sum(previous_a_dr * weight_a_dr)
    a_pre_var_dr = abs(np.max((np.sum((previous_a_dr - a_pre_mean_dr) ** 2 * weight_a_dr),
                           0.01*np.max(previous_a_dr))))
    m_pre_mean_dr = np.sum(previous_m_dr * weight_m_dr)
    m_pre_var_dr = abs(np.max((np.sum((previous_m_dr - m_pre_mean_dr) ** 2 * weight_m_dr),
                           0.01*np.max(previous_m_dr))))
    del_pre_mean_dr = np.sum(previous_del_dr * weight_del_dr)
    del_pre_var_dr = abs(np.max((np.sum((previous_del_dr - del_pre_mean_dr) ** 2 * weight_del_dr),
                             0.0001*np.max(previous_del_dr))))

    # sample parameters by the weights computed in last loop.
    population_dr = len(np.where(propose_model == 3)[0])
    sample_gamma_index_dr = np.random.choice(previous_bestfitted_index_dr, population_dr, p=weight_gamma_dr)
    sample_a_index_dr = np.random.choice(previous_bestfitted_index_dr, population_dr, p=weight_a_dr)
    sample_m_index_dr = np.random.choice(previous_bestfitted_index_dr, population_dr, p=weight_m_dr)
    sample_del_index_dr = np.random.choice(previous_bestfitted_index_dr, population_dr, p=weight_del_dr)

    # mean of the sample for gamma
    propose_gamma0_dr = params_DR[sample_gamma_index_dr, 0]
    # draw new gamma with mean and variance
    propose_gamma_dr = abs(np.random.normal(propose_gamma0_dr, np.sqrt(2 * gamma_pre_var_dr)))
    # mean of the sample for a
    propose_a0_dr = params_DR[sample_a_index_dr, 1]
    # draw new a with mean and variance
    propose_a_dr = abs(np.random.normal(propose_a0_dr, np.sqrt(2 * a_pre_var_dr)))
    # mean of the sample for nu
    propose_m0_dr = params_DR[sample_m_index_dr, 3]
    # draw new nu with mean and variance
    propose_m_dr = abs(np.random.normal(propose_m0_dr, np.sqrt(2 * m_pre_var_dr)))
    # mean of the sample for del
    propose_del0_dr = params_DR[sample_del_index_dr, 3]
    # draw new nu with mean and variance
    propose_del_dr = abs(np.random.normal(propose_del0_dr, np.sqrt(2 * del_pre_var_dr)))


    extend_weight_gamma_dr = weight_gamma_dr[sample_gamma_index_dr]
    extend_weight_a_dr = weight_a_dr[sample_a_index_dr]
    extend_weight_m_dr = weight_m_dr[sample_m_index_dr]
    extend_weight_del_dr = weight_del_dr[sample_del_index_dr]


    # compute new weights for gamma and a
    weight_gamma_denominator_dr = np.sum(extend_weight_gamma_dr * norm.pdf(propose_gamma_dr, propose_gamma0_dr,
                                                                           np.sqrt(2 * gamma_pre_var_dr)))
    weight_gamma_numerator_dr = norm.pdf(propose_gamma_dr, gamma_pre_mean_dr, np.sqrt(2 * gamma_pre_var_dr))
    weight_gamma_dr = weight_gamma_numerator_dr / weight_gamma_denominator_dr

    weight_a_denominator_dr = np.sum(extend_weight_a_dr * norm.pdf(propose_a_dr, propose_a0_dr,
                                                                   np.sqrt(2 * a_pre_var_dr)))
    weight_a_numerator_dr = norm.pdf(propose_a_dr, a_pre_mean_dr, np.sqrt(2 * a_pre_var_dr))
    weight_a_dr = weight_a_numerator_dr / weight_a_denominator_dr

    weight_m_denominator_dr = np.sum(extend_weight_m_dr * norm.pdf(propose_m_dr, propose_m0_dr,
                                                                     np.sqrt(2 * m_pre_var_dr)))
    weight_m_numerator_dr = norm.pdf(propose_m_dr, m_pre_mean_dr, np.sqrt(2 * m_pre_var_dr))
    weight_m_dr = weight_m_numerator_dr / weight_m_denominator_dr
    weight_del_denominator_dr = np.sum(extend_weight_del_dr * norm.pdf(propose_del_dr, propose_del0_dr,
                                                                     np.sqrt(2 * del_pre_var_dr)))
    weight_del_numerator_dr = norm.pdf(propose_del_dr, del_pre_mean_dr, np.sqrt(2 * del_pre_var_dr))
    weight_del_dr = weight_del_numerator_dr / weight_del_denominator_dr

    # normalize the weights
    # total_simulation[t] = sim_count
    weight_gamma_dr = weight_gamma_dr / sum(weight_gamma_dr)
    weight_a_dr = weight_a_dr / sum(weight_a_dr)
    weight_m_dr= weight_m_dr / sum(weight_m_dr)
    weight_del_dr= weight_del_dr / sum(weight_del_dr)


    return weight_gamma_dr, weight_a_dr, weight_m_dr,weight_del_dr, propose_gamma_dr, propose_a_dr, propose_m_dr,propose_del_dr
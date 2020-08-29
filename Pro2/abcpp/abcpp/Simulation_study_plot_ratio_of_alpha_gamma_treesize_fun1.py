# file : boxplot_uppertri.py
def plotratio_type2(tree_no):
    import os
    import numpy as np
    import platform
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="whitegrid")
    check = True
    gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
    a_vec =  np.array([0,0.001,0.01,0.1,0.5,1])
    gamma_list = []
    a_list = []
    nv_list = []
    fit_list = []
    treesize_list = []
    ratio_list = []
    m_list = []
    gamma_index_list = []
    a_index_list = []
    label_a = (['$\\alpha$=0', '$\\alpha$=.001', '$\\alpha$=.01', '$\\alpha$=.1', '$\\alpha$=.5', '$\\alpha$=1'])
    label_gamma = (['$\gamma$=0', '$\gamma$=.001', '$\gamma$=.01', '$\gamma$=.1', '$\gamma$=.5', '$\gamma$=1'])
    count = 0
    if platform.system() == 'Windows':
        dire = 'C:/Liang/Googlebox'
    elif platform.system() == 'Darwin':
        dire = '/Users/dudupig/GoogleDrive'
    tree_file = dire + f"/Research/Project2/treesim_newexp/example{tree_no}/Ltable.csv"
    treesize_data = pd.read_csv(tree_file)
    extant_data = treesize_data[treesize_data['V4'] == -1]
    treesize = extant_data.shape[0]
    if treesize > 65:
        print('No valid parameter combinations.')
    else:
        for gindicator in range(1,6):
            for aindicator in range(1,6):
                generating_gamma = gamma_vec[gindicator]
                generating_alpha = a_vec[aindicator]
                r_check = np.sqrt(generating_alpha/generating_gamma)/treesize
                if r_check < 1.5 and r_check>0.5:

                    file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (tree_no,gindicator,aindicator)
                    if os.path.isfile(file):
                        para_data = np.load(file,allow_pickle=True).item()
                        generation = len(para_data['gamma'])
                        population = len(para_data['gamma'][0])
                        fitness = para_data['fitness'][generation - 1]
                        if check == True:
                            q5 = np.argsort(fitness)[-int(population // 4)]  # best 25%
                            fit_index = np.where(fitness > fitness[q5])[0]
                        else:
                            fit_index = np.array(range(population))
                        gamma_est = para_data['gamma'][generation - 1][fit_index]
                        a_est = para_data['a'][generation - 1][fit_index]
                        gamma_list.append(gamma_est)
                        a_list.append(a_est)

                        ratio_value = np.sqrt(a_est / gamma_est) / treesize
                        m_measure = np.abs(ratio_value-1) + np.sqrt(np.var(ratio_value))
                        ratio_list.append(np.sqrt(a_est / gamma_est) / treesize)
                        m_list.append(m_measure)
                        gamma_index_list.append(np.repeat(gamma_vec[gindicator], len(fit_index)))
                        a_index_list.append(np.repeat(a_vec[aindicator], len(fit_index)))

                    else:
                        pass
                else:
                    pass

        upper1 = 1.5
        lower1 = 0.5
        upper2 = 10
        lower2 = 0.05

        upper_m1 = 1
        upper_m2 = 2
        rclass = []
        gamma_con = np.concatenate(gamma_list, axis=0)
        a_con = np.concatenate(a_list, axis=0)
        ratio_con = np.concatenate(ratio_list, axis = 0)
        gamma_index_con = np.concatenate(gamma_index_list, axis = 0)
        a_index_con = np.concatenate(a_index_list, axis = 0)
        m_measure_con = np.concatenate(m_list, axis = 0)

        for i in range(len(ratio_con)):
            if (ratio_con[i] < upper1 and ratio_con[i] > lower1):
                rclass.append('3')
            elif (ratio_con[i] < upper2 and ratio_con[i] > upper1 or ratio_con[i] < lower1 and
                  ratio_con[i] > lower2):
                rclass.append('2')
            else:
                rclass.append('1')

        mclass = []
        for i in range(len(m_measure_con)):
            if (m_measure_con[i] < upper_m1):
                mclass.append('Close')
            elif (m_measure_con[i] < upper_m2 and m_measure_con[i] > upper_m1):
                mclass.append('Intermediate')
            else:
                mclass.append('Far')


        datalist = {'gamma_index': gamma_index_con, 'a_index': a_index_con,\
                    'gamma_est': gamma_con, 'a_est': a_con, 'ratio': ratio_con, 'rclass': rclass,\
                    'm_measure': m_measure_con, 'mclass': mclass}
        datadf = pd.DataFrame(datalist)
        subdf = datadf[datadf['rclass'] == 'Close']

        total_est = len(gamma_con)
        est_con = np.concatenate((gamma_con, a_con))
        ratio_con2 = np.concatenate((ratio_con,ratio_con))
        m_measure_con2 = np.concatenate((m_measure_con,m_measure_con))
        paras = np.concatenate((np.repeat('gamma', total_est), np.repeat('alpha', total_est)))
        rclass_con = np.concatenate((rclass,rclass))
        mclass_con = np.concatenate((mclass, mclass))
        total_index_con = np.concatenate((gamma_index_con, a_index_con))


        datalist2 = {'index': total_index_con, 'value': est_con, 'ratio': ratio_con2, 'rclass':
            rclass_con, 'parameter': paras, 'm_measure': m_measure_con2, 'mclass': mclass_con}
        datadf2 = pd.DataFrame(datalist2)
        subdf2 = datadf2[datadf2['rclass'] == 'Close']
        datadf2.mclass = datadf2.mclass.astype('category')
        # datadf2.mclass.cat.reorder_categories(['Close', 'Intermediate', 'Far'])

        # f, axes = plt.subplots(1, 1, figsize=(9, 9), sharex=True, sharey=True)
        g = sns.catplot(x="index", y="value", hue = "mclass", col="parameter",
                           data=datadf2, alpha = 0.3, dodge=True, palette={"Far": "r",
                                                                           "Intermediate": "tab:orange",
                        "Close": "g"}, legend=False)

        g.axes[0][0].set(xlabel='Generating parameters of $\\gamma$', ylabel='Estimated values of $\\gamma$', title = "")
        g.axes[0][1].set(xlabel='Generating parameters of $\\alpha$', ylabel='Estimated values of $\\alpha$', title = "")

        title = f"Scenario {tree_no} with the tree of {treesize} species"
        plt.subplots_adjust(top=0.7)
        g.fig.suptitle(title, x=0.45, y=1)  # can also get the figure from plt.gcf()
        g.fig.tight_layout()
        g.fig.subplots_adjust(right=0.85)

        legend = plt.legend(loc='upper right',frameon=False, bbox_to_anchor=(1.4, 0.5))
        legend_title = 'Measure class'
        legend.set_title(legend_title)
        # new_labels = ['Close', 'Intermediate', 'Far']
        # for t, l in zip(legend.texts, new_labels): t.set_text(l)


        dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/Ratio_type1_m2_sce%d.png' % (
            tree_no)

        f = g.fig
        f.savefig(dir_fig)
        plt.close(f)

if __name__ == '__main__':
    plotlist1 = [i for i in range(1,23)]
    for treevector in plotlist1:
        plotratio_type2(treevector)
        print(treevector)

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

    for gindicator in range(0,6):
        for aindicator in range(0,6):
            if platform.system() == 'Windows':
                dire = 'C:/Liang/Googlebox'
            elif platform.system() == 'Darwin':
                dire = '/Users/dudupig/GoogleDrive'
            file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (tree_no,gindicator,aindicator)
            tree_file = dire + f"/Research/Project2/treesim_newexp/example{tree_no}/Ltable.csv"
            if os.path.isfile(file):
                para_data = np.load(file,allow_pickle=True).item()
                treesize_data = pd.read_csv(tree_file)
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
                extant_data = treesize_data[treesize_data['V4'] == -1]
                treesize = extant_data.shape[0]
                ratio_value = np.sqrt(a_est / gamma_est) / treesize
                m_measure = np.abs(ratio_value-1) + np.sqrt(np.var(ratio_value))
                ratio_list.append(np.sqrt(a_est / gamma_est) / treesize)
                m_list.append(m_measure)
                gamma_index_list.append(np.repeat(gamma_vec[gindicator], len(fit_index)))
                a_index_list.append(np.repeat(a_vec[aindicator], len(fit_index)))

            else:
                pass

    upper1 = 1.5
    lower1 = 0.5
    upper2 = 10
    lower2 = 0.05

    upper_m1 = 0.5
    upper_m2 = 1
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
            mclass.append('3')
        elif (m_measure_con[i] < upper_m2 and m_measure_con[i] > upper_m1):
            mclass.append('2')
        else:
            mclass.append('1')


    datalist = {'gamma_index': gamma_index_con, 'a_index': a_index_con,\
                'gamma_est': gamma_con, 'a_est': a_con, 'ratio': ratio_con, 'rclass': rclass,\
                'm_measure': m_measure_con, 'mclass': mclass}
    datadf = pd.DataFrame(datalist)
    subdf = datadf[datadf['rclass'] == '3']

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
    subdf2 = datadf2[datadf2['rclass'] == '3']



    # f, axes = plt.subplots(1, 1, figsize=(9, 9), sharex=True, sharey=True)
    g = sns.catplot(x="index", y="value", hue = "mclass", col="parameter",
                       data=datadf2, alpha = 0.3, dodge=True, palette={"1": "r", "2": "tab:orange",
                    "3": "g"}, legend=False)

    g.axes[0][0].set(xlabel='Generating parameters', ylabel='Estimated values')
    g.axes[0][1].set(xlabel='Generating parameters', ylabel='')

    title = f"Scenario {tree_no} with the tree of {treesize} species"
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(title)  # can also get the figure from plt.gcf()


    legend = plt.legend(loc='upper right',frameon=False, bbox_to_anchor=(1, 1.1))
    legend_title = 'Measure class'
    legend.set_title(legend_title)
    new_labels = ['Far', 'Intermediate', 'Close']
    for t, l in zip(legend.texts, new_labels): t.set_text(l)
    # axes.set_title(title, fontsize = 20)
    # L = plt.legend(frameon=False,loc='upper right', bbox_to_anchor=(1, 1))
    # L.get_texts()[0].set_text('$\\gamma$')
    # L.get_texts()[1].set_text('$\\alpha$')

    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/Ratio_type1_m_sce%d.png' % (
        tree_no)

    f = g.fig
    f.savefig(dir_fig)
    plt.close(f)

if __name__ == '__main__':
    plotlist1 = [i for i in range(1,23)]
    for treevector in plotlist1:
        plotratio_type2(treevector)
        print(treevector)

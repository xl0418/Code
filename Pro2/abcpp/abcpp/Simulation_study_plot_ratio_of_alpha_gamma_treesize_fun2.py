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
                ratio_list.append(np.sqrt(a_est / gamma_est) / treesize)
                gamma_index_list.append(np.repeat(gamma_vec[gindicator], len(fit_index)))
                a_index_list.append(np.repeat(a_vec[aindicator], len(fit_index)))

            else:
                pass

    upper = 1.2
    lower = 0.8
    rclass = []
    gamma_con = np.concatenate(gamma_list, axis=0)
    a_con = np.concatenate(a_list, axis=0)
    ratio_con = np.concatenate(ratio_list, axis = 0)
    gamma_index_con = np.concatenate(gamma_index_list, axis = 0)
    a_index_con = np.concatenate(a_index_list, axis = 0)

    for i in range(len(ratio_con)):
        if (ratio_con[i] < upper and ratio_con[i] > lower):
            rclass.append('g')
        else:
            rclass.append('r')

    datalist = {'gamma_index': gamma_index_con, 'a_index': a_index_con,\
                'gamma_est': gamma_con, 'a_est': a_con, 'ratio': ratio_con, 'rclass': rclass}
    datadf = pd.DataFrame(datalist)
    subdf = datadf[datadf['rclass'] == 'g']

    total_est = len(gamma_con)
    est_con = np.concatenate((gamma_con, a_con))
    ratio_con2 = np.concatenate((ratio_con,ratio_con))
    paras = np.concatenate((np.repeat('gamma', total_est), np.repeat('alpha', total_est)))
    rclass_con = np.concatenate((rclass,rclass))
    total_index_con = np.concatenate((gamma_index_con, a_index_con))


    datalist2 = {'index': total_index_con, 'value': est_con, 'ratio': ratio_con2, 'rclass':
        rclass_con, 'parameter': paras}
    datadf2 = pd.DataFrame(datalist2)
    subdf2 = datadf2[datadf2['rclass'] == 'g']



    f, axes = plt.subplots(1, 1, figsize=(9, 9), sharex=True, sharey=True)
    sns.despine(top=True, right=True)
    sns.violinplot(x="index", y="value", hue = "parameter", data=subdf2, ax=axes, alpha = 0.5,
                  palette="Set2",  inner='quartile')
    axes.set(xlabel='Generating parameters', ylabel='Estimated values')
    title = f"Scenario {tree_no} with the tree of {treesize} species"
    axes.set_title(title, fontsize = 20)
    L = plt.legend(frameon=False,loc='upper right', bbox_to_anchor=(1, 1))
    L.get_texts()[0].set_text('$\\gamma$')
    L.get_texts()[1].set_text('$\\alpha$')

    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/Ratio_type3_sce%d.png' % (
        tree_no)

    f.savefig(dir_fig)
    plt.close(f)

if __name__ == '__main__':
    plotlist1 = [i for i in range(1,23)]
    for treevector in plotlist1:
        plotratio_type2(treevector)
        print(treevector)

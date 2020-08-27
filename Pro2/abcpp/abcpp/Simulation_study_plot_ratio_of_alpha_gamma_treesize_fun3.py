# file : boxplot_uppertri.py
def plotratio_type2(tree_no):
    import os
    import numpy as np
    import platform
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="white")
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

    upper = 1.5
    lower = 0.5
    rclass = []
    gamma_con = np.concatenate(gamma_list, axis=0) * treesize**2
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
    max_value = np.ceil(np.max(subdf['gamma_est']))
    x_sequence = np.linspace(0, max_value, num = 10)

    f, axes = plt.subplots(1, 1, figsize=(9, 9), sharex=True, sharey=True)
    sns.despine(top=True, right=True)
    axes = sns.scatterplot(x="gamma_est", y="a_est", hue="a_index" ,style="gamma_index", data=subdf,
                    alpha = 0.8, palette="Set2")
    axes.plot(x_sequence, x_sequence, color = 'red', linestyle = 'dashed',
              alpha = 0.5)
    axes.set(xlabel='$\\gamma$', ylabel='$\\alpha$')
    axes.get_xticks()
    xlabels = ["{:.4f}".format(x).rstrip('0').rstrip('.') for x in axes.get_xticks() / treesize**2]
    axes.set_xticklabels(xlabels)
    title = f"Scenario {tree_no} with the tree of {treesize} species"
    axes.set_title(title, fontsize = 20)
    h, l = axes.get_legend_handles_labels()
    num_of_a = len(subdf['a_index'].unique())+1
    l[0] = '$\\alpha^\star$'
    l[num_of_a] = '$\\gamma^\star$'
    color_hl = h[:num_of_a], l[:num_of_a]
    sizes_hl = h[num_of_a:], l[num_of_a:]

    l1 = axes.legend(*color_hl, loc='upper right', bbox_to_anchor = (.9, 0.3), frameon=False)
    l2 = axes.legend(*sizes_hl, loc='upper right', bbox_to_anchor = (1.05, 0.3), frameon=False)
    axes.add_artist(l1)

    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/Ratio_type5_1_sce%d.png' % (
        tree_no)

    f.savefig(dir_fig)
    plt.close(f)

if __name__ == '__main__':
    plotlist1 = [i for i in range(1,23)]
    for treevector in plotlist1:
        plotratio_type2(treevector)
        print(treevector)

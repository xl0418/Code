# file : boxplot_uppertri.py
def plotratio(treenovector):
    import os
    import numpy as np
    import platform
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="white")
    # treeno_vec = [i for i in range(4,7)]
    treeno_vec = treenovector
    standardize = 0
    # gno_vec = [0,0,0,0,1,1,1,2,2,3]
    # ano_vec = [2,3,4,5,3,4,5,4,5,5]
    gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
    a_vec =  np.array([0,0.001,0.01,0.1,0.5,1])
    gamma_list = []
    a_list = []
    nv_list = []
    fit_list = []
    label_a = (['$\\alpha$=0', '$\\alpha$=.001', '$\\alpha$=.01', '$\\alpha$=.1', '$\\alpha$=.5', '$\\alpha$=1'])
    label_gamma = (['$\gamma$=0', '$\gamma$=.001', '$\gamma$=.01', '$\gamma$=.1', '$\gamma$=.5', '$\gamma$=1'])

    truenv=1e-11
    count = 0
    ncol_fig = 6
    nrow_fig = 6
    f, axes = plt.subplots(nrow_fig, ncol_fig,figsize = (9,9) , sharex=True, sharey=True)

    # indicator = 5
    for gindicator in range(0,6):
        for aindicator in range(0,6):
            gamma_list = []
            a_list = []
            fit_list = []
            treesize_list = []
            figrow = count//ncol_fig
            figcol = count%ncol_fig
            for tree in treeno_vec:
                # g=gno_vec[indicator]
                # a=ano_vec[indicator]
                if platform.system() == 'Windows':
                    dire = 'C:/Liang/Googlebox'
                elif platform.system() == 'Darwin':
                    dire = '/Users/dudupig/GoogleDrive'
                file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (tree,gindicator,aindicator)
                tree_file = dire + f"/Research/Project2/treesim_newexp/example{tree}/Ltable.csv"
                if os.path.isfile(file):
                    para_data = np.load(file,allow_pickle=True).item()
                    treesize_data = pd.read_csv(tree_file)
                    generation = len(para_data['gamma'])
                    population = len(para_data['gamma'][0])
                    gamma_list.append(para_data['gamma'][generation - 1])
                    a_list.append(para_data['a'][generation - 1])
                    fit_list.append(para_data['fitness'][generation - 1])
                    treesize = treesize_data.shape[0]
                    treesize_list.append(np.repeat(treesize, population))
                else:
                    print('Fail')
                    gamma_list.append([-1 for i in range(population)])
                    a_list.append([-1 for i in range(population)])
                    fit_list.append([-1 for i in range(population)])
                    treesize_list.append([-1 for i in range(population)])


            gamma_t123 = np.concatenate(gamma_list, axis=0)
            a_t123 = np.concatenate(a_list, axis=0)
            treesize_t123 = np.concatenate(treesize_list, axis=0)
            ratio = np.sqrt(a_t123/gamma_t123)/treesize_t123

            if standardize == 1 and gindicator >0 and aindicator>0:
                gamma_t123 = gamma_t123/gamma_vec[gindicator]
                a_t123 = a_t123/a_vec[aindicator]

            test_label = (['S%d' % i for i in treeno_vec])
            test_df_label = np.repeat(test_label, population)
            test_dfex_label = np.tile(test_df_label, 3)
            gamma_label = np.repeat('$\gamma$', len(treeno_vec) * population)
            a_label = np.repeat('$\\alpha$', len(treeno_vec) * population)
            ratio_label = np.repeat('$Ratio$', len(treeno_vec) * population)
            variable_label = np.concatenate([gamma_label, a_label, ratio_label])
            value = np.concatenate([gamma_t123, a_t123, ratio])
            datalist = {'value': value, 'variable': variable_label, 'test': test_dfex_label}
            datadf = pd.DataFrame(datalist)


            ax = sns.boxplot(x="test", y="value",
                    hue="variable", palette=["m", "g", "r"],
                    data=datadf, linewidth=1, ax=axes[figrow,figcol],showfliers=False)

            if standardize == 1 and gindicator >0 and aindicator>0:
                axes[figrow, figcol].axhline(y=1, color='m', linestyle='--')
                axes[figrow, figcol].axhline(y=1, color='g', linestyle='--')
            else:
                axes[figrow,figcol].axhline(y=gamma_vec[gindicator],color='m',linestyle = '--')
                axes[figrow,figcol].axhline(y=a_vec[aindicator],color = 'g',linestyle = '--')
            axes[figrow,figcol].axhline(y=1,color = 'r',linestyle = '--')
            handles, labels = ax.get_legend_handles_labels()
            ax.legend_.remove()

            ax.set_ylabel('')
            ax.set_xlabel('')
            ax.set_xticks([])
            ax.set_xticklabels([])
            # ax.set_yticks([])
            # ax.set_yticklabels([])
            ax.set(ylim=(-.2, 2.2))
            sns.despine()
            if count in range(0, nrow_fig):
                axes[figrow,figcol].title.set_text(label_a[count])


            if count in ([5, 11, 17, 23, 29, 35]):
                axes[figrow,figcol].set_ylabel(label_gamma[int(count / nrow_fig)])
                axes[figrow,figcol].yaxis.set_label_position("right")

            # if figrow == figcol:
            #     axes[figrow,figcol].set_xticks( [ i for i in range(len(treeno_vec)) ] )
            #     axes[figrow, figcol].set_xticklabels(treeno_vec)
            #     axes[figrow, figcol].set_yticks([0, 1,2])
            #     axes[figrow, figcol].set_yticklabels(['0','1','2'])

                # axes[count].legend([])
            count += 1
            figcol += 1
        figrow += 1
    x_label = f"Scenario {tree} with the tree of {treesize} species"
    f.text(0.5, 0.04, x_label, ha='center',fontsize=15)
    f.text(0.04, 0.5, 'Estimates', va='center', rotation='vertical',fontsize=15)
    l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(0.75, 7.85), loc=2, borderaxespad=0.)
    l.get_frame().set_linewidth(0.0)
    if(len(treeno_vec) == 2):
        dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/CompareS%d_%d.png' % (treeno_vec[0],treeno_vec[1])
    elif len(treeno_vec) == 1:
        dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/Ratio_sce%d.png' % (
            treeno_vec[0])
    else:
        # dir_fig = 'C:/test.png'
        str=''
        for i in range(0, len(treeno_vec)):
            str = str + '%d' % treeno_vec[i]
        dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_ratio_test/scenarioS%s.png' % str

    f.savefig(dir_fig)
    plt.close(f)

if __name__ == '__main__':
    plotlist1 = [[i] for i in range(1,23)]
    for treevector in plotlist1:
        plotratio(treevector)
        print(treevector)

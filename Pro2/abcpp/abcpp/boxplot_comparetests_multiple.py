import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="white")
# treeno_vec = [i for i in range(4,7)]
treeno_vec = [1,2,3]

# gno_vec = [0,0,0,0,1,1,1,2,2,3]
# ano_vec = [2,3,4,5,3,4,5,4,5,5]
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec =  np.array([0,0.001,0.01,0.1,0.5,1])
gamma_list = []
a_list = []
nv_list = []
fit_list = []

truenv=1e-11
count = 0
ncol_fig = 6
nrow_fig = 6
f, axes = plt.subplots(nrow_fig, ncol_fig,figsize = (9,9), sharex=True, sharey=True)

# indicator = 5
for gindicator in range(0,6):
    for aindicator in range(0,6):
        gamma_list = []
        a_list = []
        nv_list = []
        fit_list = []
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
            if os.path.isfile(file):
                para_data = np.load(file).item()
                generation = len(para_data['gamma'])
                population = len(para_data['gamma'][0])
                gamma_list.append(para_data['gamma'][generation - 1])
                a_list.append(para_data['a'][generation - 1])
                nv_list.append(para_data['nu'][generation - 1])
                fit_list.append(para_data['fitness'][generation - 1])
            else:
                print('Fail')
                gamma_list.append([-1 for i in range(population)])
                a_list.append([-1 for i in range(population)])
                nv_list.append([-1 for i in range(population)])
                fit_list.append([-1 for i in range(population)])

        gamma_t123 = np.concatenate(gamma_list, axis=0)
        a_t123 = np.concatenate(a_list, axis=0)
        nv_t123 = np.concatenate(nv_list, axis=0)
        nv_t123 = nv_t123*1e11
        test_label = (['Test %d' % i for i in treeno_vec])
        test_df_label = np.repeat(test_label, population)
        test_dfex_label = np.tile(test_df_label, 3)
        gamma_label = np.repeat('$\gamma$', len(treeno_vec) * population)
        a_label = np.repeat('$\\alpha$', len(treeno_vec) * population)
        nv_label = np.repeat('$\\nu$,1e-11', len(treeno_vec) * population)
        variable_label = np.concatenate([gamma_label, a_label, nv_label])
        value = np.concatenate([gamma_t123, a_t123, nv_t123])
        datalist = {'value': value, 'variable': variable_label, 'test': test_dfex_label}
        datadf = pd.DataFrame(datalist)

        ax = sns.boxplot(x="test", y="value",
                    hue="variable", palette=["m", "g", "b"],
                    data=datadf, linewidth=1, ax=axes[figrow,figcol],showfliers=False)
        axes[figrow,figcol].axhline(y=gamma_vec[gindicator],color='m',linestyle = '--')
        axes[figrow,figcol].axhline(y=a_vec[aindicator],color = 'g',linestyle = '--')
        axes[figrow,figcol].axhline(y=1,color = 'b',linestyle = '--')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend_.remove()

        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set(ylim=(-.2, 2.2))
        sns.despine()
        # axes[count].legend([])
        count += 1
        figcol += 1
    figrow += 1


l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
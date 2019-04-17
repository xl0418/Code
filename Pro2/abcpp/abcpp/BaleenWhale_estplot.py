import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")

gamma_list = []
a_list = []
nv_list = []
fit_list = []
label_a = (['$\\alpha$=0', '$\\alpha$=.001', '$\\alpha$=.01', '$\\alpha$=.1', '$\\alpha$=.5', '$\\alpha$=1'])
label_gamma = (['$\gamma$=0', '$\gamma$=.001', '$\gamma$=.01', '$\gamma$=.1', '$\gamma$=.5', '$\gamma$=1'])


count = 0
timescale_vec = [20000,40000,80000]
heritability_vec = [1,2]
dividing_vec = [1,4]
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'
for timescaling_index in range(3):
    for dividing_index in range(2):
        for heritability_index in range(2):
            count += 1
            timescaling = timescale_vec[timescaling_index]
            heritability = heritability_vec[heritability_index]
            dividing = dividing_vec[dividing_index]
            #
            # data_name = data_dir + 'BWest_t%i_d%i_h%i.npy' % (int(timescaling),int(dividing),int(heritability_index))
            # est_data = np.load(data_name).item()
            # generation = len(est_data['gamma'])
            # population = len(est_data['gamma'][0])
            # gamma_list.append(est_data['gamma'][generation-1])
            # a_list.append(est_data['a'][generation-1])
            # nv_list.append(est_data['nu'][generation-1])

            # random test
            gamma_list.append(np.random.normal(count,1,20000))
            a_list.append(np.random.normal(count,1,20000))
            nv_list.append(np.random.normal(count,1,20000))


gamma_flat_list = [item for sublist in gamma_list for item in sublist]
a_flat_list = [item for sublist in a_list for item in sublist]
nv_flat_list = [item for sublist in nv_list for item in sublist]


timescaling_list = np.repeat(timescale_vec,4*20000)
dividing_list = np.tile(np.repeat(dividing_vec,2*20000),3)
heritability_list = np.tile(np.repeat(heritability_vec,20000),6)

est_list = {'gamma':gamma_flat_list,'alpha':a_flat_list,'nu':nv_flat_list,'timescale':timescaling_list,
            'dividing':dividing_list,'heritability':heritability_list}
est_df = pd.DataFrame(est_list)



f, axes = plt.subplots(2, 3,figsize = (9,12))
ax1 = sns.boxplot(x="timescale", y="gamma",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[0,0], showfliers=False)
ax2 = sns.boxplot(x="timescale", y="alpha",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[0,1], showfliers=False)
ax3 = sns.boxplot(x="timescale", y="nu",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[0,2], showfliers=False)
ax4 = sns.boxplot(x="timescale", y="gamma",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[1,0], showfliers=False)
ax5 = sns.boxplot(x="timescale", y="alpha",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[1,1], showfliers=False)
ax6 = sns.boxplot(x="timescale", y="nu",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[1,2], showfliers=False)



# indicator = 5
for image_count in range(6):
    ax = sns.boxplot(x="day", y="total_bill",
                     hue="size", palette=["m", "g", "b"],
                     data=tips, linewidth=1, ax=axes[figrow, figcol], showfliers=False)        test_label = (['S%d' % i for i in treeno_vec])



        test_df_label = np.repeat(test_label, population)
        test_dfex_label = np.tile(test_df_label, 3)
        gamma_label = np.repeat('$\gamma$', len(treeno_vec) * population)
        a_label = np.repeat('$\\alpha$', len(treeno_vec) * population)
        nv_label = np.repeat('$\\nu$,$10^{-11}$', len(treeno_vec) * population)
        variable_label = np.concatenate([gamma_label, a_label, nv_label])
        value = np.concatenate([gamma_t123, a_t123, nv_t123])
        datalist = {'value': value, 'variable': variable_label, 'test': test_dfex_label}
        datadf = pd.DataFrame(datalist)

        if gindicator>aindicator:
            axes[figrow, figcol].plot([])
            axes[figrow, figcol].axis('off')
        else:
            ax = sns.boxplot(x="test", y="value",
                    hue="variable", palette=["m", "g", "b"],
                    data=datadf, linewidth=1, ax=axes[figrow,figcol],showfliers=False)

            if standardize == 1 and gindicator >0 and aindicator>0:
                axes[figrow, figcol].axhline(y=1, color='m', linestyle='--')
                axes[figrow, figcol].axhline(y=1, color='g', linestyle='--')
            else:
                axes[figrow,figcol].axhline(y=gamma_vec[gindicator],color='m',linestyle = '--')
                axes[figrow,figcol].axhline(y=a_vec[aindicator],color = 'g',linestyle = '--')
            axes[figrow,figcol].axhline(y=1,color = 'b',linestyle = '--')
            handles, labels = ax.get_legend_handles_labels()
            ax.legend_.remove()

            ax.set_ylabel('')
            ax.set_xlabel('')
            # ax.set_xticks([])
            # ax.set_xticklabels([])
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

f.text(0.84, 0.04, 'Scenarios', ha='center',fontsize=15)
f.text(0.04, 0.83, 'Estimates', va='center', rotation='vertical',fontsize=15)
l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(0.75, 7.85), loc=2, borderaxespad=0.)
l.get_frame().set_linewidth(0.0)
if(len(treeno_vec) == 2):
    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_estuppertri/CompareS%d_%d.png' % (treeno_vec[0],treeno_vec[1])
elif len(treeno_vec) == 1:
    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_estuppertri/EstS%d.png' % (treeno_vec[0])
else:
    # dir_fig = 'C:/test.png'
    str=''
    for i in range(0, len(treeno_vec)):
        str = str + '%d' % treeno_vec[i]
    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_estuppertri/CompareS%s.png' % str

f.savefig(dir_fig)
plt.close(f)




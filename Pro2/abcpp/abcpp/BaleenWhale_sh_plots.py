import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dendropy
from dendropy.model import continuous
import csv
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'
files = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'
# calculate the empirical trait and pics

with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(files+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = logTL


# reorder the traits according to the simulation order
sim_species_label = [ "Balaena_mysticetus","Balaenoptera_acutorostrata" ,"Caperea_marginata" ,    "Balaenoptera_borealis" ,
  "Balaenoptera_physalus" ,     "Eschrichtius_robustus"  ,    "Balaenoptera_musculus"    ,  "Balaenoptera_omurai"    ,
 "Eubalaena_australis" ,    "Megaptera_novaeangliae" , "Balaenoptera_bonaerensis"  , "Balaenoptera_brydei"   ,
 "Balaenoptera_edeni"     ,    "Eubalaena_glacialis"   ,  "Eubalaena_japonica" ]

obsZ_ordered_sim = length[[np.where(sim_species_label[i] == extantlabels_array)[0][0] for i in range(15)]]
obsZ = obsZ_ordered_sim

# PIC calculation
taxa1 = dendropy.TaxonNamespace()
dataset_combined = dendropy.DataSet.get(path=files+"bw_char.nex",schema="nexus")
tree_emp = dataset_combined.tree_lists[0][0]
chars_emp = dataset_combined.char_matrices[0]
pic_emp = continuous.PhylogeneticIndependentConstrasts(tree=tree_emp,
char_matrix=chars_emp)
ctree_emp = pic_emp.contrasts_tree(character_index=0,
annotate_pic_statistics=True,state_values_as_node_labels=False,corrected_edge_lengths=False)
emp_pic = []
label = []
for nd in ctree_emp.postorder_internal_node_iter():
    emp_pic.append(nd.pic_contrast_standardized)
    label.append(int(nd.label))

emp_pic_orded_node = abs(np.array(emp_pic)[np.argsort(label)])



s0h0_list = []
s0h1_list = []
s1h0_list = []
s1h1_list = []
s2h0_list = []
s2h1_list = []
s3h0_list = []
s3h1_list = []

for data_folder in ['results_0930_smtd/', 'results_0930_umtd_pics/', 'results_0930_pics/']:
    data_name = data_dir + data_folder+'trait_pics_df.npy'
    df = np.load(data_name,allow_pickle=True).item()
    trait = df['trait']
    pics = df['pics']
    scenarios = df['scenario']
    num_samples = df['pics'].shape[0]
    fitness = np.zeros([1,num_samples])
    if data_folder == 'results_0930_smtd/':
        #GOF: Goodness of fit
        s = np.argsort(obsZ)
        obsZ = obsZ[s]
        i, j = argsort2D(trait)
        trait = trait[i, j]
        fitness[0,:] += 1.0 - normalized_norm(trait, obsZ)
    elif data_folder == 'results_0930_umtd_pics/':
        fitness[0,:] += 1.0 - normalized_norm(trait, obsZ)
        fitness[0,:] += 1.0 - normalized_norm(emp_pic_orded_node, pics)
    else:
        fitness[0,:] += 1.0 - normalized_norm(emp_pic_orded_node, pics)

    q5 = np.argsort(fitness)[0][-int(num_samples // 4)]  # best 25%
    fit_index = np.where(fitness[0] > fitness[0][q5])[0]

    s0h0 = len(np.where(scenarios[fit_index] == 0)[0]) / len(fit_index)
    s0h1 = len(np.where(scenarios[fit_index] == 1)[0]) / len(fit_index)
    s1h0 = len(np.where(scenarios[fit_index] == 2)[0]) / len(fit_index)
    s1h1 = len(np.where(scenarios[fit_index] == 3)[0]) / len(fit_index)
    s2h0 = len(np.where(scenarios[fit_index] == 4)[0]) / len(fit_index)
    s2h1 = len(np.where(scenarios[fit_index] == 5)[0]) / len(fit_index)
    s3h0 = len(np.where(scenarios[fit_index] == 6)[0]) / len(fit_index)
    s3h1 = len(np.where(scenarios[fit_index] == 7)[0]) / len(fit_index)

    s0h0_list.append(s0h0)
    s0h1_list.append(s0h1)
    s1h0_list.append(s1h0)
    s1h1_list.append(s1h1)
    s2h0_list.append(s2h0)
    s2h1_list.append(s2h1)
    s3h0_list.append(s3h0)
    s3h1_list.append(s3h1)
#
#
# r = [0,1,2]
#
# # plot
# barWidth = 0.85
# names = ('SMTD', 'UMTD+PICs', 'PICs')
# # Create green Bars
# h0_colors = ['#81C7D4','#33A6B8','#0089A7','#005CAF']
# h1_colors = ['#B28FCE','#986DB2','#77428D','#4A225D']
#
# bar1 = plt.bar(r, s0h0_list, color=h0_colors[0], edgecolor='white', width=barWidth)
# # Create orange Bars
# bar2 = plt.bar(r, s0h1_list, bottom=s0h0_list, color=h1_colors[0], edgecolor='white',
#                width=barWidth)
# # Create blue Bars
# bar3 = plt.bar(r, s1h0_list, bottom=[i + j for i, j in zip(s0h1_list, s0h0_list)],
#                color=h0_colors[1],
#         edgecolor='white', width=barWidth)
#
# bar4 = plt.bar(r, s1h1_list, bottom=[i + j +k for i, j,k in zip(s0h1_list, s0h0_list,s1h0_list)],
#         color=h1_colors[1],
#         edgecolor='white', width=barWidth)
#
# bar5 = plt.bar(r, s2h0_list, bottom=[i + j +k+m for i, j,k,m in zip(s0h1_list, s0h0_list,s1h0_list,
#                                                              s1h1_list)],
#         color=h0_colors[2],
#         edgecolor='white', width=barWidth)
#
# bar6 = plt.bar(r, s2h1_list, bottom=[i + j +k+m+n for i, j,k,m,n in zip(s0h1_list, s0h0_list,
#                                                                         s1h0_list,
#                                                              s1h1_list,s2h0_list)],
#         color=h1_colors[2],
#         edgecolor='white', width=barWidth)
#
# bar7 = plt.bar(r, s3h0_list, bottom=[i + j +k+m+n+h for i, j,k,m,n,h in zip(s0h1_list, s0h0_list,
#                                                                             s1h0_list,
#                                                              s1h1_list,s2h0_list,s2h1_list)],
#         color=h0_colors[3],
#         edgecolor='white', width=barWidth)
#
# bar8 = plt.bar(r, s3h1_list, bottom=[i + j +k+m+n+h+u for i, j,k,m,n,h,u in zip(s0h1_list,
#                                                                                 s0h0_list,
#                                                                       s1h0_list,
#                                                              s1h1_list,s2h0_list,s2h1_list,s3h0_list)],
#         color=h1_colors[3],
#         edgecolor='white', width=barWidth)
# # Custom x axis
# plt.xticks(r, names)
# # plt.xlabel("group")
# plt.legend((bar1[0],bar2[0],bar3[0],bar4[0],bar5[0],bar6[0],bar7[0],bar8[0]),
#            ('s0h0','s0h1','s1h0','s1h1','s2h0','s2h1','s3h0','s3h1'))
# plt.legend(ncol=8, bbox_to_anchor=(0, 1),
#               loc='lower left', fontsize='small')
# # Show graphic
# plt.show()


category_names = ['$s=20000,h^2=0.5$','$s=20000,h^2=1$',
                  '$s=40000,h^2=0.5$','$s=40000,h^2=1$',
                  '$s=60000,h^2=0.5$','$s=60000,h^2=1$',
                  '$s=80000,h^2=0.5$','$s=80000,h^2=1$']
results = {
    'SMTD': [s0h0_list[0], s0h1_list[0], s1h0_list[0],  s1h1_list[0], s2h0_list[0],s2h1_list[0],
             s3h0_list[0],s3h1_list[0]],
    'UMTD+PICs': [s0h0_list[1], s0h1_list[1],  s1h0_list[1], s1h1_list[1], s2h0_list[1],
                  s2h1_list[1],
             s3h0_list[1],s3h1_list[1]],
    'PICs': [s0h0_list[2], s0h1_list[2],  s1h0_list[2], s1h1_list[2],s2h0_list[2],s2h1_list[2],
             s3h0_list[2],s3h1_list[2]]

}

def survey(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(9.2, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        # xcenters = starts + widths / 2

        # r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        # for y, (x, c) in enumerate(zip(xcenters, widths)):
        #     ax.text(x, y, str(int(c)), ha='center', va='center',
        #             color=text_color)
    ax.legend(ncol=4, bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return fig, ax


survey(results, category_names)
plt.show()
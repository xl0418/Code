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

    q5 = np.argsort(fitness)[0][-int(num_samples // 20)]  # best 25%
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


data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster' \
           '/results_ms_1028/'
tvp_gof = []
tv_gof = []
tvm_gof = []
for datafile in ['smtd_con40.npy', 'umtd_con40.npy', 'pics_con40.npy']:
    data_name = data_dir + datafile

    est_data = np.load(data_name,allow_pickle=True).item()
    population = int(len(est_data['model_data'][0]) / 3)
    fitness = est_data['fitness'][-1]
    total_population = len(est_data['model_data'][0])
    model_index = np.array([0, 1, 2])
    model_params = np.repeat(model_index, repeats=population)
    propose_model = model_params
    q5 = np.argsort(fitness)[-int(total_population // 20)]  # best 25%
    fit_index = np.where(fitness > fitness[q5])[0]

    modelTVPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
    modelTVperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)
    modelTVMperc = len(np.where(propose_model[fit_index] == 2)[0]) / len(fit_index)

    tvp_gof.append(modelTVPperc)
    tv_gof.append(modelTVperc)
    tvm_gof.append(modelTVMperc)



model_names = ['AWC','UWC','MWC']
models = {
    'SMTD': [tvp_gof[0], tv_gof[0], tvm_gof[0]],
    'UMTD+PICs': [tvp_gof[1], tv_gof[1],  tvm_gof[1]],
    'PICs': [tvp_gof[2], tv_gof[2],  tvm_gof[2]]
}


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
labels_1 = list(results.keys())
data_1 = np.array(list(results.values()))
data_cum_1 = data_1.cumsum(axis=1)
category_colors_1 = plt.get_cmap('RdYlGn')(
    np.linspace(0.15, 0.85, data_1.shape[1]))

labels_2 = list(models.keys())
data_2 = np.array(list(models.values()))
data_cum_2 = data_2.cumsum(axis=1)
category_colors_2 = plt.get_cmap('rainbow')(
    np.linspace(0.15, 0.85, data_2.shape[1]))



fig, axes = plt.subplots(nrows = 2, ncols = 1,figsize=(9, 9))
axes[0].invert_yaxis()
axes[0].xaxis.set_visible(False)
axes[0].set_xlim(0, np.sum(data_1, axis=1).max())

axes[1].invert_yaxis()
axes[1].xaxis.set_visible(False)
axes[1].set_xlim(0, np.sum(data_2, axis=1).max())

for i, (colname, color) in enumerate(zip(category_names, category_colors_1)):
    widths = data_1[:, i]
    starts = data_cum_1[:, i] - widths
    axes[0].barh(labels_1, widths, left=starts, height=0.5,
            label=colname, color=color)

axes[0].legend(ncol=4, bbox_to_anchor=(0, 1),
          loc='lower left', fontsize='small')

colors_three_models = ['#EE442F', '#601A4A','#85C0F9']
for i, (colname, color) in enumerate(zip(model_names, category_colors_2)):
    widths = data_2[:, i]
    starts = data_cum_2[:, i] - widths
    axes[1].barh(labels_2, widths, left=starts, height=0.5,
            label=colname, color=colors_three_models[i])

axes[1].legend(ncol=3, bbox_to_anchor=(0, 1),
          loc='lower left', fontsize='small')


plt.show()
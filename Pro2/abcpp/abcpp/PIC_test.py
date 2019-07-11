import sys
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import csv
import dendropy
from dendropy.model import continuous
from multiprocessing import Pool
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')

# gamma = 0.001
# a = 0.1
#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
num_cores = Pool(8)  # the number of cores

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err
    # calculate pic of the simulated traits

#full tree
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
    #'/home/p274981/abcpp/'

files = dir_path + 'treedata/'



time_scalar = 20000
heri_sqr = 1

td = DVTreeData(path=files, scalar=time_scalar)


K=10e5
nu=1e-3

file_result = dir_path + 'BaleenWhales/BWest_t%i_h%f.npy' % (int(time_scalar),heri_sqr)

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
length = 10**logTL


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

tree_sim = dendropy.Tree.get(
    path=files + "bw.nex", schema="nexus",
    taxon_namespace=taxa1)
tree_sim.print_plot()
simchar_dict = {}
keys = ["B.mysticetus", "B.acutorostrata", "C.marginata", "B.borealis",
        "B.physalus", "E.robustus", "B.musculus", "B.omurai",
        "E.australis", "M.novaeangliae", "B.bonaerensis", "B.brydei",
        "B.edeni", "E.glacialis", "E.japonica"]

for i in range(15):
    simchar_dict[keys[i]] = [i]

simchars = dendropy.ContinuousCharacterMatrix.from_dict(simchar_dict, taxon_namespace=taxa1)
simpic = continuous.PhylogeneticIndependentConstrasts(tree=tree_sim, char_matrix=simchars)
sim_ctree = simpic.contrasts_tree(character_index=0,
                                  annotate_pic_statistics=True,
                                  state_values_as_node_labels=False,
                                  corrected_edge_lengths=False)
sim_pic = []
sim_label = []
for nd in sim_ctree.postorder_internal_node_iter():
    sim_pic.append(nd.pic_contrast_raw)
    sim_label.append(int(nd.label))
sim_pic_ordered = np.array(sim_pic)[np.argsort(sim_label)]

import numpy as np
import dendropy
from dendropy.model import continuous
from multiprocessing import Pool
from itertools import repeat
from pic_compute import pic_compute

data_dir = '/home/p274981/abcpp/abcpp/'

# PIC calculation
#full tree

files = '/home/p274981/abcpp/BaleenWhales/treedata/'

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

num_cores = Pool(10)
tree_sim = dendropy.Tree.get(
    path=files + "bw.nex", schema="nexus",
    taxon_namespace=taxa1)

for data_folder in ['results_0930_smtd/', 'results_0930_umtd/', 'results_0930_pics/']:
    data_name = data_dir + data_folder
    file_result = data_dir+data_folder+ 'trait_pics_df.npy'
    scenario = 0
    trait_matrix = np.zeros([1,15])
    sce_vec = np.array([0])
    for data_file in ['BWest_t20000_h0', 'BWest_t20000_h1', 'BWest_t40000_h0','BWest_t40000_h1',
                      'BWest_t60000_h0', 'BWest_t60000_h1', 'BWest_t80000_h0', 'BWest_t80000_h1']:
        data = data_name + data_file + '.npy'
        est_data = np.load(data,allow_pickle=True).item()
        valid = est_data['valid']
        trait = est_data['Z']
        trait_matrix = np.concatenate([trait_matrix,trait],axis = 0)
        sce_vec = np.concatenate([sce_vec,np.repeat(scenario,len(valid))])
        scenario += 1

    trait_matrix = np.delete(trait_matrix,0,0)
    sce_vec = sce_vec[1:]
    num_valid_sims = len(sce_vec)
    pic_ordered_list = num_cores.starmap(pic_compute,
                                         zip(repeat(tree_sim), trait, repeat(taxa1),
                                             range(num_valid_sims)))

    # in case the parallel computation returns disordered output
    order_list = []
    contrast_list = []
    for i in range(num_valid_sims):
        order_list.append(pic_ordered_list[i][1])
        contrast_list.append(pic_ordered_list[i][0])
    ordered_contrast_list = [contrast_list[item] for item in np.argsort(order_list)]
    contrast_array = abs(np.vstack(ordered_contrast_list))

    trait_pics_df = {'trait':trait_matrix, 'pics':contrast_array}
    np.save(file_result, trait_pics_df)


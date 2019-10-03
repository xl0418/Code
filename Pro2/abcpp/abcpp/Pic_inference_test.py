import sys, os
sys.path.append('C:/Liang/abcpp_ms8/abcpp')
from Dvtraitsim_TVM import DVSimTVMLog10, DVSimTVM
from Dvtraitsim_TVP import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import numpy as np
import dendropy
from dendropy.model import continuous
scalar = 20000

#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

td = DVTreeData(path=files, scalar=scalar)


# Test TVP
theta = 3  # optimum of natural selection

gamma =5# intensity of natural selection
r = 1  # growth rate
a = 120  # intensity of competition
K = 10e5  # carrying capacity
nu=1.715e-3
Vmax = 1e-4

# parameter settings
obs_param_TVP = DVParamLiang(gamma=gamma, a=a, K=K,h=0.5, nu=nu, r=1, theta=theta,V00=1e-5,
                             V01=1e-5,
                          Vmax=Vmax, inittrait=theta, initpop=500,
                initpop_sigma = 10.0, break_on_mu=False)

population = 100

obs_param_TVP_batch = np.tile(obs_param_TVP, (population, 1))
simmodelTVP = dvcpp.DVSimTVP(td, obs_param_TVP_batch)
valid_TVP = np.where(simmodelTVP['sim_time'] == td.sim_evo_time)[0]
if len(simmodelTVP['Z'][valid_TVP])>0:
    print(gamma,a)
    simmodelTVP['Z'][valid_TVP].max()
    simmodelTVP['Z'][valid_TVP].min()
# else:
#     print('No valid simulations')
sorder_index =np.argsort(simmodelTVP['Z'][0])
simmodelTVP['V'][0][sorder_index]*1e5



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
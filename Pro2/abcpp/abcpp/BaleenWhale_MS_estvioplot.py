import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import seaborn as sns

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/result_0617/'

obs_file = dir_path + 'treedata/'

with open(obs_file+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
gamma_list = []
a_list = []
nu_list = []
vm_list = []
distance_list = []
ratio_dis = []
valid_length = []
timescaling_list = []
heri_list = []
count = 0
particle_size = 100
K=1e6
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
obsZ = sorted(10**logTL)

# for timescaling_index in range(len(timescale_vec)):
#     for heritability_index in range(len(heritability_vec)):
# print(count)
# count += 1
# timescaling = timescale_vec[timescaling_index]
# heritability = heritability_vec[heritability_index]

data_name = data_dir + 'modelselec4w.npy'
est_data = np.load(data_name).item()
population = int(len(est_data['model_data'][0])/3)
fitness = est_data['fitness'][-1]

q5_TVP = np.argsort(fitness[ :population])[-int(population // 200+1)]  # best 5%
q5_TV = np.argsort(fitness[ population:2*population])[-int(population // 200+1)]+population  # best 5%
q5_TVM = np.argsort(fitness[ 2*population:])[-int(population // 200+1)]+2*population  # best 5%

fit_index_TVP = np.where(fitness[ :population] > fitness[ q5_TVP])[0]
fit_index_TV = np.where(fitness[ population:2*population] > fitness[ q5_TV])[0]+population
fit_index_TVM = np.where(fitness[ 2*population:] > fitness[ q5_TVM])[0]+2*population

previous_bestfitted_index_TVP = fit_index_TVP
previous_bestfitted_index_TV = fit_index_TV - population
previous_bestfitted_index_TVM = fit_index_TVM- 2*population

top_pop = len(previous_bestfitted_index_TVP)

# last_fitness = est_data['fitness'][generation - 1]
# q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
# fit_index = np.where(last_fitness > last_fitness[q5])[0]

# mean of all samples
# gamma_mean = np.mean(est_data['gamma'][generation-1])
# a_mean = np.mean(est_data['a'][generation-1])
# nv_mean = np.mean(est_data['nu'][generation-1])
# mean of the top 5% samples

gamma_vec_TVP = est_data['gamma_data_TVP'][-2,previous_bestfitted_index_TVP]*1e8
a_vec_TVP = est_data['a_data_TVP'][-2,previous_bestfitted_index_TVP]*1e6
nu_vec_TVP = est_data['nu_data_TVP'][-2,previous_bestfitted_index_TVP]*1e4
vm_vec_TVP = est_data['vm_data_TVP'][-2,previous_bestfitted_index_TVP]
theta_vec_TVP = est_data['theta_data_TVP'][-2,previous_bestfitted_index_TVP]


gamma_vec_TV = est_data['gamma_data_TV'][-2,previous_bestfitted_index_TV]*1e8
a_vec_TV = est_data['a_data_TV'][-2,previous_bestfitted_index_TV]*1e6
nu_vec_TV = est_data['nu_data_TV'][-2,previous_bestfitted_index_TV]*1e4
vm_vec_TV = est_data['vm_data_TV'][-2,previous_bestfitted_index_TV]
theta_vec_TV = est_data['theta_data_TV'][-2,previous_bestfitted_index_TV]


gamma_vec_TVM = est_data['gamma_data_TVM'][-2,previous_bestfitted_index_TVM]*1e8
a_vec_TVM = est_data['a_data_TVM'][-2,previous_bestfitted_index_TVM]*1e6
nu_vec_TVM = est_data['nu_data_TVM'][-2,previous_bestfitted_index_TVM]*1e4
vm_vec_TVM = est_data['vm_data_TVM'][-2,previous_bestfitted_index_TVM]
theta_vec_TVM = est_data['theta_data_TVM'][-2,previous_bestfitted_index_TVM]



print('TVP 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVP)*1e-8,
    np.mean(a_vec_TVP)*1e-6, np.mean(nu_vec_TVP)*1e-4,
    np.mean(vm_vec_TVP), np.mean(theta_vec_TVP)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVP)*1e-8,
    np.var(a_vec_TVP)*1e-6, np.var(nu_vec_TVP)*1e-4,
    np.var(vm_vec_TVP), np.var(theta_vec_TVP)))

print('TV 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TV)*1e-8,
    np.mean(a_vec_TV)*1e-6, np.mean(nu_vec_TV)*1e-4,
    np.mean(vm_vec_TV), np.mean(theta_vec_TV)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TV)*1e-8,
    np.var(a_vec_TV)*1e-6, np.var(nu_vec_TV)*1e-4,
    np.var(vm_vec_TV), np.var(theta_vec_TV)))


print('TVM 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVM)*1e-8,
    np.mean(a_vec_TVM)*1e-6, np.mean(nu_vec_TVM)*1e-4,
    np.mean(vm_vec_TVM), np.mean(theta_vec_TVM)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVM)*1e-8,
    np.var(a_vec_TVM)*1e-6, np.var(nu_vec_TVM)*1e-4,
    np.var(vm_vec_TVM), np.var(theta_vec_TVM)))


est_para = ['gamma','alpha','nu','vm','theta'] # ,'vm'
model_para = ['TVP','TV','TVM']
est_array = np.concatenate([gamma_vec_TVP,a_vec_TVP,nu_vec_TVP,vm_vec_TVP,theta_vec_TVP,
                            gamma_vec_TV,a_vec_TV,nu_vec_TV,vm_vec_TV,theta_vec_TV,
                            gamma_vec_TVM,a_vec_TVM,nu_vec_TVM,vm_vec_TVM,theta_vec_TVM]) # ,vm_list_flat
est_label = np.tile(np.repeat(est_para,top_pop),len(model_para))
model_label = np.repeat(model_para,top_pop*len(est_para))

ss_list = {'est':est_array,'est_label':est_label,
           'model_label':model_label}
ss_df = pd.DataFrame(ss_list)

vioplot = sns.catplot(x="model_label", y="est",col="est_label",
     data=ss_df, kind="box",height=5, aspect=.6,sharey=False)
# vioplot.set(ylim=(-50,400))
# vioplot.set_xticklabels(["$\gamma \cdot 10^{-8}$", "$\\alpha \cdot 10^{-5}$", "$\\nu \cdot 10^{-4}$","$V_{m}$"])
vioplot.set_axis_labels("", "Estimate value")
# vioplot._legend.set_title('$h^2$')
axes = vioplot.axes.flatten()
axes[0].set_title("$\gamma \cdot 10^{-8}$")
axes[1].set_title("$\\alpha \cdot 10^{-6}$")
axes[2].set_title("$\\nu \cdot 10^{-4}$")
axes[3].set_title("$V_m $")
axes[4].set_title("$\\theta $")



Z = est_data['Z']
diff_norm = np.linalg.norm(Z - obsZ, axis=1)
plt.hist(diff_norm,bins=100)
TVP_index = np.where(fitness[:population]>0)[0]
TV_index = np.where(fitness[population:2*population]>0)[0]+population
TVM_index = np.where(fitness[2*population:]>0)[0]+2*population

plt.hist(diff_norm[:len(TVP_index)])
plt.hist(diff_norm[len(TVP_index):len(TVP_index)+len(TV_index)])
plt.hist(diff_norm[len(TVP_index)+len(TV_index):])

import sys
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
td = DVTreeData(path=obs_file, scalar=20000)

largeini = np.where(theta_vec_TVM>2000)[0]
param_TVM = DVParamLiang(gamma=1, a=1, K=1e12, h=1, nu=1, r=1, theta=1,
                         V00=.5,V01=.5, Vmax=1, inittrait=1300, initpop=1e5,
                     initpop_sigma=10.0, break_on_mu=False)
params = np.tile(param_TVM,(len(largeini),1))
params[:,0]=gamma_vec_TVM[largeini]*1e-8
params[:,1]=a_vec_TVM[largeini]*1e-6
params[:,4]=nu_vec_TVM[largeini]*1e-4
params[:,6]=theta_vec_TVM[largeini]
params[:,9]=vm_vec_TVM[largeini]

predictsim = dvcpp.DVSimTVM(td, params)
valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]
Z = predictsim['Z'][valid]
i, j = argsort2D(Z)
Z_modelTVM = Z[i, j]
diff_norm = np.linalg.norm(Z_modelTVM - obsZ, axis=1)

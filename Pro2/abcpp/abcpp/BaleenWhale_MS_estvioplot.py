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
# data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster' \
#            '/results_ms_1028/'
# data_name = data_dir + 'modelselec2w_umtd_pics.npy'

data_dir = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/results_ms_posterior/'
data_name = data_dir + 'smtd_con40.npy'

obs_file = dir_path + 'treedata/'

with open(obs_file + 'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path + 'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:, 1]
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
K = 1e6
nu = 1 / (100 * K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:, 0] == label)[0][0])

logTL = lengthdata_array[length_index, 1].astype(np.float)
obsZ = sorted(10 ** logTL)

# for timescaling_index in range(len(timescale_vec)):
#     for heritability_index in range(len(heritability_vec)):
# print(count)
# count += 1
# timescaling = timescale_vec[timescaling_index]
# heritability = heritability_vec[heritability_index]

est_data = np.load(data_name,allow_pickle=True).item()
population = int(len(est_data['model_data'][0]) / 3)
fitness = est_data['fitness'][-1]

q5_TVP = np.argsort(fitness[:population])[-int(population // 200 + 1)]  # best 5%
q5_TV = np.argsort(fitness[population:2 * population])[
            -int(population // 200 + 1)] + population  # best 5%
q5_TVM = np.argsort(fitness[2 * population:])[-int(population // 200 + 1)] + 2 * population  #
# best 5%

fit_index_TVP = np.where(fitness[:population] > fitness[q5_TVP])[0]
fit_index_TV = np.where(fitness[population:2 * population] > fitness[q5_TV])[0] + population
fit_index_TVM = np.where(fitness[2 * population:] > fitness[q5_TVM])[0] + 2 * population

previous_bestfitted_index_TVP = fit_index_TVP
previous_bestfitted_index_TV = fit_index_TV - population
previous_bestfitted_index_TVM = fit_index_TVM - 2 * population

top_pop = [len(previous_bestfitted_index_TVP),len(previous_bestfitted_index_TV),
           len(previous_bestfitted_index_TVM)]

gamma_vec_TVP = est_data['gamma_data_TVP'][-2, previous_bestfitted_index_TVP]
a_vec_TVP = est_data['a_data_TVP'][-2, previous_bestfitted_index_TVP]
nu_vec_TVP = est_data['nu_data_TVP'][-2, previous_bestfitted_index_TVP]
vm_vec_TVP = est_data['vm_data_TVP'][-2, previous_bestfitted_index_TVP]
theta_vec_TVP = est_data['theta_data_TVP'][-2, previous_bestfitted_index_TVP]
ratio_TVP = np.sqrt(a_vec_TVP/gamma_vec_TVP)

gamma_vec_TV = est_data['gamma_data_TV'][-2, previous_bestfitted_index_TV]
a_vec_TV = est_data['a_data_TV'][-2, previous_bestfitted_index_TV]
nu_vec_TV = est_data['nu_data_TV'][-2, previous_bestfitted_index_TV]
vm_vec_TV = est_data['vm_data_TV'][-2, previous_bestfitted_index_TV]
theta_vec_TV = est_data['theta_data_TV'][-2, previous_bestfitted_index_TV]
ratio_TV = np.sqrt(a_vec_TV/gamma_vec_TV)

gamma_vec_TVM = est_data['gamma_data_TVM'][-2, previous_bestfitted_index_TVM]
a_vec_TVM = est_data['a_data_TVM'][-2, previous_bestfitted_index_TVM]
nu_vec_TVM = est_data['nu_data_TVM'][-2, previous_bestfitted_index_TVM]
vm_vec_TVM = est_data['vm_data_TVM'][-2, previous_bestfitted_index_TVM]
theta_vec_TVM = est_data['theta_data_TVM'][-2, previous_bestfitted_index_TVM]
ratio_TVM = np.sqrt(a_vec_TVM/gamma_vec_TVM)

print('TVP 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVP) ,
    np.mean(a_vec_TVP) , np.mean(nu_vec_TVP) ,
    np.mean(vm_vec_TVP), np.mean(theta_vec_TVP)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVP) ,
    np.var(a_vec_TVP) , np.var(nu_vec_TVP) ,
    np.var(vm_vec_TVP), np.var(theta_vec_TVP)))

print('TV 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TV) ,
    np.mean(a_vec_TV) , np.mean(nu_vec_TV),
    np.mean(vm_vec_TV), np.mean(theta_vec_TV)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TV) ,
    np.var(a_vec_TV) , np.var(nu_vec_TV) ,
    np.var(vm_vec_TV), np.var(theta_vec_TV)))

print('TVM 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.mean(gamma_vec_TVM) ,
    np.mean(a_vec_TVM) , np.mean(nu_vec_TVM) ,
    np.mean(vm_vec_TVM), np.mean(theta_vec_TVM)))
print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f theta = %f' % (
    np.var(gamma_vec_TVM),
    np.var(a_vec_TVM) , np.var(nu_vec_TVM) ,
    np.var(vm_vec_TVM), np.var(theta_vec_TVM)))

est_para = ['gamma', 'alpha','ratio', 'nu', 'vm', 'theta']  # ,'vm'
model_para = ['AWC', 'UWC', 'MWC']
est_array = np.concatenate([gamma_vec_TVP, a_vec_TVP,ratio_TVP, nu_vec_TVP, vm_vec_TVP,
                            theta_vec_TVP,
                            gamma_vec_TV, a_vec_TV,ratio_TV, nu_vec_TV, vm_vec_TV, theta_vec_TV,
                            gamma_vec_TVM, a_vec_TVM,ratio_TVM, nu_vec_TVM, vm_vec_TVM,
                            theta_vec_TVM])
est_label = np.repeat(np.tile(est_para,len(model_para)),np.tile(top_pop,len(est_para)))
model_label = np.repeat(model_para, np.array(top_pop) * len(est_para))

ss_list = {'est': est_array, 'est_label': est_label,
           'model_label': model_label}
ss_df = pd.DataFrame(ss_list)

vioplot = sns.catplot(x="model_label", y="est", col="est_label",col_wrap=3,palette=["#CB1B45",
                                                                                    "#FAD689",
                                                                          "#0D5661"],
                      data=ss_df, kind="box", height=5, aspect=.6, sharey=False)
# vioplot.set(ylim=(-50,400))
# vioplot.set_xticklabels(["$\gamma \cdot 10^{-8}$", "$\\alpha \cdot 10^{-5}$", "$\\nu \cdot 10^{
# -4}$","$V_{m}$"])
vioplot.set_axis_labels("", "Estimate value")
# vioplot._legend.set_title('$h^2$')
axes = vioplot.axes.flatten()
axes[0].set_title("$\gamma $")
axes[1].set_title("$\\alpha $")
axes[2].set_title("$\sqrt{\\alpha/\gamma} $")

axes[3].set_title("$\\nu $")
axes[4].set_title("$V_m $")
axes[5].set_title("$\\theta $")

for no_plots in range(0,5):
    axes[no_plots].ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)


#
# GOF distribution of TVP TV TVM
Z = est_data['Z']
diff_norm = np.linalg.norm(Z - obsZ, axis=1)
plt.hist(diff_norm, bins=200)
TVP_index = np.where(fitness[:population] > 0)[0]
TV_index = np.where(fitness[population:2 * population] > 0)[0] + population
TVM_index = np.where(fitness[2 * population:] > 0)[0] + 2 * population
plot_tvm_index = np.where(fitness[TVM_index] > np.min(fitness[TVM_index]))[0] + 2 * population

plt.hist(diff_norm[:len(TVP_index)])
plt.hist(diff_norm[len(TVP_index):(len(TVP_index) + len(TV_index))])
plt.hist(diff_norm[(len(TVP_index) + len(TV_index)):])

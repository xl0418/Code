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
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'

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

timescale_vec = [20000,40000,60000,80000,100000]
heritability_vec = [0.5,1]


for timescaling_index in range(5):
    for heritability_index in range(2):
        print(count)
        count += 1
        timescaling = timescale_vec[timescaling_index]
        heritability = heritability_vec[heritability_index]

        data_name = data_dir + 'BWest_t%i_h%i.npy' % (int(timescaling), int(heritability_index))
        est_data = np.load(data_name).item()
        generation = len(est_data['gamma'])
        population = len(est_data['gamma'][0])
        fitness = est_data['fitness'][-1]
        valid = est_data['valid']
        q5 = np.argsort(fitness)[-len(valid) // 20]  # best 5%
        fit_index = np.where(fitness > fitness[q5])[0]


        # last_fitness = est_data['fitness'][generation - 1]
        # q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
        # fit_index = np.where(last_fitness > last_fitness[q5])[0]

        # mean of all samples
        # gamma_mean = np.mean(est_data['gamma'][generation-1])
        # a_mean = np.mean(est_data['a'][generation-1])
        # nv_mean = np.mean(est_data['nu'][generation-1])
        # mean of the top 5% samples

        gamma_vec = est_data['gamma'][-2]*1e7
        a_vec = est_data['a'][-2]*1e5
        nu_vec = est_data['nu'][-2]*1e5
        vm_vec = est_data['vm'][-2]
        inflow_mutation = 2*K*nu_vec*1e-5 * vm_vec/(1+4*K*nu_vec*1e-5)

        print('s=%i h2=%f 5th gamma = %.3e  a = %.3e nu = %.3e Vm = %f inflow = %f' % (timescaling,
                                                                                        heritability,
            np.mean(gamma_vec[fit_index])*1e-7,
            np.mean(a_vec[fit_index])*1e-5, np.mean(nu_vec[fit_index])*1e-5,
            np.mean(vm_vec[fit_index]), np.mean(inflow_mutation[fit_index])))
        print('Var: gamma = %.3e  a = %.3e nu = %.3e Vm = %f inflow = %f' % (
            np.var(gamma_vec[fit_index])*1e-7,
            np.var(a_vec[fit_index])*1e-5, np.var(nu_vec[fit_index])*1e-5,
            np.var(vm_vec[fit_index]), np.var(inflow_mutation[fit_index])))

        gamma_list.append(gamma_vec.tolist())
        a_list.append(a_vec.tolist())
        nu_list.append(nu_vec.tolist())
        vm_list.append(vm_vec.tolist())
        timescaling_list.append(np.repeat(timescaling,population))
        heri_list.append(np.repeat(heritability,population))



gamma_list_flat = np.array([item for sublist in gamma_list for item in sublist])
a_list_flat = np.array([item for sublist in a_list for item in sublist])
nu_list_flat = np.array([item for sublist in nu_list for item in sublist])
vm_list_flat = np.array([item for sublist in vm_list for item in sublist])

est_para = ['gamma','alpha','nu','vm'] # ,'vm'
est_array = np.concatenate([gamma_list_flat,a_list_flat,nu_list_flat,vm_list_flat]) # ,vm_list_flat
est_label = np.repeat(est_para,population*len(heritability_vec)*len(timescale_vec))
timescaling_list_flat = np.tile(np.repeat(timescale_vec,population*len(heritability_vec)),len(est_para))
heri_list_flat = np.tile(np.repeat(heritability_vec,population),len(est_para)*len(timescale_vec))


ss_list = {'est':est_array,'est_label':est_label,
           's':timescaling_list_flat,
            'heri':heri_list_flat}
ss_df = pd.DataFrame(ss_list)

sns.catplot(x="est_label", y="est", hue="heri", col="s",
     data=ss_df, kind="box", height=4, aspect=.7)

#
# g = sns.FacetGrid(ss_df, col="s",  row="heri",margin_titles=True)
# g.map(sns.boxplot, "est_label", "est")
# g.set_xlabels('$\gamma $ $(10^{-7})$')
# g.set_ylabels('$\\alpha $ $(10^{-4})$')
# axes = g.axes.flatten()
# axes[0].set_title("$h^{2}=0.5$")
# axes[1].set_title("$h^{2}=1$")

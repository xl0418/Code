import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
sys.path.append('C:/Liang/abcpp_ms/abcpp')
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
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
nv_list = []
distance_list = []
valid_length = []
timescaling_list = []
dividing_list = []
heritability_list =[]
count = 0
particle_size = 1000
K=10e8
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)

timescale_vec = [200,400,800]
heritability_vec = [1,0.5]
dividing_vec = [1,4]
for timescaling_index in range(3):
    for dividing_index in range(2):
        for heritability_index in range(2):
            print(count)
            count += 1
            timescaling = timescale_vec[timescaling_index]
            heritability = heritability_vec[heritability_index]
            dividing = dividing_vec[dividing_index]

            # data_name = data_dir + 'BWest_t%i_d%i_h%i.npy' % (int(timescaling),int(dividing),int(heritability_index))
            # est_data = np.load(data_name).item()
            # generation = len(est_data['gamma'])
            # population = len(est_data['gamma'][0])
            # gamma_list.append(np.mean(est_data['gamma'][generation-1]))
            # a_list.append(np.mean(est_data['a'][generation-1]))
            # nv_list.append(np.mean(est_data['nu'][generation-1]))

            # random test
            gamma_est = count*0.00001
            a_est = count*0.001
            nv_est = count*10**(-11)
            length = 10 ** logTL / dividing
            obsZ = length
            meantrait = np.mean(obsZ)
            td = DVTreeData(path=obs_file, scalar=timescaling)

            param = DVParamLiang(gamma=0.0, a=0.5, K=K,h=np.sqrt(heritability), nu=nu, r=1, theta=meantrait, V00 = .5, V01=.5, Vmax=1, inittrait=meantrait, initpop=1e5,
             initpop_sigma = 10.0, break_on_mu=False)
            params = np.tile(param, (particle_size, 1))  # duplicate

            predictsim = dvcpp.DVSimLiang(td, params)

            valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]

            Z = predictsim['Z'][valid]
            i, j = argsort2D(Z)
            Z = Z[i, j]
            # V = pop['V'][valid][i, j]
            Z = np.nan_to_num(Z)
            tp_distance = np.linalg.norm(Z - obsZ, axis=1)
            distance_list.append(tp_distance)
            timescaling_list.append(np.repeat(timescaling,len(valid)))
            dividing_list.append(np.repeat(dividing,len(valid)))
            heritability_list.append(np.repeat(heritability,len(valid)))


distance_list_flat = [item for sublist in distance_list for item in sublist]
timescaling_list_flat = [item for sublist in timescaling_list for item in sublist]
dividing_list_flat = [item for sublist in dividing_list for item in sublist]
heritability_list_flat = [item for sublist in heritability_list for item in sublist]


ss_list = {'distance':distance_list_flat,'timescale':timescaling_list_flat,
            'dividing':dividing_list_flat,'heritability':heritability_list_flat}
ss_df = pd.DataFrame(ss_list)


f, axes = plt.subplots(1, 2,figsize = (9,12),sharex=True,sharey=True)
ax1 = sns.boxplot(x="timescale", y="distance",
                 hue="heritability", palette=["m", "g"],
                 data=ss_df[lambda ss_df: ss_df.dividing == 1],
                  linewidth=1, ax=axes[0], showfliers=False)
ax1.title.set_text('Dividing scalar: 1')
ax2 = sns.boxplot(x="timescale", y="distance",
                 hue="heritability", palette=["m", "g"],
                 data=ss_df[lambda ss_df: ss_df.dividing == 4],
                  linewidth=1, ax=axes[1], showfliers=False)

ax2.title.set_text('Dividing scalar: 4')

handles_gamma, labels_gamma = ax1.get_legend_handles_labels()

for ax in [ax1,ax2]:
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.legend_.remove()

axes[0].set_ylabel('Distance')
axes[0].yaxis.set_label_position("left")


handles = handles_gamma #[ item for subhandle in [handles_gamma,handles_a,handles_nu] for item in subhandle]
labels = ['$h^2 = 1$', '$h^2 = 0.5$']
f.text(0.5, 0.04, 'Time scaling parameters', ha='center',fontsize=15)
# f.text(0.04, 0.5, 'Estimates', va='center', rotation='vertical',fontsize=15)
l = plt.legend(handles, labels, bbox_to_anchor=(0.85, 3.7), loc=2, borderaxespad=0.)
l.get_frame().set_linewidth(0.0)

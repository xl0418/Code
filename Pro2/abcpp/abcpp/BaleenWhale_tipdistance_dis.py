import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import seaborn as sns

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)

dir_path ='C:/Liang/Googlebox/Research/Project2/BaleenWhales/'


obs_file = dir_path + 'treedata/'
sim_file_TVP = dir_path + 'result_cluster/Est/predictsimTVP.csv'
sim_file_TV = dir_path + 'result_cluster/Est/predictsimTV.csv'
sim_file_TVM = dir_path + 'result_cluster/Est/predictsimTVM.csv'

with open(obs_file+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)


extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = 10**logTL
obsZ = sorted(length)


with open(sim_file_TVP) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TVP = list(csv2_reader)[1:]

with open(sim_file_TV) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TV = list(csv2_reader)[1:]

with open(sim_file_TVM) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TVM = list(csv2_reader)[1:]

trait_TVP_array = np.array(trait_TVP,dtype=float)
i, j = argsort2D(trait_TVP_array)
trait_TVP_array = trait_TVP_array[i, j]

trait_TV_array = np.array(trait_TV,dtype=float)
i, j = argsort2D(trait_TV_array)
trait_TV_array = trait_TV_array[i, j]

trait_TVM_array = np.array(trait_TVM,dtype=float)
i, j = argsort2D(trait_TVM_array)
trait_TVM_array = trait_TVM_array[i, j]

obsZ_distance = np.linalg.norm(obsZ)

tvp_distance = np.linalg.norm(trait_TVP_array - obsZ, axis=1)/obsZ_distance
tv_distance = np.linalg.norm(trait_TV_array - obsZ, axis=1)/obsZ_distance
tvm_distance = np.linalg.norm(trait_TVM_array - obsZ, axis=1)/obsZ_distance

# Make a separate list for each airline
x1 = list(tvp_distance)
x2 = list(tv_distance)
x3 = list(tvm_distance)

# Assign colors for each airline and the names
colors = ['red', 'green','blue']
names = ['TVP','TV','TVM']

# Make the histogram using a list of lists
# Normalize the flights and assign colors and names
plt.hist([x1, x2,x3], bins=20, color=colors, label=names)
plt.legend()
plt.xlabel('Relatvie Distance')
plt.ylabel('Frequence')
plt.title('Fitness measure')


def gap(targetlist):
    return np.array([targetlist[i+1]-targetlist[i] for i in range(len(targetlist)-1)])

gap_tvp = np.apply_along_axis(gap, 1, trait_TVP_array)
gap_tv =np.apply_along_axis(gap, 1, trait_TV_array)
gap_tvm =np.apply_along_axis(gap, 1, trait_TVM_array)

gap_obs = gap(obsZ)

gap_tvp_flatten = gap_tvp.flatten()
gap_tv_flatten = gap_tv.flatten()
gap_tvm_flatten = gap_tvm.flatten()

gap_whole = np.concatenate([gap_tvp_flatten,gap_tv_flatten,gap_tvm_flatten])
gap_label = np.tile(range(14),3000)
model_label = np.repeat(['TVP','TV','TVM'],14000)
df_list = {'gap':gap_whole, 'gap_label':gap_label,'model':model_label}
df_pd = pd.DataFrame(df_list)
species_gap = range(14)

fig, axes = plt.subplots(1, 3,sharey=True,sharex=True)
axes = axes.flatten()


ax_tvp = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model']=='TVP'],
                 showfliers=False,ax= axes[0],color="skyblue")
axes[0].scatter(species_gap,gap_obs,c='red')
ax_tvp.set_title('TVP')

ax_tv = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model']=='TV'],
                 showfliers=False,ax= axes[1],color="skyblue")
axes[1].scatter(species_gap,gap_obs,c='red')
ax_tv.set_title('TV')

ax_tvm = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model']=='TVM'],
                 showfliers=False,ax= axes[2],color="skyblue")
axes[2].scatter(species_gap,gap_obs,c='red')
ax_tvm.set_title('TVM')

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv


dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

obs_file = dir_path + 'treedata/'
sim_file = dir_path + '50itertest/BWMS1q_d40_fm1.npy'

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
length = 10**logTL/40
obsZ = length
assert os.path.isfile(sim_file),"%s doesn't exist!" % sim_file
para_data = np.load(sim_file).item()
valid = para_data['valid']
trait = para_data['Z']

tp_valid = np.where(valid<1000)[0]
dr_valid = np.where(valid>=1000)[0]

tp_distance = np.linalg.norm(trait[tp_valid,:] - obsZ, axis=1)
dr_distance = np.linalg.norm(trait[dr_valid,:] - obsZ, axis=1)


# Make a separate list for each airline
x1 = list(tp_distance)
x2 = list(dr_distance)

# Assign colors for each airline and the names
colors = ['#E69F00', '#56B4E9']
names = ['TP','DR']

# Make the histogram using a list of lists
# Normalize the flights and assign colors and names
plt.hist([x1, x2], bins=50, color=colors,normed=True, label=names)

# Plot formatting
plt.legend()
plt.xlabel('Distance')
plt.ylabel('Frequence')
plt.title('Fitness measure')
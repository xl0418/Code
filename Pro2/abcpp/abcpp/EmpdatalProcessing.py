import numpy as np
import csv
import pandas as pd
#full tree
empiricaldata = 'BW'
assert empiricaldata=='BW' or empiricaldata=='AS', 'Pls specify the empirical data...'
if empiricaldata == 'BW':
    dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

    files = dir_path + 'treedata/'
else:
    dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

    files = dir_path + 'treedata/'


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
length = 10**logTL/40
trait_labels_df = {'species': extantlabels_array,'TL':length}
tldf = pd.DataFrame(trait_labels_df)
tldf.to_csv(files+"trait_labels.csv", sep=',',index=False)




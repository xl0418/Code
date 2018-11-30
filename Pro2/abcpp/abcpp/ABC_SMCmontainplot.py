#%%
import os
import numpy as np
import platform
import pandas as pd

# treeno_vec = [i for i in range(4,7)]
tree = 2
gindicator = 1
aindicator = 5


if platform.system() == 'Windows':
    dire = 'C:/Liang/Googlebox'
elif platform.system() == 'Darwin':
    dire = '/Users/dudupig/GoogleDrive'
file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (tree,gindicator,aindicator)
assert os.path.isfile(file),"%s doesn't exist!" % file
para_data = np.load(file).item()

para_data_gamma_df = pd.DataFrame(para_data['gamma'])
para_data_a_df = pd.DataFrame(para_data['a'])

filesmcg_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdatag%d.csv' % tree
filesmca_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdataa%d.csv' % tree

para_data_gamma_df.to_csv(filesmcg_name, encoding='utf-8', index=False)
para_data_a_df.to_csv(filesmca_name, encoding='utf-8', index=False)

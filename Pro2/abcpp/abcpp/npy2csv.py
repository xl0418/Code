#%%
import os
import numpy as np
import platform
import pandas as pd

# treeno_vec = [i for i in range(4,7)]
test = 4
gindicator = 4
aindicator = 5


if platform.system() == 'Windows':
    dire = 'C:/Liang/Googlebox'
elif platform.system() == 'Darwin':
    dire = '/Users/dudupig/GoogleDrive'
file = dire + '/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (test,gindicator,aindicator)
assert os.path.isfile(file),"%s doesn't exist!" % file
para_data = np.load(file).item()

para_data_gamma_df = pd.DataFrame(para_data['gamma'])
para_data_a_df = pd.DataFrame(para_data['a'])

filesmcg_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdatag%d.csv' % test
filesmca_name = 'C:\\Liang\\PhdIntroProject2\\smcplot\\smcdataa%d.csv' % test

para_data_gamma_df.to_csv(filesmcg_name, encoding='utf-8', index=False)
para_data_a_df.to_csv(filesmca_name, encoding='utf-8', index=False)



# For model selection
iterations = 20
generating = 'TP'
fileMS = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/modelsele%s.npy' % generating
fileMS = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/BaleenWhalesMS.npy'

assert os.path.isfile(fileMS),"%s doesn't exist!" % fileMS
ms_data = np.load(fileMS).item()
modeldata = ms_data['model_data']
fitness = ms_data['fitness']
total_population = modeldata.shape[1]
bestpercent_index = int(total_population // 10)
bestmodel = np.zeros(shape=[iterations,bestpercent_index-1])
for g in range(iterations):
    q5 = np.argsort(fitness[g, :])[-bestpercent_index]  # best 25%
    fit_index = np.where(fitness[g, :] > fitness[g, q5])[0]
    bestmodel[g,:] = modeldata[g,fit_index]


bestmodel_df = pd.DataFrame(bestmodel)


modeldata_df = pd.DataFrame(ms_data['model_data'])
fitness_df = pd.DataFrame(ms_data['fitness'])

filesmcms_name = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/modeldata%s.csv' % generating
filesmcfit_name = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/modelfit%s.csv' % generating
# bestmodel_name = 'C:/Liang/Googlebox/Research/Project2/modelsele/example1/bestmodel%s.csv' % generating
bestmodel_name = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/BaleenWhalesMS.csv'
modeldata_df.to_csv(filesmcms_name, encoding='utf-8', index=False)
fitness_df.to_csv(filesmcfit_name, encoding='utf-8', index=False)
bestmodel_df.to_csv(bestmodel_name, encoding='utf-8', index=False)
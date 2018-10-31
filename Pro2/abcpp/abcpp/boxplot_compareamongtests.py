import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sns.set(style="white")
treeno_vec = [4,5,6]
gno_vec = [0,0,0,0,1,1,1,2,2,3]
ano_vec = [2,3,4,5,3,4,5,4,5,5]
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec =  np.array([0,0.001,0.01,0.1,0.5,1])
gamma_list = []
a_list = []
nv_list = []
fit_list = []

truenv=1e-11
indicator = 6

for tree in treeno_vec:
    g=gno_vec[indicator]
    a=ano_vec[indicator]
    if platform.system() == 'Windows':
        dire = 'C:/Liang/'
    elif platform.system() == 'Darwin':
        dire = '/Users/dudupig/'
    file = dire + 'Googlebox/Research/Project2/smc_newdata/test%d/smc%dg%da.npy' % (tree,g,a)
    if os.path.isfile(file):
        para_data = np.load(file).item()
        generation = len(para_data['gamma'])
        population = len(para_data['gamma'][0])
        gamma_list.append(para_data['gamma'][generation - 1])
        a_list.append(para_data['a'][generation - 1])
        nv_list.append(para_data['nu'][generation - 1])
        fit_list.append(para_data['fitness'][generation - 1])
    else:
        print('Fail')
        gamma_list.append([0])
        a_list.append([0])
        nv_list.append([0])
        fit_list.append([0])

gamma_t123 = np.concatenate(gamma_list, axis=0)
a_t123 = np.concatenate(a_list, axis=0)
nv_t123 = np.concatenate(nv_list, axis=0)
test_label = (['Test %d' % i for i in treeno_vec])
test_df_label = np.repeat(test_label, population)
test_dfex_label = np.tile(test_df_label, len(treeno_vec))
gamma_label = np.repeat('$\gamma$', len(treeno_vec) * population)
a_label = np.repeat('$\\alpha$', len(treeno_vec) * population)
nv_label = np.repeat('$\\nu$', len(treeno_vec) * population)
variable_label = np.concatenate([gamma_label, a_label, nv_label])
value = np.concatenate([gamma_t123, a_t123, nv_t123])
datalist = {'value': value, 'variable': variable_label, 'test': test_dfex_label}
datadf = pd.DataFrame(datalist)

sns.boxplot(x="test", y="value",
            hue="variable", palette=["m", "g", "b"],
            data=datadf, linewidth=3)
plt.axhline(y=gamma_vec[g],color='m',linestyle = '--')
plt.axhline(y=a_vec[a],color = 'g',linestyle = '--')
sns.despine(offset=10, trim=True)
plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'

smtd_gof = []
umtd_gof = []
pics_gof = []
for data_folder in ['results_0930_smtd/', 'results_0930_umtd/', 'results_0930_pics/']:
    data_name = data_dir + data_folder
    for data_file in ['BWest_t20000_h0', 'BWest_t20000_h1', 'BWest_t40000_h0','BWest_t40000_h1',
                      'BWest_t60000_h0', 'BWest_t60000_h1', 'BWest_t80000_h0', 'BWest_t80000_h1']:
        data = data_name + data_file + '.npy'
        est_data = np.load(data,allow_pickle=True).item()
        valid = est_data['valid']
        trait = est_data['Z'][valid,:]

        population = int(len(est_data['model_data'][0]) / 3)
        fitness = est_data['fitness'][-1]
        total_population = len(est_data['model_data'][0])
        model_index = np.array([0, 1, 2])
        model_params = np.repeat(model_index, repeats=population)
        propose_model = model_params
        q5 = np.argsort(fitness)[-int(total_population // 4)]  # best 25%
        fit_index = np.where(fitness > fitness[q5])[0]

        modelTVPperc = len(np.where(propose_model[fit_index] == 0)[0]) / len(fit_index)
        modelTVperc = len(np.where(propose_model[fit_index] == 1)[0]) / len(fit_index)
        modelTVMperc = len(np.where(propose_model[fit_index] == 2)[0]) / len(fit_index)

        tvp_gof.append(modelTVPperc)
        tv_gof.append(modelTVperc)
        tvm_gof.append(modelTVMperc)

r = [0,1,2]

# plot
barWidth = 0.85
names = ('SMTD', 'UMTD+PICs', 'PICs')
# Create green Bars
plt.bar(r, tvp_gof, color='#b5ffb9', edgecolor='white', width=barWidth)
# Create orange Bars
plt.bar(r, tv_gof, bottom=tvp_gof, color='#f9bc86', edgecolor='white', width=barWidth)
# Create blue Bars
plt.bar(r, tvm_gof, bottom=[i + j for i, j in zip(tvp_gof, tv_gof)], color='#a3acff',
        edgecolor='white', width=barWidth)

# Custom x axis
plt.xticks(r, names)
plt.xlabel("group")

# Show graphic
plt.show()

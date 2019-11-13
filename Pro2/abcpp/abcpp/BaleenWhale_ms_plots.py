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


data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster' \
           '/results_ms_1028/'
tvp_gof = []
tv_gof = []
tvm_gof = []
for datafile in ['smtd_con40.npy', 'umtd_con40.npy', 'modelselec2w_pics.npy']:
    data_name = data_dir + datafile

    est_data = np.load(data_name,allow_pickle=True).item()
    population = int(len(est_data['model_data'][0]) / 3)
    fitness = est_data['fitness'][-1]
    total_population = len(est_data['model_data'][0])
    model_index = np.array([0, 1, 2])
    model_params = np.repeat(model_index, repeats=population)
    propose_model = model_params
    q5 = np.argsort(fitness)[-int(total_population // 20)]  # best 25%
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
tvp = plt.bar(r, tvp_gof, color='#F17C67', edgecolor='white', width=barWidth)
# Create orange Bars
tv = plt.bar(r, tv_gof, bottom=tvp_gof, color='#FAD689', edgecolor='white', width=barWidth)
# Create blue Bars
tvm = plt.bar(r, tvm_gof, bottom=[i + j for i, j in zip(tvp_gof, tv_gof)], color='#2EA9DF',
        edgecolor='white', width=barWidth)

# Custom x axis
plt.xticks(r, names)
plt.legend((tvp,tv,tvm),
           ('AWC','UWC','MWC'),loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True)
# Show graphic
plt.show()





model_names = ['AWC','UWC','MWC']
models = {
    'SMTD': [tvp_gof[0], tv_gof[0], tvm_gof[0]],
    'UMTD+PICs': [tvp_gof[1], tv_gof[1],  tvm_gof[1]],
    'PICs': [tvp_gof[2], tv_gof[2],  tvm_gof[2]]
}
def survey(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(9.2, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        # xcenters = starts + widths / 2

        # r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        # for y, (x, c) in enumerate(zip(xcenters, widths)):
        #     ax.text(x, y, str(int(c)), ha='center', va='center',
        #             color=text_color)
    ax.legend(ncol=3, bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return fig, ax


survey(models, model_names)
plt.show()
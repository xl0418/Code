import sys
sys.path.append('C:/Liang/Code/Pro2/Python_p2')
from DV_model_sim_along_phy import DVtraitsim_tree
from ABC_SMC_DVmodel import calibration
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec

file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'
count = 0
junk = np.zeros((36,100))
for i in range(len(gamma_vec)):
    for j in range(len(a_vec)):
        for loop in range(1,101):
            simresult = DVtraitsim_tree(file, replicate=loop, theta=0, gamma1=gamma_vec[i], r=1, a=a_vec[j], scalar=1000, K=100000, nu=0.0001,
                                Vmax=1)
            if simresult[2]:
                junk[count,loop-1] = 1
            print('count %d; loop %d' % (count,loop))
        count += 1

junk_ratio = 1- np.sum(junk,axis=1)/100
junk_ratio.reshape((6,6))

piror = [0,1,0,1]
samplesize = 5000
K_vec = [100000,1000000,10000000]
for i in range(len(K_vec)):
    datafile = file + 'cali%dK' % (i+1)
    cali= calibration(samplesize=samplesize, priorpar=piror,K=K_vec[i], treefile=file,calidata_file=datafile)

plt.figure(figsize=(5,15))
for i in range(len(K_vec)):
    datafile = file + 'cali%dK' % (i+1)
    cali1 = np.load(datafile+'.npz')
    plt.subplot(3, 1, i+1)
    plt.title(('K=%d' % K_vec[i]))
    if i in range(2):
        plt.xticks([])
    calipara1 = cali1['calipar']
    # sns.jointplot(calipara1[:, 0], calipara1[:, 1], kind="kde", stat_func=None, color="#4CB391")
    sns.kdeplot(calipara1[:, 0],calipara1[:, 1], cmap="Greens", shade=True, shade_lowest=True )


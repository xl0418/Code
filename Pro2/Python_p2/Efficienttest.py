from DV_model_sim_along_phy import DVtraitsim_tree
import timeit

starttime = timeit.default_timer()
file = 'C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\'

for loop in range(0,100):
    DVtraitsim_tree(file, replicate = loop,theta = 0, gamma1 = 0.01, r = 1, a = 0.1,scalar = 1000, K = 100000, nu = 0.0001, Vmax = 1)
    print(loop)
endtime = timeit.default_timer()
elapse = endtime - starttime
timetext = 'Elapsed time: %.2f' % elapse
print(timetext)
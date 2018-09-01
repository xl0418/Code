import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/Code/Pro2/Python_p2')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from DV_model_sim_along_phy import DVtraitsim_tree
if platform.system()=='Windows':
    file = 'C:\\Liang\\Code\\Pro2\\abcpp\\tree_data\\example1\\'
elif platform.system()=='Darwin':
    file = '/Users/dudupig/Documents/GitHub/Code/Pro2/abcpp/tree_data/example3/'
gamma1=0
a=1
scalor=10000
K=1000000000
nu=1/(100*K)

simresult = DVtraitsim_tree(file=file, gamma1=gamma1, a=a, K=K, scalar=scalor,nu=nu,Vmax=1)
simresult[3]

for r in range(1000):
    print(r)
    simresult = DVtraitsim_tree(file=file, gamma1=gamma1, a=a, K=K, nu=nu, scalar=scalor)
    if simresult[2]:
        pic = 0
        break
    else:
        pic = 1
print(pic)
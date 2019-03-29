import sys
import platform
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from PhyloDiff_model_sim import Candimodels

if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_emp/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_shared import DVTreeData, DVParam

def Candi_para(para):
    dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
    files = dir_path + 'treedata/'
    scalar = 1000
    td = DVTreeData(path=files, scalar=scalar)
    return Candimodels(td=td,param=para)

import sys
from multiprocessing import Pool
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from ParallelComOrderFun import square
import numpy as np
num_cores = Pool(2)  # the number of cores



input = np.array([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0])

pic_ordered_list = num_cores.starmap(square, zip(input))

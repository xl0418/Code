import sys
import platform
import numpy as np
sys.path.append('C:/Liang/Code/Pro2/abcpp/abcpp/')
from PhyloDiff_model_sim import Candimodels

paramstest = np.array([0.0, 0.0, 0, 1.0, 0, 1.0])
testsim = Candimodels(td,paramstest)
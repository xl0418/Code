#!/bin/bash
module load icc		# intel c++ compiler 
module load Python  # intel's python3 
icc -ffast-math -fopenmp -O3 -march=native -shared -fPIC -std=c++14 ./module/module.cpp -I./module/pybind11/include `python3-config --includes --cflags --ldflags` -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`

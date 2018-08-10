#!/bin/bash
g++ -ffast-math -fopenmp -O3 -mtune=native -shared -fPIC -std=c++14 ./module/module.cpp -I./module/pybind11/include `python3-config --includes --cflags --ldflags` -o abcpp/dvtraitsim_cpp`python3-config --extension-suffix`

#!/bin/bash
# Change accorsingly. Tested for python 2.7 as well.
PYTHON=python3.5
#PYTHON=python2.7

# Make an object file. "../" is the path to the EIGEN library
g++  -std=c++14 -march=native -mfpmath=sse -ffast-math -O3 -c -fPIC -I ../ sumc.cpp

# Generate a C++ interface source code. Needs SWIG to be installed
swig -builtin -c++ -python sumc.i

# Compile the interface source file and link everything
g++ -std=c++14 -march=native -mfpmath=sse -ffast-math -O3 -I ../ -I/usr/include/${PYTHON} -I/usr/lib/${PYTHON} -c -fPIC sumc_wrap.cxx
g++ -shared -Wl,-soname,_sumc.so -o _sumc.so sumc.o sumc_wrap.o

# Comple a c++ version too.
g++ -std=c++14 -march=native -mfpmath=sse -ffast-math -O3 -I ../ sumc.cpp test.cpp -o test

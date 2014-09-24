#!/bin/bash

# steps to build the ppi_networkit module
git submodule init
git submodule update

# build NetworKit and the Subgraphs library
echo "#######################################################################"
echo "#              Building NetworKit and Subgraphs library               #"
echo "#######################################################################"

cd ppi_networkit
mkdir build
cd build
cmake ../
make
cd ..

# build the cython module `ppi_networkit`
echo "#######################################################################"
echo "#                      Building cython extension                      #"
echo "#######################################################################"

cd cython
python3 setup.py build_ext --inplace


cd ..
cd ..

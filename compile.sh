#!/bin/bash

# user parameter
TARGET="PlasmaScale"

# cmake parameters
export CC=gcc
export CXX=g++

# build directory
cd /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/build

# automatic makefile generation and compilation
echo "====================================================================================="
echo "                                      CMAKE                                          "
echo "====================================================================================="
cmake ..

echo "====================================================================================="
echo "                                       MAKE                                          " 
echo "====================================================================================="
make ${TARGET}
# make


# additional run
# echo "====================================================================================="
# ./app/${TARGET} -i /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/inp/2stream.inp
# echo "====================================================================================="

#!/bin/bash

# user parameter
TARGET="PlasmaScale"

# cmake parameters
export CC=clang
export CXX=clang++

# build directory
cd /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/build2

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
echo "====================================================================================="
#rm /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/test_gradientinfo.dmp
./app/${TARGET} -i /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/ionwave.inp #-d /Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/test_gradientinfo.dmp -dp 1
echo "====================================================================================="

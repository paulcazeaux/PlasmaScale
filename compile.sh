#!/bin/bash

# user parameter
TARGET="PlasmaScale"
#TARGET="t-HaarTools"

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
InputFile=/Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/app/cfg/ionwave.inp
ExportFile=/Users/cazeaux/Dropbox/Postdocs/Plasma/Output/Binning/t-test.dmp
rm ${ExportFile}
./app/${TARGET} -i ${InputFile}  -d ${ExportFile} -dp 1
echo "====================================================================================="

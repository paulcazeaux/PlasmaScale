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
InputFile=/Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/ionwave.inp
ExportFile=/Users/cazeaux/Dropbox/Postdocs/Plasma/PlasmaScale/TestIonWave/EFPIt3.dmp
rm ${ExportFile}
./app/${TARGET} -i ${InputFile}  -d ${ExportFile} -dp 1
echo "====================================================================================="

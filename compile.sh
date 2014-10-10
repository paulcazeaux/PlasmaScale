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
for i in {170..170}
do
	ExportFile=~/Dropbox/Postdocs/Plasma/Output/TestExpansion/t-$i.dmp
	rm ${ExportFile}
	./app/${TARGET} -i ${InputFile} -d ${ExportFile} -dp 1 -nox -s 10 -nm $i
done
echo "====================================================================================="

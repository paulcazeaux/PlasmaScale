#!/bin/bash

# user parameter
TARGET="PlasmaScale"
#TARGET="t-HaarTools"

# cmake parameters
export CC=clang
export CXX=clang++

# build directory
cd /Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/build

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
echo "                                     EXECUTION                                       "
echo "====================================================================================="

InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp

ExportFile=/Users/cazeaux/Desktop/Output/PlasmaExpansion/Size_2000/null
rm ${ExportFile}.dmp ${ExportFile}.out
./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 30 -nox > ${ExportFile}.out
open ${ExportFile}.out

# ExportFolder=/Users/cazeaux/Desktop/Output/PlasmaExpansion/Test_ndt_Size_500/

# ndt=1
# while [  $ndt -lt 45  ]; do
# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_coframe_512
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 25 -s 100 -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_coframe_512
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 25 -s 100 -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_coframe_512
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 25 -s 100  -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_coframe_512
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 25 -s 100  -nox > ${ExportFile}.out

# 	let ndt=ndt+1
# 	sleep 5
# 	#open ${ExportFile}.out
# done
# wait
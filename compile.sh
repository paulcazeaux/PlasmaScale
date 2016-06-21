#!/bin/bash

# user parameter
TARGET="PlasmaScale"
#TARGET="t-HaarTools"
#TARGET="t-Eigen"

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


# additional run
echo "====================================================================================="
echo "                                     EXECUTION                                       "
echo "====================================================================================="

InputFile=~/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp

# ExportFile=~/Gaussian_20/Reference
# rm ${ExportFile}.dmp ${ExportFile}.out
# ./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 30 -nox > ${ExportFile}.out
# open ${ExportFile}.out

InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp

ExportFile=~/Gaussian_20/EFPI_x2048_250
rm ${ExportFile}.dmp ${ExportFile}.out
./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 30 -nox > ${ExportFile}.out
open ${ExportFile}.out

# ExportFolder=/Users/cazeaux/Desktop/Output/PlasmaExpansion/Gaussian/

# ndt=1
# while [  $ndt -lt 45  ]; do
# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_test_ndt
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 10 -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_test_ndt
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 10 -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_test_ndt
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 10 -nox > ${ExportFile}.out &

# 	let ndt=ndt+1
# 	sleep 5

# 	InputFile=/Users/cazeaux/Dropbox/Workplace/Archive/EPFL/Plasma/PlasmaScale/app/cfg/plasmaexpansion.inp
# 	echo ":6 s/\(^\s*[-+]\=\d\+[.]\=[-+eE0-9]*\s\+[-+]\=\d\+[.]\=[-+eE0-9]*\s\+\)\(\d\+\)/\1${ndt}/g
# 	w
# 	q" | ex ${InputFile}

# 	ExportFile=${ExportFolder}${ndt}_test_ndt
# 	rm ${ExportFile}.dmp ${ExportFile}.out
# 	./app/${TARGET} -i ${InputFile} -d ${ExportFile}.dmp -dp 1 -s 10 -nox > ${ExportFile}.out

# 	let ndt=ndt+1
# 	sleep 5
# 	#open ${ExportFile}.out
# done
# wait
cmake_minimum_required(VERSION 2.6) # Pour éviter un warning au lancement de CMake, on doit préciser la version minimale

######################################################################################
### 0) Set compiler
######################################################################################

#Set(CMAKE_C_COMPILER mpicc)
#Set(CMAKE_CXX_COMPILER mpicxx)

######################################################################################

######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/LA)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/JSON)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/MOPT)

######################################################################################


######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

######################################################################################
# Boost inclusion
SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_MULTITHREAD OFF)
INCLUDE(FindBoost)
FIND_PACKAGE( Boost 1.47.0 COMPONENTS system filesystem REQUIRED )
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIR})

######################################################################################
### 3) Define your target(s): just the flow-shop lib here
######################################################################################

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

FILE( GLOB_RECURSE source_files *.cpp )

ADD_LIBRARY(MOPT STATIC ${source_files})

cmake_minimum_required(VERSION 2.6) # Pour éviter un warning au lancement de CMake, on doit préciser la version minimale
######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)

######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)


######################################################################################
### 3) Define your target: just an executable here
######################################################################################

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/app)

ADD_EXECUTABLE(PlasmaScale  PlasmaScale.cpp)
ADD_DEPENDENCIES(PlasmaScale PSCALE)


######################################################################################
### 4) Link the librairies for your executable
######################################################################################
TARGET_LINK_LIBRARIES( PlasmaScale PSCALE ${XGRAFIX_LIBRARIES} ${FFTW_LIBRARIES})

######################################################################################
cmake_minimum_required(VERSION 2.6) # Pour éviter un warning au lancement de CMake, on doit préciser la version minimale
######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)
######################################################################################
### 2) Define your target(s)
######################################################################################

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/app)

FILE( GLOB_RECURSE sourcefiles *.cpp )

ADD_EXECUTABLE(PlasmaScale  ${sourcefiles})


######################################################################################
### 3) Link the librairies for your executable
######################################################################################
TARGET_LINK_LIBRARIES( PlasmaScale ${XGRAFIX_LIBRARIES} ${FFTW_LIBRARIES})

######################################################################################
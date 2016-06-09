cmake_minimum_required(VERSION 2.6) # Pour éviter un warning au lancement de CMake, on doit préciser la version minimale
######################################################################################
### 1) Include the sources
######################################################################################

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)

######################################################################################
### 3) Define your target(s): just the flow-shop lib here
######################################################################################

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

FILE( GLOB_RECURSE source_files *.cpp )

ADD_LIBRARY(PSCALE STATIC ${source_files})


######################################################################################
### 3) Link the librairies for your library
######################################################################################
TARGET_LINK_LIBRARIES( PSCALE ${XGRAFIX_LIBRARIES} ${FFTW_LIBRARIES} ${UMFPACK_LIBRARIES})

######################################################################################
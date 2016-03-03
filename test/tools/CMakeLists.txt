###############################################################################
##
## CMakeLists file for representations/test
##
###############################################################################


######################################################################################
### 1) Include the sources
######################################################################################

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/tools)

######################################################################################


######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
######################################################################################



######################################################################################
### 3) Define your targets and link the librairies
######################################################################################

SET (TEST_LIST
	t-HaarTools
	#t-Tools
)

FOREACH (test ${TEST_LIST})
	  SET ("T_${test}_SOURCES" "${test}.cpp")
ENDFOREACH (test)

#CREATE EXECUTABLES AND ADD ITS AS TEST 
FOREACH (test ${TEST_LIST})
	  ADD_EXECUTABLE(${test} ${T_${test}_SOURCES}) 
	  ADD_TEST(${test} ${test}) 
ENDFOREACH (test)
	
# Link the librairies
FOREACH (test ${TEST_LIST})
 TARGET_LINK_LIBRARIES(${test} PSCALE)
ENDFOREACH (test)
	
######################################################################################


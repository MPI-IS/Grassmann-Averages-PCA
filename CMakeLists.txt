cmake_minimum_required(VERSION 2.8)
project(RobustPCA)

include(CTest)



# thirdparties
# testing using boost
set(Boost_USE_MULTITHREADED ON)

set(Boost_DEBUG ON)
set(Boost_DETAILED_FAILURE_MSG ON)
if(NOT DEFINED BOOST_ROOT)
  set(Boost_NO_SYSTEM_PATHS ON)
else()
  set(Boost_NO_SYSTEM_PATHS OFF)
endif()
set(Boost_NO_BOOST_CMAKE ON)



find_package(Boost COMPONENTS unit_test_framework)
if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()



include_directories(${Boost_INCLUDE_DIRS})
link_directories (${Boost_LIBRARY_DIRS})
file(TO_NATIVE_PATH ${Boost_LIBRARY_DIRS} Boost_LIBRARY_DIRS_NATIVE)

# avoiding boost.test autolink
add_definitions(-DBOOST_TEST_NO_LIB)

if(WIN32)
  set(TEST_DYNAMIC_LIBRARY_PATH_CMD "PATH")
elseif(APPLE)
  set(TEST_DYNAMIC_LIBRARY_PATH_CMD "DYLD_LIBRARY_PATH")
else()
  set(TEST_DYNAMIC_LIBRARY_PATH_CMD "LD_LIBRARY_PATH")
endif()





include_directories(${RobustPCA_SOURCE_DIR})
set(RobustPCA_LIB
     include/private/boost_ublas_matlab_helper.hpp
     include/private/boost_ublas_matrix_helper.hpp

     include/robust_pca.hpp
     src/robust_pca.cpp)
add_library(robust_pca ${RobustPCA_LIB})









set(RobustPCA_tests
    test/test_main.hpp
    test/test_main.cpp

    test/test_robust_pca.cpp
    test/test_robust_pca_trimmed.cpp
    test/test_row_proxy.cpp)
add_executable(robustpca_test ${RobustPCA_tests})
target_link_libraries(robustpca_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# the main definition is in there
set_source_files_properties(test/test_main.cpp PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_MAIN")

# the link is with dynamic 
set_target_properties(robustpca_test PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")


add_test(NAME robustpca_test-1
         COMMAND robustpca_test)
# basically, we cannot concatenate PATH with something else. I tried different configuration, and the spaces and ; in PATH do
# not play well with the command itself. Even the PATH=%PATH%;blablabla with space does not seem to work.
set_tests_properties(
          robustpca_test-1 
          PROPERTIES 
            ENVIRONMENT 
            ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${Boost_LIBRARY_DIRS})
            #"${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${TEST_DYNAMIC_LIBRARY_PATH_CMD_PREVIOUS}${TEST_DYNAMIC_LIBRARY_PATH_CMD_SEP}${Boost_LIBRARY_DIRS_NATIVE}")
            #"${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${TEST_DYNAMIC_LIBRARY_PATH_CMD_INV_PREFIX}${TEST_DYNAMIC_LIBRARY_PATH_CMD}${TEST_DYNAMIC_LIBRARY_PATH_CMD_INV_SUFFIX}${TEST_DYNAMIC_LIBRARY_PATH_CMD_SEP}${Boost_LIBRARY_DIRS_NATIVE}")



# currently under maintainance for the cmake community
add_custom_target(MatlabCMakePackage
  SOURCES 
  ${CMAKE_SOURCE_DIR}/extensions/FindMatlab.cmake
  )

# simple wart to the current state of the Matlab release into the community
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/extensions)
find_package(Matlab)


  
message("MATLAB_INCLUDE_DIR = ${MATLAB_INCLUDE_DIR}")
if(NOT MATLAB_INCLUDE_DIR)
  message(FATAL_ERROR "Unable to find Matlab include directory")
endif()

include_directories(${MATLAB_INCLUDE_DIR})
set(RobustPCAMatlab_LIB
  extensions/matlab.cpp)
add_library(robust_pca_m SHARED ${RobustPCAMatlab_LIB})
target_link_libraries(robust_pca_m ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY})

set_target_properties(robust_pca_m PROPERTIES PREFIX "")
set_target_properties(robust_pca_m PROPERTIES SUFFIX ".${MEX_EXTENSION}")
if(WIN32)
  set_target_properties(robust_pca_m PROPERTIES DEFINE_SYMBOL "DLL_EXPORT_SYM=__declspec(dllexport)")
endif()



# unit tests in matlab
add_custom_target(MatlabUnitTests
  SOURCES 
  test/matlab_unit_tests.m
  )


set(MATLAB_UNIT_TESTS_CMD -nosplash -nojvm -nodesktop -nodisplay)
if(WIN32)
  set(MATLAB_UNIT_TESTS_CMD ${MATLAB_UNIT_TESTS_CMD} -wait)
endif()
#        COMMAND ${MATLAB_PROGRAM}  -logfile ${CMAKE_BINARY_DIR}/matlabVersionLog.cmaketmp  -r "version, exit" 
#        OUTPUT_VARIABLE MATLAB_VERSION_FROM_CMD_dummy
#        RESULT_VARIABLE resultMatlabVersionCall
#        TIMEOUT 30


add_test(NAME robustpca_matlabtest-1
         COMMAND ${MATLAB_PROGRAM} ${MATLAB_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/robustpca_matlabtest_1.log 
         -r "addpath('$<TARGET_FILE_DIR:robust_pca_m>', '${CMAKE_SOURCE_DIR}/test');runtests('matlab_unit_tests'), exit")
set_tests_properties(robustpca_matlabtest-1 PROPERTIES TIMEOUT 30)
set_tests_properties(robustpca_matlabtest-1 PROPERTIES FAIL_REGULAR_EXPRESSION "[1-9][0-9]*[\r\n\t ]*[fF]ailed")



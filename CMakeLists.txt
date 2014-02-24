cmake_minimum_required(VERSION 2.8)
project(RobustPCA)

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_DETAILED_FAILURE_MSG ON)

set(Boost_USE_MULTITHREADED ON)
set(Boost_DEBUG             TRUE)



include_directories(${RobustPCA_SOURCE_DIR})
set(RobustPCA_LIB
     include/robust_pca.hpp
     src/robust_pca.cpp)
     
     
set(RobustPCA_Matlab
     extension/matlab/robust_pca.mex)
     

add_library(robust_pca ${RobustPCA_LIB})

set(BOOST_ROOT  "D:/Code/thirdparties/visual2012x64/boost_1_55" CACHE PATH "boost root directory")
find_package(boost 1.55.0 COMPONENTS test_exec_monitor unit_test_framework)
if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "currently boost is needed")
endif()
message("Boost_TEST_EXEC_MONITOR_LIBRARY = ${Boost_TEST_EXEC_MONITOR_LIBRARY}")
message("Boost_UNIT_TEST_FRAMEWORK_LIBRARY = ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")

include_directories(${Boost_INCLUDE_DIRS})
link_directories (${Boost_LIBRARY_DIRS})

set(RobustPCA_tests
    test/test_client_function.cpp)

add_executable(robustpca_test ${RobustPCA_tests})
target_link_libraries(robustpca_test ${Boost_TEST_EXEC_MONITOR_LIBRARY})




#find_package(Matlab)
find_library(MATLAB_MEX_LIBRARY
  libmex
  ${MATLAB_ROOT}
)
find_library(MATLAB_MX_LIBRARY
  libmx
  ${MATLAB_ROOT}
)
find_library(MATLAB_ENG_LIBRARY
  libeng
  ${MATLAB_ROOT}
)

find_path(MATLAB_INCLUDE_DIR
  "mex.h"
  ${MATLAB_ROOT})


if(NOT ${MATLAB_INCLUDE_DIR})
  message(FATAL_ERROR "Unable to find Matlab") # todo make it a 
endif()

include_directories(${MATLAB_INCLUDE_DIR})
set(RobustPCAMatlab_LIB
  extensions/matlab.cpp)
add_library(robust_pca_m SHARED ${RobustPCAMatlab_LIB})
target_link_libraries(robust_pca_m ${MATLAB_MEX_LIBRARY})

# Raffi/ I think this one may be done automatically
set_target_properties(robust_pca_m PROPERTIES SUFFIX ".mexw64")
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

set(BOOST_ROOT  "D:/Code/externals" CACHE PATH "boost root directory")
find_package(boost COMPONENTS test_exec_monitor)
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
# Copyright 2014, Max Planck Institute for Intelligent Systems.
# Distributed under the BSD 3-Clause license.
# (See accompanying file LICENSE.txt or copy at
# http://opensource.org/licenses/BSD-3-Clause)

cmake_minimum_required(VERSION 2.8)
project(RobustPCA)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)


if(NOT WITHOUT_TESTS)
  # ctest sets BUILD_TESTING automatically, but does not seem to serve its purpose.
  include(CTest)
  enable_testing()
endif()



# thirdparties

# Boost, uBlas and some other header only libraries are needed in any case. The boost.test framework
# is optional.

# disable auto link
add_definitions(-DBOOST_ALL_NO_LIB)
if(NOT WITHOUT_TESTS)
  set(Boost_USE_MULTITHREADED ON)

  set(Boost_DEBUG ON)
  set(Boost_DETAILED_FAILURE_MSG ON)
  if(NOT DEFINED BOOST_ROOT)
    set(Boost_NO_SYSTEM_PATHS ON)
  else()
    set(Boost_NO_SYSTEM_PATHS OFF)
  endif()
  set(Boost_NO_BOOST_CMAKE ON)

  find_package(Boost COMPONENTS unit_test_framework system thread)
else()
  find_package(Boost COMPONENTS system thread)
endif()

if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()



# Matlab bindings
if(NOT WITHOUT_MATLAB)
  # simple wart to the current state of the Matlab release into the community
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/extensions)
  set(MATLAB_FIND_DEBUG TRUE)
  find_package(MATLAB 8.2 EXACT REQUIRED COMPONENTS MX_LIBRARY MAIN_PROGRAM)
  
  # the following lines should not be needed
  if(NOT ${MATLAB_FOUND})
    message(FATAL_ERROR "Unable to find Matlab include directory")
  else()
    message(STATUS "[MATLAB] - version ${MATLAB_VERSION}")
    message(STATUS "[MATLAB] - include directory ${MATLAB_INCLUDE_DIR}")
    message(STATUS "[MATLAB] - mex library ${MATLAB_MEX_LIBRARY}")
    message(STATUS "[MATLAB] - mex extension ${MATLAB_MEX_EXTENSION}")
    message(STATUS "[MATLAB] - mx library ${MATLAB_MX_LIBRARY}")
    message(STATUS "[MATLAB] - matlab ${MATLAB_PROGRAM}")
  endif()

  # Adding a project for the script that finds Matlab, the script being under development currently
  add_custom_target(MatlabCMakePackage
    SOURCES 
    ${CMAKE_SOURCE_DIR}/extensions/FindMatlab.cmake
    )
  set_target_properties(MatlabCMakePackage PROPERTIES FOLDER "Community")
endif()







include_directories(${Boost_INCLUDE_DIRS})



# main library
include_directories(${RobustPCA_SOURCE_DIR})
set(RobustPCA_LIB
     include/private/boost_ublas_matlab_helper.hpp
     include/private/boost_ublas_matrix_helper.hpp

     include/robust_pca.hpp
     src/robust_pca.cpp)
add_library(robustpca ${RobustPCA_LIB})
target_link_libraries(robustpca ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})







# Unit tests using boost
if(NOT WITHOUT_TESTS)
  # avoiding boost.test autolink
  add_definitions(-DBOOST_TEST_NO_LIB)

  if(WIN32)
    set(TEST_DYNAMIC_LIBRARY_PATH_CMD "PATH")
  elseif(APPLE)
    set(TEST_DYNAMIC_LIBRARY_PATH_CMD "DYLD_LIBRARY_PATH")
  else()
    set(TEST_DYNAMIC_LIBRARY_PATH_CMD "LD_LIBRARY_PATH")
  endif()



  set(RobustPCA_tests
      test/test_main.hpp
      test/test_main.cpp

      test/test_robust_pca.cpp
      test/test_robust_pca_trimmed.cpp
      test/test_row_proxy.cpp)
  add_executable(robustpca_test ${RobustPCA_tests})
  target_link_libraries(robustpca_test robustpca ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  # the main definition is in there
  set_source_files_properties(test/test_main.cpp PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_MAIN")

  # the link is with dynamic 
  set_target_properties(robustpca_test PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
  set_target_properties(robustpca_test PROPERTIES FOLDER "UnitTests")


  add_test(
    NAME robustpca_test-1
    COMMAND robustpca_test)

  # basically, we cannot concatenate PATH with something else. I tried different configuration, and the spaces and ; in PATH do
  # not play well with the command itself. Even the PATH=%PATH%;blablabla with space does not seem to work.
  set_tests_properties(
    robustpca_test-1 
    PROPERTIES ENVIRONMENT ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${Boost_LIBRARY_DIRS})
endif()





# adding the matlab MEX extensions
if(NOT WITHOUT_MATLAB)
  

  include_directories(${MATLAB_INCLUDE_DIR})
  set(RobustPCAMEX_LIB
      extensions/matlab.cpp)
  add_library(robustpca_m SHARED ${RobustPCAMEX_LIB})
  target_link_libraries(robustpca_m robustpca ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY})

  set_target_properties(robustpca_m PROPERTIES PREFIX "")
  set_target_properties(robustpca_m PROPERTIES SUFFIX ".${MEX_EXTENSION}")
  if(WIN32)
    set_target_properties(robustpca_m PROPERTIES DEFINE_SYMBOL "DLL_EXPORT_SYM=__declspec(dllexport)")
  endif()
  
  set_target_properties(robustpca_m PROPERTIES FOLDER "Matlab")
endif()



# adding the matlab unit tests
if(NOT WITHOUT_MATLAB AND NOT WITHOUT_TESTS)
  add_custom_target(
    robustpca_matlab_unittest_scripts 
    SOURCES test/robustpca_matlab_unit_tests.m)
  set_target_properties(robustpca_matlab_unittest_scripts PROPERTIES FOLDER "UnitTests")

  set(MATLAB_UNIT_TESTS_CMD -nosplash -nojvm -nodesktop -nodisplay)
  if(WIN32)
    set(MATLAB_UNIT_TESTS_CMD ${MATLAB_UNIT_TESTS_CMD} -wait)
  endif()

  macro(copy_dependency_with_config target_name dependency_name)
    set(dependency_name_debug ${${dependency_name}_DEBUG})
    set(dependency_name_non_debug ${${dependency_name}_RELEASE})
    add_custom_command(
      TARGET ${target_name}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo copy ${dependency_name} (
	      $<$<CONFIG:Debug>:${dependency_name_debug}>
	      $<$<NOT:$<CONFIG:Debug>>:${dependency_name_non_debug}>) 
	      into 
	      $<TARGET_FILE_DIR:${target_name}>/.

      COMMAND ${CMAKE_COMMAND} -E copy_if_different 
	      $<$<CONFIG:Debug>:${dependency_name_debug}>
	      $<$<NOT:$<CONFIG:Debug>>:${dependency_name_non_debug}> 
	      $<TARGET_FILE_DIR:${target_name}>/.
      COMMENT "Copy ${target_name} dependencies into the output folder"
    )    
  


  endmacro(copy_dependency_with_config)
  
  copy_dependency_with_config(robustpca_m Boost_SYSTEM_LIBRARY)
  copy_dependency_with_config(robustpca_m Boost_THREAD_LIBRARY)




  add_test(NAME robustpca_matlabtest-1
            COMMAND ${MATLAB_PROGRAM} ${MATLAB_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/robustpca_matlabtest_1.log 
           -r "addpath('$<TARGET_FILE_DIR:robustpca_m>', '${CMAKE_SOURCE_DIR}/test', '${Boost_LIBRARY_DIRS}'); path, runtests('robustpca_matlab_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(robustpca_matlabtest-1 PROPERTIES TIMEOUT 30)
  # check to see if the returned error works on windows. The line below does definitely not work on windows, a driver should
  # check the content of the log file.
  #set_tests_properties(robustpca_matlabtest-1 PROPERTIES FAIL_REGULAR_EXPRESSION "[1-9][0-9]*[\r\n\t ]*[fF]ailed")
endif()


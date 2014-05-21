# Copyright 2014, Max Planck Institute for Intelligent Systems.
# Distributed under the BSD 3-Clause license.
# (See accompanying file LICENSE.txt or copy at
# http://opensource.org/licenses/BSD-3-Clause)

cmake_minimum_required(VERSION 2.8.11)
project(RobustPCA)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)


if(NOT WITHOUT_TESTS)
  # ctest sets BUILD_TESTING automatically, but does not seem to serve its purpose.
  include(CTest)
  enable_testing()
endif()

#set(CMAKE_MACOSX_RPATH 1)

# thirdparties

# Boost, uBlas and some other header only libraries are needed in any case. The boost.test framework
# is optional.

# disable auto link
add_definitions(-DBOOST_ALL_NO_LIB)
set(Boost_REALPATH ON)
if(NOT WITHOUT_TESTS)
  set(Boost_USE_MULTITHREADED ON)

  set(Boost_DEBUG ON)
  set(Boost_DETAILED_FAILURE_MSG ON)
  if(DEFINED BOOST_ROOT)
    set(Boost_NO_SYSTEM_PATHS ON)
  else()
    set(Boost_NO_SYSTEM_PATHS OFF)
  endif()
  set(Boost_NO_BOOST_CMAKE ON)

  find_package(Boost COMPONENTS unit_test_framework system thread chrono)
else()
  find_package(Boost COMPONENTS system thread chrono)
endif()

if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()



# Matlab bindings
if(NOT WITHOUT_MATLAB)
  # simple wart to the current state of the Matlab release into the community
  set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/extensions ${CMAKE_MODULE_PATH})
  set(MATLAB_FIND_DEBUG TRUE)
  find_package(Matlab 8.2 REQUIRED COMPONENTS MX_LIBRARY MAIN_PROGRAM)
  
  # the following lines should not be needed
  if(NOT ${Matlab_FOUND})
    message(FATAL_ERROR "Unable to find Matlab include directory")
  else()
    message(STATUS "[MATLAB] - version ${Matlab_VERSION_STRING}")
    message(STATUS "[MATLAB] - include directory ${Matlab_INCLUDE_DIRS}")
    message(STATUS "[MATLAB] - mex library ${Matlab_MEX_LIBRARY}")
    message(STATUS "[MATLAB] - mex extension ${Matlab_MEX_EXTENSION}")
    message(STATUS "[MATLAB] - mx library ${Matlab_MX_LIBRARY}")
    message(STATUS "[MATLAB] - matlab ${Matlab_PROGRAM}")
  endif()

  # Adding a project for the script that finds Matlab, the script being under development currently
  add_custom_target(MatlabCMakePackage
    SOURCES 
    ${CMAKE_SOURCE_DIR}/extensions/FindMatlab.cmake
    )
  set_target_properties(MatlabCMakePackage PROPERTIES FOLDER "Community")
endif()


# general defines
if(WIN32)
  add_definitions(-D_WIN32_WINNT=0x0501)
else()
  # this is mainly because of Boost to which the mex file is linked. 
  # Basically, it will change the rpath of the produced .mex file to point to either $ORIGIN (ldd variants)
  # or @loader_path (otool variants). Then dependant .so/.dylib will be found relatively to the .mex file, which
  # is in this case the current directory. 
  # This means also that the boost dependencies should be copied into the target directory of the .mex file, which is 
  # performed below in this file. However, this is not enough and in most case a wrapper should be implemented in order
  # to load the .mex with the necessary load flags (load flags that honors the .mex symbols from its dependencies prior to the ones of
  # the process).
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  if(NOT APPLE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
  elseif(UNIX) # APPLE is Unix
    message("tototototot loader path")
    
    #set(CMAKE_INSTALL_RPATH "@loader_path") # to be tested on osx
    set(CMAKE_INSTALL_NAME_DIR "@loader_path") # to be tested on osx
  endif()
  #set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()



# compilation options
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++11 HAS_CXX11_FLAG)
check_cxx_compiler_flag(-std=c++0x HAS_CXX0X_FLAG)
if (HAS_CXX11_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(HAS_CXX0X_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()




include_directories(${Boost_INCLUDE_DIRS})



# main library
include_directories(${RobustPCA_SOURCE_DIR})
set(RobustPCA_LIB
     include/robust_pca.hpp
     include/robust_pca_trimming.hpp
     include/robust_pca_stable_trimming.hpp
     include/private/utilities.hpp

     include/private/boost_ublas_matlab_helper.hpp
     include/private/boost_ublas_matrix_helper.hpp

     src/robust_pca.cpp)
add_library(robustpca ${RobustPCA_LIB})
target_link_libraries(robustpca ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})
set_target_properties(robustpca PROPERTIES INSTALL_RPATH "@loader_path/.")







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
      test/test_robust_pca_trimming.cpp
      test/test_row_proxy.cpp
      test/test_k_first.cpp)
  add_executable(robustpca_test ${RobustPCA_tests})
  target_link_libraries(robustpca_test robustpca ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_CHRONO_LIBRARY})

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
  list(LENGTH Boost_LIBRARY_DIRS _boost_lib_length)
  if(${_boost_lib_length} GREATER 1)
    list(GET Boost_LIBRARY_DIRS 0 _boost_lib_env)
  else()
    list(GET Boost_LIBRARY_DIRS 0 _boost_lib_env)
  endif()
  set_tests_properties(
    robustpca_test-1 
    PROPERTIES ENVIRONMENT ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${_boost_lib_env})
endif()





# adding the matlab MEX extensions
if(NOT WITHOUT_MATLAB)
  

  include_directories(${Matlab_INCLUDE_DIRS})
  set(RobustPCAMEX_LIB
      extensions/matlab.cpp)
  add_library(robustpca_m SHARED ${RobustPCAMEX_LIB})
  target_link_libraries(robustpca_m robustpca ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})

  set_target_properties(robustpca_m PROPERTIES PREFIX "")
  #set_target_properties(robustpca_m PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  set(Matlab_NEED_WRAPPER TRUE)

  # just for convenience, does not produce any artifacts.
  add_custom_target(
    matlab_mex_dynamic_loader
    SOURCES ${CMAKE_SOURCE_DIR}/extensions/matlab_wrapper.c.template)
  set_target_properties(matlab_mex_dynamic_loader PROPERTIES FOLDER "Matlab")



  macro(create_matlab_mex_wrapper library_name)

    set(_matlab_wrapper_template ${CMAKE_SOURCE_DIR}/extensions/matlab_wrapper.c.template)
    set(_matlab_wrapper_template_output ${CMAKE_BINARY_DIR}/Matlab/matlab_wrapper_${library_name}.c)
    set_target_properties(${library_name} 
      PROPERTIES 
        PREFIX ""
        SUFFIX ".${Matlab_MEX_EXTENSION}_wrapped") # to avoid the possibility to call it from matlab*
    
    # the template contains all the necessary code. The library to load is defined by a #define
    configure_file(${_matlab_wrapper_template} ${_matlab_wrapper_template_output} COPYONLY)
    

    add_library(${library_name}_matlab_loader SHARED ${_matlab_wrapper_template_output})
    target_compile_definitions(${library_name}_matlab_loader PUBLIC "MATLAB_SO_WRAP=\"./$<TARGET_FILE_NAME:${library_name}>\"") 
    target_link_libraries(${library_name}_matlab_loader dl)
    
    # this order is to be kept: we construct the loader before the real mex, since there is not link dependency here, it should work.
    # this also enables to honor the dependencies of the primary target mex file. 
    add_dependencies(${library_name} ${library_name}_matlab_loader) # do not change to target_link_library here
    set_target_properties(${library_name}_matlab_loader 
      PROPERTIES 
        PREFIX ""
        OUTPUT_NAME ${library_name}
        SUFFIX ".${Matlab_MEX_EXTENSION}"
        FOLDER "Matlab")

    unset(_matlab_wrapper_template)
    unset(_matlab_mex_to_load)

  endmacro(create_matlab_mex_wrapper)


  if(UNIX AND Matlab_NEED_WRAPPER)
    create_matlab_mex_wrapper(robustpca_m)
  else()
    set_target_properties(robustpca_m PROPERTIES SUFFIX ".${Matlab_MEX_EXTENSION}")
  endif()



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
  add_dependencies(robustpca_matlab_unittest_scripts robustpca_m)

  set(Matlab_UNIT_TESTS_CMD -nosplash -nojvm -nodesktop -nodisplay)
  if(WIN32)
    set(Matlab_UNIT_TESTS_CMD ${Matlab_UNIT_TESTS_CMD} -wait)
  endif()

  macro(copy_dependency_with_config target_name dependency_name)

    set(dependency_name_debug ${${dependency_name}_DEBUG})
    set(dependency_name_non_debug ${${dependency_name}_RELEASE})

    # check to see if the function GetPrerequisites can do the job

    if(WIN32)
      # on windows, the dependencies are .lib but we should copy .dll files
      get_filename_component(dependency_name_debug1 ${dependency_name_debug} NAME_WE)
      get_filename_component(dependency_name_debug_dir ${dependency_name_debug} DIRECTORY)
      set(dependency_name_debug ${dependency_name_debug_dir}/${dependency_name_debug1}${CMAKE_SHARED_LIBRARY_SUFFIX})

      get_filename_component(dependency_name_non_debug1 ${dependency_name_non_debug} NAME_WE)
      get_filename_component(dependency_name_non_debug_dir ${dependency_name_non_debug} DIRECTORY)
      set(dependency_name_non_debug ${dependency_name_non_debug_dir}/${dependency_name_non_debug1}${CMAKE_SHARED_LIBRARY_SUFFIX})

    endif()


    add_custom_command(
      TARGET ${target_name}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo copy ${dependency_name} 
        $<$<CONFIG:Debug>:${dependency_name_debug}>
        $<$<NOT:$<CONFIG:Debug>>:${dependency_name_non_debug}>
        into 
        $<TARGET_FILE_DIR:${target_name}>/.

      COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        $<$<CONFIG:Debug>:${dependency_name_debug}>
        $<$<NOT:$<CONFIG:Debug>>:${dependency_name_non_debug}> 
        $<TARGET_FILE_DIR:${target_name}>/.
      COMMENT "Copy ${target_name} dependencies into the output folder"
    )    
  


  endmacro(copy_dependency_with_config)
  
  
  if(NOT Boost_USE_STATIC_LIBS)
    copy_dependency_with_config(robustpca_m Boost_SYSTEM_LIBRARY)
    copy_dependency_with_config(robustpca_m Boost_THREAD_LIBRARY)
    copy_dependency_with_config(robustpca_m Boost_CHRONO_LIBRARY)
  endif()



  add_test(NAME robustpca_matlabtest-1
            COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/robustpca_matlabtest_1.log 
           -r "addpath('$<TARGET_FILE_DIR:robustpca_m>', '${CMAKE_SOURCE_DIR}/test'); path, runtests('robustpca_matlab_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(robustpca_matlabtest-1 PROPERTIES TIMEOUT 180)
  # check to see if the returned error works on windows. The line below does definitely not work on windows, a driver should
  # check the content of the log file.
  #set_tests_properties(robustpca_matlabtest-1 PROPERTIES FAIL_REGULAR_EXPRESSION "[1-9][0-9]*[\r\n\t ]*[fF]ailed")
  
  add_test(NAME robustpca_matlabtest-2
            CONFIGURATIONS Release
            COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/robustpca_matlabtest_2.log 
           -r "addpath('$<TARGET_FILE_DIR:robustpca_m>', '${CMAKE_SOURCE_DIR}/test'); path, runtests('robustpca_matlab_performance_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(robustpca_matlabtest-2 PROPERTIES TIMEOUT 240)

  
endif()


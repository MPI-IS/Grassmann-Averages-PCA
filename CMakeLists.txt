# Copyright 2014, Max Planck Institute for Intelligent Systems.
# Distributed under the BSD 3-Clause license.
# (See accompanying file LICENSE.txt or copy at
# http://opensource.org/licenses/BSD-3-Clause)

cmake_minimum_required(VERSION 2.8.11)
project(GrassmannAveragePCA)


# build type, by default to release (with optimisations) 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
#set(Boost_USE_STATIC_LIBS ON)

if(NOT Boost_USE_STATIC_LIBS)
  # link against dynamic libraries
  add_definitions(-DBOOST_ALL_DYN_LINK)
endif()

set(Boost_REALPATH ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_DEBUG ON)
set(Boost_DETAILED_FAILURE_MSG ON)
if(DEFINED BOOST_ROOT)
  set(Boost_NO_SYSTEM_PATHS ON)
else()
  set(Boost_NO_SYSTEM_PATHS OFF)
endif()
set(Boost_NO_BOOST_CMAKE ON)

if(NOT WITHOUT_TESTS)
  find_package(Boost COMPONENTS unit_test_framework system thread chrono date_time)
else()
  find_package(Boost COMPONENTS system thread chrono date_time)
endif()

if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()



# Matlab bindings
if(NOT WITHOUT_MATLAB)
  # simple wart to the current state of the Matlab release into the community
  set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/extensions ${CMAKE_MODULE_PATH})
  set(MATLAB_FIND_DEBUG TRUE)
  find_package(Matlab 8.1 REQUIRED COMPONENTS MX_LIBRARY MAIN_PROGRAM)
  
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
  add_custom_target(
    MatlabCMakePackage
    SOURCES ${CMAKE_SOURCE_DIR}/extensions/FindMatlab.cmake)
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
include(CheckIncludeFileCXX)

set(HAS_AVX FALSE)
check_include_file_cxx(smmintrin.h HAS_SSE41_INSTRINSICS)

# this variable contains the libraries to which the tests should be linked against. In case              
# we are using the static libraries of boost, its content will be expanded.
set(BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_CHRONO_LIBRARY})


if(NOT MSVC)
  # c++11 options
  check_cxx_compiler_flag(-std=c++11 HAS_CXX11_FLAG)
  check_cxx_compiler_flag(-std=c++0x HAS_CXX0X_FLAG)
  # pthread options
  check_cxx_compiler_flag(-pthread   HAS_MINUS_PTHREAD)
  # visbility options 
  check_cxx_compiler_flag(-fvisibility=hidden         HAS_VISIBILITY_HIDDEN)
  check_cxx_compiler_flag(-fvisibility-inlines-hidden HAS_VISIBILITY_INLINE_HIDDEN)
  if(HAS_CXX11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(HAS_CXX0X_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()


  # if we do not do that, the symbols linked from boost remain weak and then clash with the ones
  # defined in the matlab process. So by default the symbols are hidden.
  # this also means that for shared libraries (like MEX), the entry point should be explicitely
  # declared with default visibility, otherwise Matlab cannot find the entry point.i
  if(HAS_VISIBILITY_INLINE_HIDDEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
  endif()
  if(HAS_VISIBILITY_HIDDEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
  endif()

  if(Boost_USE_STATIC_LIBS)
    # does not work, libstdc++ should be compiled with -fPIC
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++") 
    if(HAS_MINUS_PTHREAD AND NOT APPLE)
      # Apparently, compiling with -pthread generated the proper link flags and some defines at compilation
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    endif()

    # update to include rt library (from chrono)
    if(NOT APPLE)
      set(BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS} rt)
    endif()
  endif()
else()
  check_cxx_compiler_flag(/AVX HAS_AVX)
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  set(MSVC_Additional_flags "/fp:fast /GF /Oy /GT /Ox /Ob2 /Oi /Os")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_Additional_flags}")
  if(HAS_SSE41_INSTRINSICS AND HAS_AVX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX ")
    add_definitions(-DROBUSTPCA_INNER_PROD_INTRINSICS_AVX)
  endif()
endif()



include_directories(${Boost_INCLUDE_DIRS})



# main library
include_directories(${GrassmannAveragePCA_SOURCE_DIR})
set(${PROJECT_NAME}_LIB
     include/robust_pca.hpp
     include/robust_pca_trimming.hpp
     include/private/utilities.hpp

     include/private/boost_ublas_external_storage.hpp
     include/private/boost_ublas_row_iterator.hpp

     src/robust_pca.cpp)
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB})
target_link_libraries(${PROJECT_NAME} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
#set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "@loader_path/.")


# examples, checks the compilation only
set(${PROJECT_NAME}Examples_LIB
    test/example_robustpca.cpp)
add_library(${PROJECT_NAME}_examples ${${PROJECT_NAME}Examples_LIB})
target_link_libraries(${PROJECT_NAME}_examples ${PROJECT_NAME})


# documentation of the main doxygen page
add_custom_target(
  doxygen_files
  SOURCES 
    doc/main.md
    doc/Doxyfile
  )
set_target_properties(doxygen_files PROPERTIES FOLDER "Documentation")



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

  # basically, we cannot concatenate PATH with something else. I tried different configuration, and the spaces and ; in PATH do
  # not play well with the command itself. Even the PATH=%PATH%;blablabla with space does not seem to work.
  list(LENGTH Boost_LIBRARY_DIRS _boost_lib_length)
  if(${_boost_lib_length} GREATER 1)
    list(GET Boost_LIBRARY_DIRS 0 _boost_lib_env)
  else()
    list(GET Boost_LIBRARY_DIRS 0 _boost_lib_env)
  endif()


  set(${PROJECT_NAME}_tests
      test/test_main.hpp
      test/test_main.cpp

      test/test_grassmannpca.cpp
      test/test_grassmannpca_trimming.cpp
      test/test_row_proxy.cpp)

  # this file is used for some internal tests related to ordering
  if(${Boost_MINOR_VERSION} GREATER 48)
    set(${PROJECT_NAME}_tests ${${PROJECT_NAME}_tests} test/test_k_first.cpp)
  endif()

  add_executable(${PROJECT_NAME}_test ${${PROJECT_NAME}_tests})
  target_link_libraries(${PROJECT_NAME}_test ${PROJECT_NAME} ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS})

  # the main definition is in there
  set_source_files_properties(test/test_main.cpp PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_MAIN")

  # the link is with dynamic 
  if(NOT Boost_USE_STATIC_LIBS)
    set_target_properties(${PROJECT_NAME}_test PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
  endif()
  set_target_properties(${PROJECT_NAME}_test PROPERTIES FOLDER "UnitTests")

  add_test(
    NAME ${PROJECT_NAME}_test-1
    COMMAND ${PROJECT_NAME}_test)
  if(NOT Boost_USE_STATIC_LIBS)
    set_tests_properties(
      ${PROJECT_NAME}_test-1 
      PROPERTIES ENVIRONMENT ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${_boost_lib_env})
  endif()

  # second test based on files
  add_executable(${PROJECT_NAME}_test_with_files test/test_robust_pca_trimming_from_file.cpp)
  set_source_files_properties(test/test_robust_pca_trimming_from_file.cpp PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_MAIN")
  target_link_libraries(${PROJECT_NAME}_test_with_files ${PROJECT_NAME} ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS}) 
  if(NOT Boost_USE_STATIC_LIBS)
    set_target_properties(${PROJECT_NAME}_test_with_files PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
  endif()
  set_target_properties(${PROJECT_NAME}_test_with_files PROPERTIES FOLDER "UnitTests")

  add_test(
    NAME ${PROJECT_NAME}_test-2
    CONFIGURATIONS Release
    COMMAND ${PROJECT_NAME}_test_with_files --data ${CMAKE_SOURCE_DIR}/test/mat_test.csv --eigen_vectors ${CMAKE_SOURCE_DIR}/test/mat_test_init_vectors.csv --expected_result ${CMAKE_SOURCE_DIR}/test/mat_test_desired_output.csv)
  if(NOT Boost_USE_STATIC_LIBS)
    set_tests_properties(
      ${PROJECT_NAME}_test-2
      PROPERTIES ENVIRONMENT ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${_boost_lib_env})
  endif()

endif()





# adding the matlab MEX extensions
if(NOT WITHOUT_MATLAB)
  
  # just for convenience, does not produce any artifacts.
  # this target is for a dynamic loader which copes with the multiple symbol clashes that may occur
  # with matlab.
  add_custom_target(
    matlab_mex_dynamic_loader
    SOURCES ${CMAKE_SOURCE_DIR}/extensions/matlab_wrapper.c.template)
  set_target_properties(matlab_mex_dynamic_loader PROPERTIES FOLDER "Matlab") 
  
  
  #
  # the MEX file project
  #
  set(GAPCA_MEX_Project ${PROJECT_NAME}_mex)

  include_directories(${Matlab_INCLUDE_DIRS})
  add_library(${GAPCA_MEX_Project} SHARED extensions/matlab.cpp)
  target_link_libraries(${GAPCA_MEX_Project} ${PROJECT_NAME} ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})
  set_target_properties(${GAPCA_MEX_Project} 
      PROPERTIES 
        PREFIX ""
        OUTPUT_NAME ${PROJECT_NAME}
        SUFFIX ".${Matlab_MEX_EXTENSION}"
        FOLDER "Matlab")
  
  # for debug
  #set_target_properties(${GAPCA_MEX_Project} PROPERTIES LINK_FLAGS "-t")
  
  #set_target_properties(${GAPCA_MEX_Project} PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  if(Boost_USE_STATIC_LIBS)
    set(Matlab_NEED_WRAPPER FALSE)
  else()
    set(Matlab_NEED_WRAPPER TRUE)
  endif()


  # entry point in the mex file
  if(WIN32)
    set_target_properties(${GAPCA_MEX_Project} 
      PROPERTIES 
        DEFINE_SYMBOL "DLL_EXPORT_SYM=__declspec(dllexport)")
  else()
    set_target_properties(${GAPCA_MEX_Project}
      PROPERTIES
        DEFINE_SYMBOL "DLL_EXPORT_SYM=__attribute__ ((visibility (\"default\")))"
    )
  endif()


  #
  # create a matlab loader to prevent some symbol crashes
  #
  macro(create_matlab_mex_wrapper library_name)

    set(_matlab_wrapper_template        ${CMAKE_SOURCE_DIR}/extensions/matlab_wrapper.c.template)
    set(_matlab_wrapper_template_output ${CMAKE_BINARY_DIR}/Matlab/matlab_wrapper_${library_name}.c)

    # get the previous output name
    get_target_property(_matlab_library_output_name ${library_name} OUTPUT_NAME)

    # change the suffix to avoid the possibility of matlab calling it
    set_target_properties(${library_name} 
      PROPERTIES 
        PREFIX ""
        SUFFIX ".${Matlab_MEX_EXTENSION}_wrapped") 
    if(UNIX)
      set_target_properties(${library_name}
        PROPERTIES DEFINE_SYMBOL "DLL_EXPORT_SYM=__attribute__ ((visibility (\"default\")))"
      )
    endif()
    
    
    
    # create the loader 
    # the template contains all the necessary code. The library to load is defined by a #define MATLAB_SO_WRAP
    configure_file(${_matlab_wrapper_template} ${_matlab_wrapper_template_output} COPYONLY)
    
    add_library(${library_name}_matlab_loader SHARED ${_matlab_wrapper_template_output})
    target_compile_definitions(${library_name}_matlab_loader PUBLIC "MATLAB_SO_WRAP=\"./$<TARGET_FILE_NAME:${library_name}>\"") 
    target_link_libraries(${library_name}_matlab_loader dl ${Matlab_MEX_LIBRARY})
    
    
    # this order is to be kept: we construct the loader before the real mex, since there is not link dependency here, it should work.
    # this also enables to honour the dependencies of the primary target mex file. 
    add_dependencies(${library_name} ${library_name}_matlab_loader) # do not change to target_link_library here
    set_target_properties(${library_name}_matlab_loader 
      PROPERTIES 
        PREFIX ""
        OUTPUT_NAME ${_matlab_library_output_name}
        SUFFIX ".${Matlab_MEX_EXTENSION}"
        FOLDER "Matlab")

    unset(_matlab_wrapper_template)
    unset(_matlab_mex_to_load)
    unset(_matlab_library_output_name)

  endmacro(create_matlab_mex_wrapper)



  #
  # copies the dependency files into the target output directory
  #
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
      COMMAND ${CMAKE_COMMAND} -E echo Copy ${dependency_name} 
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
  
  
  
  
  # copying the boost dependencies
  if(NOT Boost_USE_STATIC_LIBS)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_SYSTEM_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_THREAD_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_CHRONO_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_DATE_TIME_LIBRARY)
  endif()

  # copying the matlab documentation
  get_target_property(${GAPCA_MEX_Project}_name ${GAPCA_MEX_Project} OUTPUT_NAME)
  set(${GAPCA_MEX_Project}_help_file ${CMAKE_SOURCE_DIR}/extensions/robustpca_m.m)
  add_custom_command(
    TARGET ${GAPCA_MEX_Project}
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo Copy Matlab wrapper documentation into target directory $<TARGET_FILE_DIR:${GAPCA_MEX_Project}>/.
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${GAPCA_MEX_Project}_help_file} $<TARGET_FILE_DIR:${GAPCA_MEX_Project}>/${${GAPCA_MEX_Project}_name}.m
    COMMENT "Copy ${GAPCA_MEX_Project} help into the output folder"
  )


  # create wrappers if needed
  if(UNIX AND Matlab_NEED_WRAPPER)
    create_matlab_mex_wrapper(${GAPCA_MEX_Project})
  endif()
  
  


endif()



# adding the matlab unit tests
if(NOT WITHOUT_MATLAB AND NOT WITHOUT_TESTS)


  set(GAPCA_MEXUnitTest_Project ${GAPCA_MEX_Project}_matlab_unittest_scripts)
  add_custom_target(
    ${GAPCA_MEXUnitTest_Project} 
    SOURCES 
      test/grassmannpca_matlab_unit_tests.m
      test/grassmannpca_matlab_performance_unit_tests.m)

  set_target_properties(${GAPCA_MEXUnitTest_Project} PROPERTIES FOLDER "UnitTests")
  add_dependencies(${GAPCA_MEXUnitTest_Project} ${GAPCA_MEX_Project})

  set(Matlab_UNIT_TESTS_CMD -nosplash -nojvm -nodesktop -nodisplay)
  if(WIN32)
    set(Matlab_UNIT_TESTS_CMD ${Matlab_UNIT_TESTS_CMD} -wait)
  endif()
    
  add_test(NAME ${PROJECT_NAME}_matlabtest-1
           COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/Matlab/${PROJECT_NAME}_matlabtest_1.log 
           -r "addpath('$<TARGET_FILE_DIR:${GAPCA_MEX_Project}>', '${CMAKE_SOURCE_DIR}/test'); path, runtests('grassmannpca_matlab_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(${PROJECT_NAME}_matlabtest-1 PROPERTIES TIMEOUT 180)
  # check to see if the returned error works on windows. The line below does definitely not work on windows, a driver should
  # check the content of the log file.
  #set_tests_properties(${PROJECT_NAME}_matlabtest-1 PROPERTIES FAIL_REGULAR_EXPRESSION "[1-9][0-9]*[\r\n\t ]*[fF]ailed")
  
  add_test(NAME ${PROJECT_NAME}_matlabtest-2
            CONFIGURATIONS Release
            COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/Matlab/${PROJECT_NAME}_matlabtest_2.log 
           -r "addpath('$<TARGET_FILE_DIR:${GAPCA_MEX_Project}>', '${CMAKE_SOURCE_DIR}/test'); path, runtests('grassmannpca_matlab_performance_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(${PROJECT_NAME}_matlabtest-2 PROPERTIES TIMEOUT 600)

  
endif()


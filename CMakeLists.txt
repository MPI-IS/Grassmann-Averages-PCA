# Copyright 2014, Max Planck Society.
# Distributed under the BSD 3-Clause license.
# (See accompanying file LICENSE.txt or copy at
# http://opensource.org/licenses/BSD-3-Clause)

cmake_minimum_required(VERSION 2.8.11)
project(GrassmannAveragesPCA)


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






##############################################################################################
# Thirdparties and tools
#

find_package(Doxygen)


# Boost, uBlas and some other libraries are needed in any case. 

# disable auto link
add_definitions(-DBOOST_ALL_NO_LIB)
set(Boost_USE_STATIC_LIBS OFF)

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
  find_package(Boost COMPONENTS unit_test_framework system thread chrono date_time program_options)
else()
  find_package(Boost COMPONENTS system thread chrono date_time program_options)
endif()

if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()



# Matlab bindings
if(NOT WITHOUT_MATLAB)
  # simple wart to the current state of the Matlab release into the community
  #set(CMAKE_MODULE_PATH /Users/raffi/Personnel/CMake/Modules ${CMAKE_MODULE_PATH})
  set(MATLAB_FIND_DEBUG TRUE)
  find_package(Matlab 8.1 REQUIRED COMPONENTS MX_LIBRARY MAIN_PROGRAM)
  
  # the following lines should not be needed
  if(NOT ${Matlab_FOUND})
    message(FATAL_ERROR "Unable to find Matlab include directory")
  else()
    set(MATLAB_RELEASE_NAME)
    matlab_get_release_name_from_version(${Matlab_VERSION_STRING} MATLAB_RELEASE_NAME)
    message(STATUS "[MATLAB] - version ${Matlab_VERSION_STRING}")
    message(STATUS "[MATLAB] - release ${MATLAB_RELEASE_NAME}")
    message(STATUS "[MATLAB] - include directory ${Matlab_INCLUDE_DIRS}")
    message(STATUS "[MATLAB] - mex library ${Matlab_MEX_LIBRARY}")
    message(STATUS "[MATLAB] - mex extension ${Matlab_MEX_EXTENSION}")
    message(STATUS "[MATLAB] - mx library ${Matlab_MX_LIBRARY}")
    message(STATUS "[MATLAB] - matlab ${Matlab_PROGRAM}")
  endif()

  # Adding a project for the script that finds Matlab, the script being under development currently
  add_custom_target(
    MatlabCMakePackage
    SOURCES 
      ${CMAKE_SOURCE_DIR}/extensions/FindMatlab.cmake
      ${CMAKE_SOURCE_DIR}/extensions/FindMatlab_TestsRedirect.cmake)
  set_target_properties(MatlabCMakePackage PROPERTIES FOLDER "Community")
endif()


# openCV for videos reading applications
# openCV should be installed somewhere defined on the command line. If this is not the case, an error message is printed and 
# the build is aborted.
set(OpenCV_LIB_COMPONENTS opencv_core)
if(NOT OpenCVRoot)
  
  set(OpenCV_CUDA OFF)
  set(OpenCV_STATIC OFF)

 
  find_package(OpenCV)

  if(OpenCV_INCLUDE_DIRS)
    set(OPENCV_AVAILABLE TRUE CACHE INTERNAL "Setting the availability of opencv")
    message(STATUS "[OPENCV] Open CV found at location ${OpenCV_INCLUDE_DIRS}")
  else()
    set(OPENCV_AVAILABLE FALSE CACHE INTERNAL "Setting the availability of opencv")
    message(WARNING "[OPENCV] OpenCVRoot is not defined. OpenCVRoot should be defined with the option -DOpenCVRoot=<root-to-opencv> 
                     in order to benefit from the video applications")
  endif()

else()
  set(opencv_root ${OpenCVRoot})
  if(EXISTS ${opencv_root}/opencv)
    set(opencv_root ${opencv_root}/opencv)
  endif()

  if(NOT EXISTS ${opencv_root}/include)
    message(FATAL_ERROR "[OPENCV] Cannot find the header directory of open cv. Please ensure you have decompressed the version for windows")
  endif()


  # apparently this is the way cmake works... did not know, the OpenCVConfig.cmake file is enough for the configuration
  set(OpenCV_DIR ${opencv_root}/ CACHE PATH "Location of the OpenCV configuration directory")
  set(OpenCV_SHARED ON)
  set(OpenCV_STATIC OFF)
  set(OpenCV_CUDA OFF)
  #set(BUILD_SHARED_LIBS OFF)
  find_package(OpenCV REQUIRED)

  if(NOT OpenCV_INCLUDE_DIRS)
    message(FATAL_ERROR "[OPENCV] Cannot add the OpenCV include directories")
  endif()

  list(REMOVE_DUPLICATES OpenCV_LIB_DIR)
  list(LENGTH OpenCV_LIB_DIR list_lenght)
  if(${list_lenght} GREATER 1)
    list(GET OpenCV_LIB_DIR 0 OpenCV_LIB_DIR)
  endif()
  set(OpenCV_BIN_DIR ${OpenCV_LIB_DIR}/../bin)
  get_filename_component(OpenCV_BIN_DIR ${OpenCV_BIN_DIR} ABSOLUTE)

  set(OPENCV_AVAILABLE TRUE CACHE INTERNAL "Setting the availability of opencv")


endif()






#
# general defines and compilation options
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
  if(HAS_CXX11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(HAS_CXX0X_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()


  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")


  if(Boost_USE_STATIC_LIBS)
    # update to include rt library (from chrono)
    if(NOT APPLE)
      set(BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS} rt)
    endif()
  endif()
else()
  check_cxx_compiler_flag(/AVX HAS_AVX)
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
  set(MSVC_Additional_flags "/fp:fast /GF /Oy /GT /Ox /Ob2 /Oi /Os")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_Additional_flags}")
  if(HAS_SSE41_INSTRINSICS AND HAS_AVX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX ")
    add_definitions(-DGRASSMANNPCA_INNER_PROD_INTRINSICS_AVX)
  endif()
endif()






#
# the project starts here
#

include_directories(${Boost_INCLUDE_DIRS})


##############################################################################################
# Main library
# This library is header only in fact, but we add an empty cpp file for convenience with CMake.
include_directories(${GrassmannAveragesPCA_SOURCE_DIR})
set(${PROJECT_NAME}_LIB
     include/grassmann_pca.hpp
     include/grassmann_pca_with_trimming.hpp
     include/private/utilities.hpp
     include/private/em_pca.hpp
     include/private/chunks_processor.hpp

     include/private/boost_ublas_external_storage.hpp
     include/private/boost_ublas_row_iterator.hpp

     src/grassmann_pca.cpp)
add_library(grassmann_averages ${${PROJECT_NAME}_LIB})
target_link_libraries(grassmann_averages ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
install(
  DIRECTORY include
  DESTINATION "."
  PATTERN "*bak" EXCLUDE)



# examples, checks the compilation only
set(${PROJECT_NAME}Examples_LIB
    test/example_grassmannpca.cpp)
add_library(grassmann_averages_examples ${${PROJECT_NAME}Examples_LIB})
target_link_libraries(grassmann_averages_examples grassmann_averages)

##############################################################################################
# Applications



if(OPENCV_AVAILABLE)

  # not linking against everything in opencv 
  set(opencv_core_library "opencv_core" "opencv_highgui")

  set(${PROJECT_NAME}_ga_movie_runner_SRC
      applications/video_processing_grassmann_pca.cpp)
  add_executable(${PROJECT_NAME}_ga_movie_runner ${${PROJECT_NAME}_ga_movie_runner_SRC})
 
  target_link_libraries(${PROJECT_NAME}_ga_movie_runner grassmann_averages ${opencv_core_library} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(${PROJECT_NAME}_ga_movie_runner PRIVATE ${OpenCV_INCLUDE_DIRS})
  set_target_properties(${PROJECT_NAME}_ga_movie_runner PROPERTIES FOLDER "Applications")



  # Reads a full movie and runs the Trimmed Grassman PCA on it
  set(${PROJECT_NAME}_trimmed_ga_movie_runner_SRC
      applications/video_processing_trimmed_grassmann_pca.cpp)
  add_executable(${PROJECT_NAME}_trimmed_ga_movie_runner ${${PROJECT_NAME}_trimmed_ga_movie_runner_SRC})
 
  target_link_libraries(${PROJECT_NAME}_trimmed_ga_movie_runner grassmann_averages ${opencv_core_library} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(${PROJECT_NAME}_trimmed_ga_movie_runner PRIVATE ${OpenCV_INCLUDE_DIRS})
  set_target_properties(${PROJECT_NAME}_trimmed_ga_movie_runner PROPERTIES FOLDER "Applications")


  # Reads a full movie and runs the EM PCA on it
  set(${PROJECT_NAME}_empca_movie_runner_SRC
      applications/video_processing_empca.cpp)
  add_executable(${PROJECT_NAME}_empca_movie_runner ${${PROJECT_NAME}_empca_movie_runner_SRC})

  target_link_libraries(${PROJECT_NAME}_empca_movie_runner  grassmann_averages ${opencv_core_library} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(${PROJECT_NAME}_empca_movie_runner  PRIVATE ${OpenCV_INCLUDE_DIRS})
  set_target_properties(${PROJECT_NAME}_empca_movie_runner  PROPERTIES FOLDER "Applications")
endif()

##############################################################################################
# Documentation


# documentation of the main doxygen page
# part of the doxygen configuration includes variables that are expanded by cmake.
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/documentation)
configure_file(doc/Doxyfile ${DOXYGEN_OUTPUT_DIRECTORY}/DoxyfileConfigured)
set(doxygen_sources 
    doc/main.md
    doc/Doxyfile
    ${DOXYGEN_OUTPUT_DIRECTORY}/DoxyfileConfigured
  )


if(DOXYGEN_FOUND)
  # generation and "installation" of the doxygen documentation 
  add_custom_target(
    doxygen_documentation
    ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT_DIRECTORY}/DoxyfileConfigured
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    COMMENT "Generating Doxygen documentation"
    SOURCES ${doxygen_sources}
  )
  
  install(
    DIRECTORY ${CMAKE_BINARY_DIR}/documentation/html/
    DESTINATION documentation
    PATTERN "html/*")
  #set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} documentation;documentation/index.html)
else()
  add_custom_target(
    doxygen_documentation
    SOURCES ${doxygen_sources}
  )
endif()
set_target_properties(doxygen_documentation PROPERTIES FOLDER "Documentation")

# the readme file
add_custom_target(
  readme
  SOURCES 
    README.txt
    LICENSE.txt
  )
set_target_properties(readme PROPERTIES FOLDER "Documentation")
install(
  FILES README.txt LICENSE.txt
  DESTINATION ".")
#set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} readme.txt;readme.txt)



##############################################################################################
# Tests

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
      test/test_empca.cpp
      test/test_row_proxy.cpp
      test/test_mus_proximity.cpp
      test/test_faster_acos.cpp
      )

  # this file is used for some internal tests related to ordering
  if(${Boost_MINOR_VERSION} GREATER 48)
    set(${PROJECT_NAME}_tests ${${PROJECT_NAME}_tests} test/test_k_first.cpp)
  endif()

  add_executable(${PROJECT_NAME}_test ${${PROJECT_NAME}_tests})
  target_link_libraries(${PROJECT_NAME}_test grassmann_averages ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS})
  
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
  add_executable(${PROJECT_NAME}_test_with_files test/test_grassmannpca_trimming_from_file.cpp)
  target_link_libraries(${PROJECT_NAME}_test_with_files grassmann_averages ${BOOST_ADDITIONAL_LIBRARIES_FOR_TESTS}) 
  if(NOT Boost_USE_STATIC_LIBS)
    set_target_properties(${PROJECT_NAME}_test_with_files PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
  endif()
  set_target_properties(${PROJECT_NAME}_test_with_files PROPERTIES FOLDER "UnitTests")

  add_test(
    NAME ${PROJECT_NAME}_test-2
    CONFIGURATIONS Release
    COMMAND ${PROJECT_NAME}_test_with_files 
      --data ${CMAKE_SOURCE_DIR}/test/mat_test.csv 
      --basis_vectors ${CMAKE_SOURCE_DIR}/test/mat_test_init_vectors.csv 
      --expected_result ${CMAKE_SOURCE_DIR}/test/mat_test_desired_output.csv)
  if(NOT Boost_USE_STATIC_LIBS)
    set_tests_properties(
      ${PROJECT_NAME}_test-2
      PROPERTIES ENVIRONMENT ${TEST_DYNAMIC_LIBRARY_PATH_CMD}=${_boost_lib_env})
  endif()

endif()



##############################################################################################
# Matlab extensions

# adding the matlab MEX extensions
if(NOT WITHOUT_MATLAB)

  if(WIN32)
    set(MATLAB_INSTALL_DIRECTORY lib/${MATLAB_RELEASE_NAME})
  else()
    set(MATLAB_INSTALL_DIRECTORY lib)
  endif()
  
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
  
  matlab_add_mex(
    NAME ${GAPCA_MEX_Project} 
    REDUCE_VISIBILITY
    OUTPUT_NAME ${PROJECT_NAME}
    SRC extensions/matlab.cpp
    DOCUMENTATION ${${GAPCA_MEX_Project}_help_file}
    LINK_TO grassmann_averages
  )
  set_target_properties(${GAPCA_MEX_Project} PROPERTIES FOLDER "Matlab")  
  
  install(TARGETS ${GAPCA_MEX_Project} 
           DESTINATION ${MATLAB_INSTALL_DIRECTORY})  
  
  install(FILES ${${GAPCA_MEX_Project}_help_file}
           DESTINATION ${MATLAB_INSTALL_DIRECTORY}
           RENAME ${${GAPCA_MEX_Project}_name}.m)  
          


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

      unset(dependency_name_non_debug1)
      unset(dependency_name_non_debug_dir)
      unset(dependency_name_debug1)
      unset(dependency_name_debug_dir)
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
    
    install( FILES ${dependency_name_debug}
              DESTINATION ${MATLAB_INSTALL_DIRECTORY}
              CONFIGURATIONS Debug)
    install( FILES ${dependency_name_non_debug}
              DESTINATION ${MATLAB_INSTALL_DIRECTORY}
              CONFIGURATIONS Release)
              
    unset(dependency_name_debug)
    unset(dependency_name_non_debug)
  
  endmacro(copy_dependency_with_config)
  
  
  
  
  # copying the boost dependencies
  if(NOT Boost_USE_STATIC_LIBS)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_SYSTEM_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_THREAD_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_CHRONO_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_DATE_TIME_LIBRARY)
    copy_dependency_with_config(${GAPCA_MEX_Project} Boost_PROGRAM_OPTIONS_LIBRARY)
  endif()


endif()




##############################################################################################
# Matlab unit tests

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

  matlab_add_unit_test(
    NAME ${PROJECT_NAME}_matlabtest-1
    TIMEOUT 180
    UNITTEST_FILE ${CMAKE_SOURCE_DIR}/test/grassmannpca_matlab_unit_tests.m
    ADDITIONAL_PATH $<TARGET_FILE_DIR:${GAPCA_MEX_Project}>
    )
  
  matlab_add_unit_test(
    NAME ${PROJECT_NAME}_matlabtest-2
    TIMEOUT 600
    UNITTEST_FILE ${CMAKE_SOURCE_DIR}/test/grassmannpca_matlab_performance_unit_tests.m
    ADDITIONAL_PATH $<TARGET_FILE_DIR:${GAPCA_MEX_Project}>
    TEST_ARGS CONFIGURATIONS Release
    )
  
  
endif()




##############################################################################################
# Installation rules

# CPack rules
set(CPACK_PACKAGE_VENDOR "Max Planck Institute for Intelligent Systems")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Grassmann Averages for computing a scalable and robust PCA, Matlab extensions")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_STRIP_FILES TRUE)
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/extensions/installer_welcome.txt")

if(WIN32)
  set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
  set(CPACK_NSIS_URL_INFO_ABOUT http://ps.is.tuebingen.mpg.de/project/Robust_PCA)
  set(CPACK_NSIS_MENU_LINKS
       "documentation/index.html" "C++ library documentation"
       "README.txt" "README.txt"
        "http://ps.is.tuebingen.mpg.de/project/Robust_PCA" "Grassmann Averages Web Site")  
else()
  set(CPACK_GENERATOR TBZ2)
endif()

include(CPack)

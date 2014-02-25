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
add_library(robust_pca ${RobustPCA_LIB})




# testing using boost
set(Boost_DEBUG ON)
set(Boost_DETAILED_FAILURE_MSG ON)

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





set(RobustPCA_Matlab
     extension/matlab/robust_pca.mex)

set(MATLAB_VERSIONS_MAPPING
  "R2013b" "8.2"
  "R2013a" "8.1"
  "R2012b" "8.0"
  "R2012a" "7.14"

  "R2011b" "7.13"
  "R2011a" "7.12"
  "R2010b" "7.11"
  
  ${MATLAB_ADDITIONAL_VERSIONS}
  
  CACHE STRING "Supported matlab versions")



# get the version of matlab (17.58) from a release name (R2017k)
macro (get_matlab_version_from_release_name release_name version_name)
  list(FIND MATLAB_VERSIONS_MAPPING ${release_name} index)
  if(${index} EQUAL -1)
    message(WARNING "The release name ${release_name} is not registered")
  endif()
  math(EXPR index "${index}+1")
  list(GET MATLAB_VERSIONS_MAPPING ${index} version)
  set(${version_name} ${version})
endmacro(get_matlab_version_from_release_name)

# extracts all the supported release names (R2017k...) of matlab
macro(extract_matlab_releases list_releases)
  list(LENGTH MATLAB_VERSIONS_MAPPING versions_length)
  math(EXPR versions_length "${versions_length}-1")
  set(${list_releases})
  foreach(matlab_release RANGE 0 ${versions_length} 2)
    list(GET MATLAB_VERSIONS_MAPPING ${matlab_release} current)
    list(APPEND ${list_releases} ${current})
  endforeach(matlab_release)
endmacro(extract_matlab_releases)

# extracts all the supported versions of matlab
macro(extract_matlab_versions list_versions)
  list(LENGTH MATLAB_VERSIONS_MAPPING versions_length)
  set(${list_versions})
  foreach(matlab_version RANGE 1 ${versions_length} 2)
    list(GET MATLAB_VERSIONS_MAPPING ${matlab_version} current)
    list(APPEND ${list_versions} ${current})
  endforeach(matlab_version)
endmacro(extract_matlab_versions)



# this function parses the registry and founds the matlab versions that are "really" installed
# the other approach is to use a brute force search
# set win64 to TRUE if the 64 bit version of matlab should be looked for
# the returned list contains all versions under HKLM\\SOFTWARE\\Mathworks\\MATLAB or an empty list in case an error occurred (or nothing found)
#
# Only the version is provided: no path, no test for existence
macro(extract_matlab_versions_from_registry win64 matlab_versions)
  
  if(NOT CMAKE_HOST_WIN32)
    message(FATAL_ERROR "This macro can only be called by a windows host (call to reg.exe")
  endif()
  
  #message(STATUS "System processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
  
  
  
  # list the keys under HKEY_LOCAL_MACHINE\SOFTWARE\mathworks but the call to reg does not work
  # from cmake, curiously, as is. The command provides the desired result under the command line though.
  # Fix: this is because "/reg:64" should appended to the command, otherwise it gets on the 32 bit software key (curiously again)
  find_program(REG_LOCATION "reg")
  file(TO_NATIVE_PATH ${REG_LOCATION} REG_LOCATION)
  
  # if reg.exe is not found, then it is impossible to use this method.
  if(NOT REG_LOCATION)
    message(STATUS "[MATLAB] reg.exe not found")
    set(${matlab_versions})
    return()
  endif()
  
  
  
  if(${win64} AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "64")
    set(APPEND_REG "/reg:64")
  else()
    set(APPEND_REG "/reg:32")
  endif()
  
  #message(STATUS "APPEND_REG ${APPEND_REG}")
  
  # /reg:64 should be added on 64 bits capable OSs in order to enable the redirection of 64 bits applications
  execute_process(
    COMMAND ${REG_LOCATION} query HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\MATLAB /f * /k ${APPEND_REG}
    RESULT_VARIABLE resultMatlab
    OUTPUT_VARIABLE varMatlab)
  #message("Matlabs = ${varMatlab} | ${resultMatlab}")
  
  
  set(matlabs_from_registry)
  if(${resultMatlab} EQUAL 0)
    # just for tests
    #set(varMatlab "${varMatlab} HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\MATLAB\\7.47 HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\MATLAB\\9")
    string(
      REGEX MATCHALL "MATLAB\\\\([0-9]+(\\.[0-9]+)?)"
      matlab_versions_regex ${varMatlab})
    #message(STATUS "regex = ${matlab_versions_regex}")
    foreach(match IN LISTS matlab_versions_regex)
      string(
        REGEX MATCH "MATLAB\\\\(([0-9]+)(\\.([0-9]+))?)"
        current_match ${match})
      #message(STATUS "current match is ${CMAKE_MATCH_0} ${CMAKE_MATCH_1} ${CMAKE_MATCH_2} ${CMAKE_MATCH_3} ${CMAKE_MATCH_4}")
      
      set(current_matlab_version ${CMAKE_MATCH_1})
      set(current_matlab_version_major ${CMAKE_MATCH_2})
      set(current_matlab_version_minor ${CMAKE_MATCH_4})
      if(NOT current_matlab_version_minor)
        set(current_matlab_version_minor "0")
      endif()
      
      #message(STATUS "Matlab in registry ${current_matlab_version_major}.${current_matlab_version_minor}")

      list(APPEND matlabs_from_registry ${current_matlab_version})
    endforeach(match)
    
  endif()
  
  set(${matlab_versions} ${matlabs_from_registry})

endmacro(extract_matlab_versions_from_registry)

macro(extract_matlab_versions_from_registry_brute_force matlab_versions)
  # get the supported versions
  set(matlab_supported_versions)
  extract_matlab_versions(matlab_supported_versions)
  
  
  # this is a manual population of the versions we want to look for
  # this can be done as is, but preferably with the call to 
  # extract_matlab_versions and variable 
  
  # populating the versions we want to look for
  # set(matlab_supported_versions)
  
  # # matlab 7
  # set(matlab_major 7)
  # foreach(current_matlab_minor RANGE 4 20)
    # list(APPEND matlab_supported_versions "${matlab_major}.${current_matlab_minor}")
  # endforeach(current_matlab_minor)

  # # matlab 8
  # set(matlab_major 8)
  # foreach(current_matlab_minor RANGE 0 5)
    # list(APPEND matlab_supported_versions "${matlab_major}.${current_matlab_minor}")
  # endforeach(current_matlab_minor)
  
  # # taking into account the possible additional versions provided by the user
  # if(DEFINED MATLAB_ADDITIONAL_VERSIONS)
    # list(APPEND matlab_supported_versions MATLAB_ADDITIONAL_VERSIONS)
  # endif()
  
  
  # we order from more recent to older
  list(REMOVE_DUPLICATES matlab_supported_versions)
  list(SORT matlab_supported_versions)
  list(REVERSE matlab_supported_versions)
  
  
  set(${matlab_versions} ${matlab_supported_versions})
  

endmacro(extract_matlab_versions_from_registry_brute_force)


# populates the matlab root with valid versions of matlab. The matlab_versions comes either from
# extract_matlab_versions_from_registry_brute_force or extract_matlab_versions_from_registry.
# The returned matlab_roots is organized in pairs version_number,matlab_root_path.  
macro(get_all_valid_matlab_roots matlab_versions matlab_roots)
  
  set(matlab_roots_list )
  foreach(current_matlab_version ${matlab_versions})
    get_filename_component(
      current_MATLAB_ROOT
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${current_matlab_version};MATLABROOT]"
      ABSOLUTE)
      
    if(EXISTS ${current_MATLAB_ROOT})
      list(APPEND matlab_roots_list ${current_matlab_version} ${current_MATLAB_ROOT})
    else()
      #message(WARNING "[MATLAB] The Matlab version found in the registry (${current_matlab_version}) does not exist on the drives. Skipping")
    endif()
  
  endforeach(current_matlab_version)
  set(matlab_roots ${matlab_roots_list})
endmacro(get_all_valid_matlab_roots)




set(toto)
get_matlab_version_from_release_name("R2012a" toto)
message("version is ${toto}")

extract_matlab_releases(toto)
message("ReleasesS is ${toto}")

extract_matlab_versions(toto)
message("versionSS is ${toto}")

# listing the matlab versions installed on the WIN machine if MATLAB_ROOT is not set
# does not work
if(WIN32)
  
  # for windows
  
  if(NOT DEFINED MATLAB_ROOT OR NOT ${MATLAB_ROOT})
    # testing if we are able to extract the needed information from the registry
    set(matlab_versions_from_registry)
    extract_matlab_versions_from_registry(TRUE matlab_versions_from_registry)
    
    
    if(NOT matlab_versions_from_registry)
      # the returned list is empty, doing the search on all known versions
      message(STATUS "[MATLAB] Smart search of Matlab not found in the registry, testing brute force")
      extract_matlab_versions_from_registry_brute_force(matlab_versions_from_registry)
    endif()
    
    #message(FATAL_ERROR "matlab_versions_from_registry ${matlab_versions_from_registry}")
    
    # filtering the results
    set(matlab_roots)
    get_all_valid_matlab_roots(${matlab_versions_from_registry} matlab_roots)
    message(STATUS "Matlab ROOTS ${matlab_roots}")
    
    
  elseif(NOT EXISTS ${MATLAB_ROOT})
    message(FATAL_ERROR "[MATLAB] the specified path for MATLAB_ROOT does not exist (${MATLAB_ROOT})")
  endif()
 
  
else()
  # for linux/osx

  if(NOT DEFINED MATLAB_ROOT OR NOT ${MATLAB_ROOT})
  
    set(matlab_roots)

    # maybe using CMAKE_PROGRAM_PATH to add some more hints
    find_program(
      MATLAB_PROGRAM
      "matlab")
    if(${MATLAB_PROGRAM})
      message(STATUS "[MATLAB] found from the command line at ${MATLAB_PROGRAM}")
      
      # resolve symlinks
      get_filename_component(current_matlab_location ${MATLAB_PROGRAM} REALPATH)
      
      # get the directory
      get_filename_component(current_matlab_location ${current_matlab_location} DIRECTORY)
      get_filename_component(current_matlab_location ${current_matlab_location} DIRECTORY) # matlab should be in bin
      list(APPEND matlab_roots ${current_matlab_location})
    else()
      #message(STATUS "[MATLAB] not found from the command line")
    endif()
    
    # on mac, we look for the /Application paths
    # this corresponds to the behaviour on Windows
    if((NOT matlab_roots) AND APPLE)
      set(matlab_versions)
      extract_matlab_releases(matlab_versions)
      message(STATUS "[MATLAB] Matlab supported versions ${matlab_versions}. If more version should be supported "
                     "the variable MATLAB_ADDITIONAL_VERSIONS can be set according to the documentation")
      

      foreach(current_matlab_version IN LISTS matlab_versions)
        set(matlab_full_string "/Applications/MATLAB_${current_matlab_version}.app")
        if(EXISTS ${matlab_full_string})
          list(APPEND matlab_roots ${matlab_full_string})
        endif()
      endforeach(current_matlab_version)
    
    endif()
    
    message(STATUS "[MATLAB] Matlab roots found: ${matlab_roots}")
    

    
  elseif(NOT EXISTS ${MATLAB_ROOT})
    message(FATAL_ERROR "[MATLAB] the specified path for MATLAB_ROOT does not exist (${MATLAB_ROOT})")
  endif()
endif()


# take the first possible matlab root
if(NOT MATLAB_ROOT AND matlab_roots)
  list(GET matlab_roots 0 MATLAB_ROOT)
  list(LENGTH matlab_roots numbers_of_matlab_roots)
  if(numbers_of_matlab_roots GREATER 1)
    message(WARNING "[MATLAB] Found several distributions of Matlab. Setting the current version to ${MATLAB_ROOT}. If this is not"
                    "a desirable behaviour, provide the -DMATLAB_ROOT on the command line")
  endif()
endif()

file(TO_CMAKE_PATH ${MATLAB_ROOT} MATLAB_ROOT)

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(64Build FALSE)
else()
  set(64Build TRUE)
endif()

if(APPLE)
  set(MATLAB_BIN1 "mac") # i should be for intel
  set(MATLAB_suffix32 "i")
  set(MATLAB_suffix64 "i64")
elseif(UNIX)
  set(MATLAB_BIN1 "gln")
  set(MATLAB_suffix32 "x86")
  set(MATLAB_suffix64 "xa64")
else()
  set(MATLAB_BIN1 "win")
  set(MATLAB_suffix32 "32")
  set(MATLAB_suffix64 "64")
endif()

set(MATLAB_BIN_DIR ${MATLAB_ROOT}/bin)
set(MATLAB_INCLUDE_DIR ${MATLAB_ROOT}/extern/include)
if(64Build)
  set(MATLAB_BIN_DIR_ARCH ${MATLAB_ROOT}/bin/${MATLAB_BIN1}${MATLAB_suffix64}  CACHE PATH "Matlab directory for architecture specific binaries")
  set(MATLAB_EXTERN_LIB_DIR ${MATLAB_ROOT}/extern/lib/${MATLAB_BIN1}${MATLAB_suffix64} CACHE PATH "Matlab directory for link")
else()
  set(MATLAB_BIN_DIR_ARCH ${MATLAB_ROOT}/bin/${MATLAB_BIN1}${MATLAB_suffix32} CACHE PATH "Matlab directory for architecture specific binaries" )
  set(MATLAB_EXTERN_LIB_DIR ${MATLAB_ROOT}/extern/lib/${MATLAB_BIN1}${MATLAB_suffix64} CACHE PATH "Matlab directory for link")
endif()


# todo setup the extension properly
set(mexext_suffix "")
if(WIN32)
  list(APPEND mexext_suffix ".bat")
endif()

foreach(current_mexext_suffix IN LISTS mexext_suffix)
  if(NOT DEFINED MATLAB_MEXEXTENSIONS_PROG OR NOT MATLAB_MEXEXTENSIONS_PROG)
    find_program(
      MATLAB_MEXEXTENSIONS_PROG
      "mexext${current_mexext_suffix}"
      PATHS ${MATLAB_BIN_DIR}
      DOC "Matlab MEX extension provider"
      NO_DEFAULT_PATH
    )
  endif()
endforeach(current_mexext_suffix)

# Raffi/ I think this one may be done automatically
set(MEX_EXTENSION "" CACHE STRING "Extensions for the mex targets")

if(NOT MEX_EXTENSION)
  execute_process(
    COMMAND ${MATLAB_MEXEXTENSIONS_PROG} 
    OUTPUT_VARIABLE MEX_EXTENSION)
  string(STRIP ${MEX_EXTENSION} MEX_EXTENSION)
  set(MEX_EXTENSION ${MEX_EXTENSION} CACHE STRING "Extensions for the mex targets (automatically given by Matlab)")
  message(STATUS "MEX_EXTENSION ${MEX_EXTENSION}")
endif()


find_program(
  MATLAB_MEX_PROG
  "mex"
  PATHS ${MATLAB_BIN_DIR_ARCH}
  DOC "Matlab MEX compiler"
  NO_DEFAULT_PATH
)




message("MATLAB_MEXEXTENSIONS_PROG = ${MATLAB_MEXEXTENSIONS_PROG}")


#find_package(Matlab)
find_library(MATLAB_MEX_LIBRARY
  mex
  PATHS ${MATLAB_EXTERN_LIB_DIR}
  PATH_SUFFIXES "microsoft" "."
  NO_DEFAULT_PATH
)
find_library(MATLAB_MX_LIBRARY
  mx
  PATHS ${MATLAB_EXTERN_LIB_DIR}
  PATH_SUFFIXES "microsoft" "."
  NO_DEFAULT_PATH
)
find_library(MATLAB_ENG_LIBRARY
  eng
  PATHS ${MATLAB_EXTERN_LIB_DIR}
  PATH_SUFFIXES "microsoft" "."
  NO_DEFAULT_PATH
)

find_path(MATLAB_INCLUDE_DIR2
  "mex.h"
  PATHS ${MATLAB_INCLUDE_DIR}
  NO_DEFAULT_PATH
  )

message("MATLAB_INCLUDE_DIR = ${MATLAB_INCLUDE_DIR} | ${MATLAB_INCLUDE_DIR2}")
if(NOT MATLAB_INCLUDE_DIR)
  message(FATAL_ERROR "Unable to find Matlab") # todo make it optional
endif()

include_directories(${MATLAB_INCLUDE_DIR})
set(RobustPCAMatlab_LIB
  extensions/matlab.cpp)
add_library(robust_pca_m SHARED ${RobustPCAMatlab_LIB})
target_link_libraries(robust_pca_m ${MATLAB_MEX_LIBRARY})


set_target_properties(robust_pca_m PROPERTIES SUFFIX ${MEX_EXTENSION})

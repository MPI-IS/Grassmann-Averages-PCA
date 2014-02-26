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


# returns the extension of the mex (the suffixes) 
macro(get_matlab_suffix matlab_root mex_suffix)

  # todo setup the extension properly
  set(mexext_suffix "")
  if(WIN32)
    list(APPEND mexext_suffix ".bat")
  endif()

  # we first try without suffix, since cmake does not understand a list with one empty string element
  find_program(
    MATLAB_MEXEXTENSIONS_PROG
    "mexext"
    PATHS ${matlab_root}/bin
    DOC "Matlab MEX extension provider"
    NO_DEFAULT_PATH
  )

  foreach(current_mexext_suffix IN LISTS mexext_suffix)
    if(NOT DEFINED MATLAB_MEXEXTENSIONS_PROG OR NOT MATLAB_MEXEXTENSIONS_PROG)
      # this call should populate the cache automatically
      find_program(
        MATLAB_MEXEXTENSIONS_PROG
        "mexext${current_mexext_suffix}"
        PATHS ${matlab_root}/bin
        DOC "Matlab MEX extension provider"
        NO_DEFAULT_PATH
      )
    endif()
  endforeach(current_mexext_suffix)
  
  if((NOT MATLAB_MEXEXTENSIONS_PROG) OR (NOT EXISTS ${MATLAB_MEXEXTENSIONS_PROG}))
    message(FATAL_ERROR "[MATLAB] Cannot found mexext program. Matlab root is ${matlab_root}")
  endif()

  set(MEX_EXTENSION "" CACHE STRING "Extensions for the mex targets")

  if(NOT MEX_EXTENSION)
    execute_process(
      COMMAND ${MATLAB_MEXEXTENSIONS_PROG} 
      OUTPUT_VARIABLE MEX_EXTENSION)
    string(STRIP ${MEX_EXTENSION} MEX_EXTENSION)
    set(MEX_EXTENSION ${MEX_EXTENSION} CACHE STRING "Extensions for the mex targets (automatically given by Matlab)")
    #message(STATUS "MEX_EXTENSION ${MEX_EXTENSION}")
  endif()  
  
  set(${mex_suffix} ${MEX_EXTENSION})
endmacro(get_matlab_suffix)



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
      extract_matlab_releases(matlab_releases)
      message(STATUS "[MATLAB] Matlab supported versions ${matlab_releases}. If more version should be supported "
                     "the variable MATLAB_ADDITIONAL_VERSIONS can be set according to the documentation")
      

      foreach(current_matlab_release IN LISTS matlab_releases)
        set(matlab_full_string "/Applications/MATLAB_${current_matlab_release}.app")
        if(EXISTS ${matlab_full_string})
          set(current_matlab_version)
          get_matlab_version_from_release_name(${current_matlab_release} current_matlab_version)
          list(APPEND matlab_roots ${current_matlab_version})
          list(APPEND matlab_roots ${matlab_full_string})
        endif()
      endforeach(current_matlab_release)
    
    endif()
    
    message(STATUS "[MATLAB] Matlab roots found: ${matlab_roots}")
    

    
  elseif(NOT EXISTS ${MATLAB_ROOT})
    message(FATAL_ERROR "[MATLAB] the specified path for MATLAB_ROOT does not exist (${MATLAB_ROOT})")
  endif()
endif()


# take the first possible matlab root
if(NOT MATLAB_ROOT AND matlab_roots)
  list(GET matlab_roots 0 MATLAB_VERSION)
  list(GET matlab_roots 1 MATLAB_ROOT)
  list(LENGTH matlab_roots numbers_of_matlab_roots)
  
  # adding a warning in case of ambiguity
  if(numbers_of_matlab_roots GREATER 2)
    message(WARNING "[MATLAB] Found several distributions of Matlab. Setting the current version to ${MATLAB_VERSION} (located ${MATLAB_ROOT})."
                    " If this is not the desired behaviour, provide the -DMATLAB_ROOT on the command line")
  endif()
endif()

# todo: if the version is empty, we should determine it from the matlab executable
if(NOT MATLAB_VERSION)
  if((NOT DEFINED MATLAB_PROGRAM) OR (NOT ${MATLAB_PROGRAM}) OR (NOT EXISTS ${MATLAB_PROGRAM}))
    message(STATUS "[MATLAB] Looking for Matlab under ${MATLAB_ROOT}")
    find_program(
      MATLAB_PROGRAM
      matlab
      PATHS ${MATLAB_ROOT} ${MATLAB_ROOT}/bin
      DOC "Matlab main program"
      NO_DEFAULT_PATH
    )
    
    
    # todo put this in an appropriate macro, since it is (theoretically) platform independent.
    if(MATLAB_PROGRAM)
      # timeout set to 30 seconds, in case it does not start
      message(STATUS "[MATLAB] Determining the version of Matlab")
      
      #  
      # the log file is needed since on windows the command executes in a new window and it is not possible 
      # to get back the answer of Matlab
      # the -wait command is needed on windows, otherwise the call returns immediately after the program launches itself.
      if(WIN32)
        set(matlab_additional_commands "-wait")
      endif()
      
      # note as said before OUTPUT_VARIABLE cannot be used in a platform independent manner
      # however, not setting it would flush the output of matlab in the current console (unix variant)
      execute_process(
        COMMAND ${MATLAB_PROGRAM} -nosplash -nojvm ${matlab_additional_commands} -logfile ${CMAKE_BINARY_DIR}/matlabVersionLog.cmaketmp -nodesktop -nodisplay -r "version, exit" 
        OUTPUT_VARIABLE MATLAB_VERSION_FROM_CMD_dummy
        RESULT_VARIABLE resultMatlabVersionCall
        TIMEOUT 30
        )
      # read back the log
      file(READ ${CMAKE_BINARY_DIR}/matlabVersionLog.cmaketmp MATLAB_VERSION_FROM_CMD)
      
      if(${resultMatlabVersionCall})
        message(WARNING "[MATLAB] Unable to determine the version of Matlab. Matlab output is ${MATLAB_VERSION_FROM_CMD}. Returned with error ${resultMatlabVersionCall}.")
      else()
        set(index -1)
        string(FIND ${MATLAB_VERSION_FROM_CMD} "ans" index)
        if(index EQUAL -1)
          message(WARNING "[MATLAB] Cannot find the version of Matlab answered from Matlab")
        else()
          set(list_of_all_versions)
          string(SUBSTRING ${MATLAB_VERSION_FROM_CMD} ${index} -1 substring_ans)
          string(
            REGEX MATCHALL "ans[\r\n\t ]*=[\r\n\t ]*([0-9]+(\\.[0-9]+)?)"
            matlab_versions_regex 
            ${substring_ans})
          foreach(match IN LISTS matlab_versions_regex)
            string(
              REGEX MATCH "ans[\r\n\t ]*=[\r\n\t ]*(([0-9]+)(\\.([0-9]+))?)"
              current_match ${match})
            #message(STATUS "current match is ${CMAKE_MATCH_0} ${CMAKE_MATCH_1} ${CMAKE_MATCH_2} ${CMAKE_MATCH_3} ${CMAKE_MATCH_4}")          
            list(APPEND list_of_all_versions ${CMAKE_MATCH_1})
          endforeach()
          list(REMOVE_DUPLICATES list_of_all_versions)
          list(GET list_of_all_versions 0 MATLAB_VERSION_tmp)
          
          # set the version into the cache
          set(MATLAB_VERSION ${MATLAB_VERSION_tmp} CACHE STRING "Matlab version (automatically determined)")
          list(LENGTH list_of_all_versions list_of_all_versions_length)
          if(${list_of_all_versions_length} GREATER 1)
            message(WARNING "[MATLAB] Found several versions, taking the first one (versions found ${list_of_all_versions})")
          endif()
        endif()
        
        #message(STATUS "[MATLAB] Matlab found in ${MATLAB_PROGRAM} ${substring_ans} / result = ${resultMatlabVersionCall} / regex = ${matlab_versions_regex}")
      endif()
    endif()
  endif()
endif()

message(STATUS "[MATLAB] Current version is ${MATLAB_VERSION} located ${MATLAB_ROOT}")
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
set(MATLAB_INCLUDE_DIR_TO_LOOK ${MATLAB_ROOT}/extern/include)
if(64Build)
  set(MATLAB_BIN_DIR_ARCH ${MATLAB_ROOT}/bin/${MATLAB_BIN1}${MATLAB_suffix64}  CACHE PATH "Matlab directory for architecture specific binaries")
  set(MATLAB_EXTERN_LIB_DIR ${MATLAB_ROOT}/extern/lib/${MATLAB_BIN1}${MATLAB_suffix64} CACHE PATH "Matlab directory for link")
else()
  set(MATLAB_BIN_DIR_ARCH ${MATLAB_ROOT}/bin/${MATLAB_BIN1}${MATLAB_suffix32} CACHE PATH "Matlab directory for architecture specific binaries" )
  set(MATLAB_EXTERN_LIB_DIR ${MATLAB_ROOT}/extern/lib/${MATLAB_BIN1}${MATLAB_suffix64} CACHE PATH "Matlab directory for link")
endif()

if(WIN32)
  set(MATLAB_LIB_DIR_FOR_LOOKUP ${MATLAB_EXTERN_LIB_DIR}/microsoft)
  set(MATLAB_LIB_PREFIX_FOR_LOOKUP "lib")
else()
  set(MATLAB_LIB_DIR_FOR_LOOKUP ${MATLAB_BIN_DIR_ARCH})
  set(MATLAB_LIB_PREFIX_FOR_LOOKUP "lib")
endif()


set(MEX_EXTENSION)
get_matlab_suffix(${MATLAB_ROOT} current_mex_extension)



find_program(
  MATLAB_MEX_PROG
  "mex"
  PATHS ${MATLAB_BIN_DIR_ARCH}
  DOC "Matlab MEX compiler"
  NO_DEFAULT_PATH
)


set(toto)
get_matlab_version_from_release_name("R2012a" toto)
message("version is ${toto}")

extract_matlab_releases(toto)
message("ReleasesS is ${toto}")

extract_matlab_versions(toto)
message("versionSS is ${toto}")



message("MATLAB_LIB_PREFIX_FOR_LOOKUP ${MATLAB_LIB_PREFIX_FOR_LOOKUP} | MATLAB_LIB_DIR_FOR_LOOKUP = ${MATLAB_LIB_DIR_FOR_LOOKUP}")

set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} ${MATLAB_LIB_PREFIX_FOR_LOOKUP})

#find_package(Matlab)
find_library(
  MATLAB_MEX_LIBRARY
  mex
  PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
  NO_DEFAULT_PATH
)
find_library(
  MATLAB_MX_LIBRARY
  mx
  PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
  NO_DEFAULT_PATH
)
find_library(
  MATLAB_ENG_LIBRARY
  eng
  PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
  NO_DEFAULT_PATH
)

find_path(
  MATLAB_INCLUDE_DIR
  mex.h
  PATHS ${MATLAB_INCLUDE_DIR_TO_LOOK}
  NO_DEFAULT_PATH
  )



set(MATLAB_LIBRARIES ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY})

message(STATUS "[MATLAB] libraries found ${MATLAB_LIBRARIES}")

if(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)
  mark_as_advanced(
    MATLAB_LIBRARIES
    MATLAB_MEX_LIBRARY
    MATLAB_MX_LIBRARY
    MATLAB_ENG_LIBRARY
    MATLAB_INCLUDE_DIR
    MATLAB_FOUND
    MATLAB_ROOT
    MATLAB_VERSION
  )

  set(MATLAB_FOUND TRUE)
endif()




  
message("MATLAB_INCLUDE_DIR = ${MATLAB_INCLUDE_DIR}")
if(NOT MATLAB_INCLUDE_DIR)
  message(FATAL_ERROR "Unable to find Matlab include directory")
endif()

include_directories(${MATLAB_INCLUDE_DIR})
set(RobustPCAMatlab_LIB
  extensions/matlab.cpp)
add_library(robust_pca_m SHARED ${RobustPCAMatlab_LIB})
target_link_libraries(robust_pca_m ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY})


set_target_properties(robust_pca_m PROPERTIES SUFFIX ".${MEX_EXTENSION}")
set_target_properties(robust_pca_m PROPERTIES COMPILE_DEFINITIONS "DLL_EXPORT_SYM=__declspec(dllexport)")
#.rst:
# FindMatlab
# ----------------
#
# Finds Matlab installations and provides Matlab tools and libraries to cmake.
# 
# This package first intention is finding the libraries associated with Matlab in order
# to be able to compile Matlab extensions (mex files). It can also be used to run unit test on these mex extensions,
# and run Matlab.
#
# The variable ``MATLAB_ROOT`` may be specified in order to give the path of the desired Matlab version. Also, 
# additionnal information is provided when ``MATLAB_FIND_DEBUG`` is set.
# Otherwise, the behaviour is platform dependant:
#
# - on Windows, the installed versions of Matlab are retrieved from the Windows registry
# - on Mac, the installed versions of Matlab are given by the MATLAB pathes in /Application
# - on Unix, the desired Matlab should be accessible from the PATH.
#
# When a Matlab binary is found and the ``MATLAB_VERSION`` is not given, then the
# version is queried from Matlab directly. On Windows, it can make a window running Matlab appear.
#
# The mapping of the release names and the version of Matlab is performed by defining pairs (name, version). 
# The variable ``MATLAB_ADDITIONAL_VERSIONS`` may be provided in order to handle additional versions, in the following form:
#
# ::
#
#    set(MATLAB_ADDITIONAL_VERSIONS 
#       "release_name1" "corresponding_version1"
#       "release_name2" "corresponding_version2"
#       ...
#       )
# 
# such as 
# ::
# 
#    set(MATLAB_ADDITIONAL_VERSIONS
#        "R2013b" "8.2"
#        "R2013a" "8.1"
#        "R2012b" "8.0")
#
#
# Defined variables
# -----------------
# * ``MATLAB_FOUND`` true if the Matlab installation is found.
# * ``MATLAB_ROOT`` the root of the Matlab installation
# * ``MATLAB_VERSION`` the version of the Matlab installation
# * ``MATLAB_PROGRAM`` the Matlab binary program. Available only if the component ``MAIN_PROGRAM`` is asked
# * ``MATLAB_INCLUDE_DIR`` the path of the Matlab libraries headers
# * ``MATLAB_MEX_LIBRARY`` library for mex
# * ``MATLAB_MX_LIBRARY`` mx library of Matlab (arrays). Available only if the component ``MX_LIBRARY`` is asked
# * ``MATLAB_ENG_LIBRARY`` Matlab engine library. Available only if the component ``ENG_LIBRARY`` is asked
# * ``MATLAB_LIBRARIES`` the whole set of libraries of Matlab
# * ``MATLAB_MEX_COMPILER`` the mex compiler of Matlab. Currently not used internally. Available only if the component ``MEX_COMPILER`` is asked
#
#
# Defined macros
# --------------
# * ``matlab_get_version_from_release_name`` returns the version from the release name
# * ``matlab_extract_all_installed_versions_from_registry`` parses the registry for all Matlab versions. Available on Windows only. 
#   The part of the registry parsed is dependent on the host processor 
# * ``matlab_get_all_valid_matlab_roots_from_registry`` returns all the possible Matlab paths, according to a previously given list. Only the
#   existing/accessible paths are kept. This is mainly useful for the "brute force" search of Matlab installation.
# * ``matlab_get_mex_suffix`` returns the suffix to be used for the mex files (platform/architecture dependant)
# * ``matlab_get_version_from_matlab_run`` returns the version of Matlab, given the full directory of the Matlab program.
#
# 
# Future work
# -----------
# - win32:an additional variable telling that the registry is x86 or x64, maybe depending on the target build.
# - add a unit test method for Matlab versions >= 8.1 (R2013+)

#=============================================================================
# Copyright 2009-2014 Raffi Enficiaud
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

include(FindPackageHandleStandardArgs)


# The currently supported versions. Other version can be added by the user by providing MATLAB_ADDITIONAL_VERSIONS
if(NOT MATLAB_ADDITIONAL_VERSIONS)
  set(MATLAB_ADDITIONAL_VERSIONS)
endif()

set(MATLAB_VERSIONS_MAPPING
  "R2013b" "8.2"
  "R2013a" "8.1"
  "R2012b" "8.0"
  "R2012a" "7.14"

  "R2011b" "7.13"
  "R2011a" "7.12"
  "R2010b" "7.11"
  
  ${MATLAB_ADDITIONAL_VERSIONS}
  )


# temporary folder for all Matlab runs
set(_matlab_temporary_folder ${CMAKE_BINARY_DIR}/Matlab)
  
if(NOT EXISTS ${_matlab_temporary_folder})
  file(MAKE_DIRECTORY ${_matlab_temporary_folder})
endif()

# get the version of Matlab (17.58) from a release name (R2017k)
macro (matlab_get_version_from_release_name release_name version_name)
  list(FIND MATLAB_VERSIONS_MAPPING ${release_name} index)
  if(${index} EQUAL -1)
    message(WARNING "The release name ${release_name} is not registered")
  endif()
  math(EXPR index "${index}+1")
  list(GET MATLAB_VERSIONS_MAPPING ${index} version)
  set(${version_name} ${version})
endmacro(matlab_get_version_from_release_name)

# extracts all the supported release names (R2017k...) of Matlab
macro(matlab_get_supported_releases list_releases)
  list(LENGTH MATLAB_VERSIONS_MAPPING versions_length)
  math(EXPR versions_length "${versions_length}-1")
  set(${list_releases})
  foreach(matlab_release RANGE 0 ${versions_length} 2)
    list(GET MATLAB_VERSIONS_MAPPING ${matlab_release} current)
    list(APPEND ${list_releases} ${current})
  endforeach(matlab_release)
endmacro(matlab_get_supported_releases)

# extracts all the supported versions of Matlab
macro(matlab_get_supported_versions list_versions)
  list(LENGTH MATLAB_VERSIONS_MAPPING versions_length)
  set(${list_versions})
  foreach(matlab_version RANGE 1 ${versions_length} 2)
    list(GET MATLAB_VERSIONS_MAPPING ${matlab_version} current)
    list(APPEND ${list_versions} ${current})
  endforeach(matlab_version)
endmacro(matlab_get_supported_versions)



# This function parses the registry and founds the Matlab versions that are "really" installed
# the other approach is to use a brute force search
# set win64 to TRUE if the 64 bit version of Matlab should be looked for
# the returned list contains all versions under HKLM\\SOFTWARE\\Mathworks\\MATLAB or an empty list in case an error occurred (or nothing found)
#
# Only the version is provided: no path, no test for existence
function(matlab_extract_all_installed_versions_from_registry win64 matlab_versions)
  
  if(NOT CMAKE_HOST_WIN32)
    message(FATAL_ERROR "This macro can only be called by a windows host (call to reg.exe")
  endif()
  
  #message(STATUS "System processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
  
  
  
  # list the keys under HKEY_LOCAL_MACHINE\SOFTWARE\mathworks but the call to reg does not work
  # from cmake, curiously, as is. The command provides the desired result under the command line though.
  # Fix: this is because "/reg:64" should appended to the command, otherwise it gets on the 32 bits software key (curiously again)
  find_program(MATLAB_REG_EXE_LOCATION "reg")
  file(TO_NATIVE_PATH ${MATLAB_REG_EXE_LOCATION} MATLAB_REG_EXE_LOCATION)
  
  # if reg.exe is not found, then it is impossible to use this method.
  if(NOT MATLAB_REG_EXE_LOCATION)
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] reg.exe not found")
    endif()
    set(${matlab_versions} "" PARENT_SCOPE)
    return()
  endif()
  
  
  
  if(${win64} AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "64")
    set(APPEND_REG "/reg:64")
  else()
    set(APPEND_REG "/reg:32")
  endif()
  
  # /reg:64 should be added on 64 bits capable OSs in order to enable the redirection of 64 bits applications
  execute_process(
    COMMAND ${MATLAB_REG_EXE_LOCATION} query HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\MATLAB /f * /k ${APPEND_REG}
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
  
  set(${matlab_versions} ${matlabs_from_registry} PARENT_SCOPE)

endfunction(matlab_extract_all_installed_versions_from_registry)

macro(extract_matlab_versions_from_registry_brute_force matlab_versions)
  # get the supported versions
  set(matlab_supported_versions)
  matlab_get_supported_versions(matlab_supported_versions)
  
  
  # this is a manual population of the versions we want to look for
  # this can be done as is, but preferably with the call to 
  # matlab_get_supported_versions and variable 
  
  # populating the versions we want to look for
  # set(matlab_supported_versions)
  
  # # Matlab 7
  # set(matlab_major 7)
  # foreach(current_matlab_minor RANGE 4 20)
    # list(APPEND matlab_supported_versions "${matlab_major}.${current_matlab_minor}")
  # endforeach(current_matlab_minor)

  # # Matlab 8
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


# populates the Matlab root with valid versions of Matlab. The matlab_versions comes either from
# extract_matlab_versions_from_registry_brute_force or matlab_extract_all_installed_versions_from_registry.
# The returned matlab_roots is organized in pairs version_number,matlab_root_path.  
function(matlab_get_all_valid_matlab_roots_from_registry matlab_versions matlab_roots)
  
  set(_matlab_roots_list )
  foreach(current_matlab_version ${matlab_versions})
    get_filename_component(
      current_MATLAB_ROOT
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${current_matlab_version};MATLABROOT]"
      ABSOLUTE)
      
    if(EXISTS ${current_MATLAB_ROOT})
      list(APPEND _matlab_roots_list ${current_matlab_version} ${current_MATLAB_ROOT})
    endif()
  
  endforeach(current_matlab_version)
  set(${matlab_roots} ${_matlab_roots_list} PARENT_SCOPE)
endfunction(matlab_get_all_valid_matlab_roots_from_registry)


# returns the extension of the mex files (the suffixes).
# This function should not be called before the appropriate Matlab root has been found.
function(matlab_get_mex_suffix matlab_root mex_suffix)

  # todo setup the extension properly. Currently I do not know if this is sufficient for all win32 distributions.
  # there is also CMAKE_EXECUTABLE_SUFFIX that could be tweaked
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
  
  
  # the program has been found?
  if((NOT MATLAB_MEXEXTENSIONS_PROG) OR (NOT EXISTS ${MATLAB_MEXEXTENSIONS_PROG}))
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] Cannot found mexext program. Matlab root is ${matlab_root}")
    endif()
    return()
  endif()

  set(_matlab_mex_extension)

  
  execute_process(
    COMMAND ${MATLAB_MEXEXTENSIONS_PROG} 
    OUTPUT_VARIABLE _matlab_mex_extension)
  string(STRIP ${_matlab_mex_extension} _matlab_mex_extension)

  set(${mex_suffix} ${_matlab_mex_extension} PARENT_SCOPE)
endfunction(matlab_get_mex_suffix)



# .. command:: matlab_get_version_from_matlab_run
#
#   This function runs Matlab specified on arguments and extracts its version.
#   
#   matlab_get_version_from_matlab_run(matlab_binary_path matlab_version)
function(matlab_get_version_from_matlab_run matlab_binary_program matlab_list_versions)

  set(${matlab_list_versions} "" PARENT_SCOPE)

  
  if(MATLAB_FIND_DEBUG)
    message(STATUS "[MATLAB] Determining the version of Matlab from ${matlab_binary_program}")
  endif()

  if(EXISTS ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
    if(MATLAB_FIND_DEBUG)
      message(STATUS "[MATLAB] Removing previous ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp file")
    endif()
    file(REMOVE ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
  endif()

  
  # the log file is needed since on windows the command executes in a new window and it is not possible 
  # to get back the answer of Matlab
  # the -wait command is needed on windows, otherwise the call returns immediately after the program launches itself.
  if(WIN32)
    set(_matlab_additional_commands "-wait")
  endif()
  
  # timeout set to 30 seconds, in case it does not start
  # note as said before OUTPUT_VARIABLE cannot be used in a platform independent manner
  # however, not setting it would flush the output of Matlab in the current console (unix variant)
  execute_process(
    COMMAND ${matlab_binary_program} -nosplash -nojvm ${_matlab_additional_commands} -logfile ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp -nodesktop -nodisplay -r "version, exit" 
    OUTPUT_VARIABLE _matlab_version_from_cmd_dummy
    RESULT_VARIABLE _matlab_result_version_call
    TIMEOUT 30
    )
  
  
  if(${_matlab_result_version_call})
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] Unable to determine the version of Matlab. Matlab call returned with error ${_matlab_result_version_call}.")
    endif()
    return()
  elseif(NOT EXISTS ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] Unable to determine the version of Matlab. The log file does not exist.")
    endif()
    return()
  endif()

  # if successful, read back the log
  file(READ ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp _matlab_version_from_cmd)
  file(REMOVE ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)

  set(index -1)
  string(FIND ${_matlab_version_from_cmd} "ans" index)
  if(index EQUAL -1)
    
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] Cannot find the version of Matlab returned by the run.")
    endif()
    
  else()
    set(matlab_list_of_all_versions_tmp)
    
    string(SUBSTRING ${_matlab_version_from_cmd} ${index} -1 substring_ans)
    string(
      REGEX MATCHALL "ans[\r\n\t ]*=[\r\n\t ]*([0-9]+(\\.[0-9]+)?)"
      matlab_versions_regex 
      ${substring_ans})
    foreach(match IN LISTS matlab_versions_regex)
      string(
        REGEX MATCH "ans[\r\n\t ]*=[\r\n\t ]*(([0-9]+)(\\.([0-9]+))?)"
        current_match ${match})
      
      list(APPEND matlab_list_of_all_versions_tmp ${CMAKE_MATCH_1})
    endforeach()
    list(REMOVE_DUPLICATES matlab_list_of_all_versions_tmp)
    set(${matlab_list_versions} ${matlab_list_of_all_versions_tmp} PARENT_SCOPE)
    
  endif()
    
endfunction(matlab_get_version_from_matlab_run)


# this variable will get all Matlab installations found in the current system.
set(_matlab_possible_roots)


# listing the Matlab versions installed on the WIN machine if MATLAB_ROOT is not set
if(WIN32)
  
  # for windows
  
  
  if(NOT DEFINED MATLAB_ROOT OR NOT ${MATLAB_ROOT})
    # if MATLAB_ROOT not specified, we look for Matlab installation in the registry
    # if unsuccessful, we look for all known revision and filter the existing ones. 
  
    # testing if we are able to extract the needed information from the registry
    set(matlab_versions_from_registry)
    matlab_extract_all_installed_versions_from_registry(TRUE matlab_versions_from_registry)
    
    # the returned list is empty, doing the search on all known versions
    if(NOT matlab_versions_from_registry)
      
      if(MATLAB_FIND_DEBUG)
        message(STATUS "[MATLAB] Search for Matlab from the registry unsuccessful, testing all supported versions")
      endif()
      
      extract_matlab_versions_from_registry_brute_force(matlab_versions_from_registry)
    endif()
        
    # filtering the results with the registry keys
    matlab_get_all_valid_matlab_roots_from_registry(${matlab_versions_from_registry} _matlab_possible_roots)
    
    
    
  elseif(NOT EXISTS ${MATLAB_ROOT})
  
    # if MATLAB_ROOT specified but erroneous
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] the specified path for MATLAB_ROOT does not exist (${MATLAB_ROOT})")
    endif()
  endif()
 
  
else()
  # for linux/osx

  if(NOT DEFINED MATLAB_ROOT OR NOT ${MATLAB_ROOT})
  
    # if MATLAB_ROOT not specified, we look for Matlab from the command line PATH
    # maybe using CMAKE_PROGRAM_PATH to add some more hints
    find_program(
      MATLAB_PROGRAM
      "matlab")
  
    if(NOT ${MATLAB_PROGRAM})
      #execute_process(COMMAND which matlab OUTPUT_VARIABLE _which_matlab RESULT_VARIABLE _which_matlab_result)
      get_filename_component(MATLAB_PROGRAM "matlab" PROGRAM) 
      if(MATLAB_FIND_DEBUG)
        message(STATUS "[MATLAB] matlab program result from the command line ${MATLAB_PROGRAM}")
      endif()

    endif()
  
    if(MATLAB_PROGRAM AND EXISTS ${MATLAB_PROGRAM})
      if(MATLAB_FIND_DEBUG)
        message(STATUS "[MATLAB] found from the command line at ${MATLAB_PROGRAM}")
      endif()

      # resolve symlinks
      get_filename_component(_matlab_current_location ${MATLAB_PROGRAM} REALPATH)
      if(${CMAKE_VERSION} VERSION_LESS "2.8.12")
        set(_directory_alias PATH)
      else()
        set(_directory_alias DIRECTORY)
      endif()
      # get the directory (the command below has to be run twice)
      get_filename_component(_matlab_current_location ${_matlab_current_location} ${_directory_alias})
      get_filename_component(_matlab_current_location ${_matlab_current_location} ${_directory_alias}) # Matlab should be in bin
      list(APPEND _matlab_possible_roots "NOT-FOUND" ${_matlab_current_location}) # empty version
    endif()
    
    # on mac, we look for the /Application paths
    # this corresponds to the behaviour on Windows. On Linux, we do not have any other guess.
    if((NOT _matlab_possible_roots) AND APPLE)
      
      matlab_get_supported_releases(_matlab_releases)
      if(MATLAB_FIND_DEBUG)
        message(STATUS "[MATLAB] Matlab supported versions ${_matlab_releases}. If more version should be supported "
                       "the variable MATLAB_ADDITIONAL_VERSIONS can be set according to the documentation")
      endif()

      foreach(_current_matlab_release IN LISTS _matlab_releases)
        set(_matlab_full_string "/Applications/MATLAB_${_current_matlab_release}.app")
        if(EXISTS ${_matlab_full_string})
          set(current_matlab_version)
          matlab_get_version_from_release_name(${_current_matlab_release} current_matlab_version)
          list(APPEND _matlab_possible_roots ${_current_matlab_version})
          list(APPEND _matlab_possible_roots ${_matlab_full_string})
        endif()
        
        unset(_matlab_full_string)
      endforeach(_current_matlab_release)
      unset(_current_matlab_release)
      unset(_matlab_releases) 
    endif()

    # we need to clear MATLAB_PROGRAM here
    unset(MATLAB_PROGRAM CACHE)
    unset(MATLAB_PROGRAM) 
    
  elseif(NOT EXISTS ${MATLAB_ROOT})
    # if MATLAB_ROOT specified but erroneous
    if(MATLAB_FIND_DEBUG)
      message(WARNING "[MATLAB] the specified path for MATLAB_ROOT does not exist (${MATLAB_ROOT})")
    endif()
  endif()
endif()


if(MATLAB_FIND_DEBUG)
  message(STATUS "[MATLAB] Matlab root folders are ${_matlab_possible_roots}")
endif()




# take the first possible Matlab root
if(NOT MATLAB_ROOT AND _matlab_possible_roots)
  list(GET _matlab_possible_roots 0 MATLAB_VERSION)
  list(GET _matlab_possible_roots 1 MATLAB_ROOT)
  list(LENGTH _matlab_possible_roots numbers_of_matlab_roots)
  
  # adding a warning in case of ambiguity
  if(numbers_of_matlab_roots GREATER 2)
    message(WARNING "[MATLAB] Found several distributions of Matlab. Setting the current version to ${MATLAB_VERSION} (located ${MATLAB_ROOT})."
                    " If this is not the desired behaviour, provide the -DMATLAB_ROOT on the command line")
  endif()
endif()


if(NOT MATLAB_VERSION OR ${MATLAB_VERSION} STREQUAL "NOT-FOUND")
  if((NOT DEFINED MATLAB_PROGRAM) OR (NOT ${MATLAB_PROGRAM}) OR (NOT EXISTS ${MATLAB_PROGRAM}))
    if(MATLAB_FIND_DEBUG)
      message(STATUS "[MATLAB] - Unknown version, looking for Matlab under ${MATLAB_ROOT}")
    endif()
    message("MATLAB_PROGRAM ${MATLAB_PROGRAM}")
    find_program(
      MATLAB_PROGRAM
      matlab
      PATHS ${MATLAB_ROOT} ${MATLAB_ROOT}/bin
      DOC "Matlab main program"
      NO_DEFAULT_PATH
    )
    message("MATLAB_PROGRAM ${MATLAB_PROGRAM}")

     
    
    if(MATLAB_PROGRAM)
      set(matlab_list_of_all_versions)
      matlab_get_version_from_matlab_run(${MATLAB_PROGRAM} matlab_list_of_all_versions)
    
      list(GET matlab_list_of_all_versions 0 MATLAB_VERSION_tmp)
          
      # set the version into the cache
      set(MATLAB_VERSION ${MATLAB_VERSION_tmp})# CACHE STRING "Matlab version (automatically determined)")
      list(LENGTH list_of_all_versions list_of_all_versions_length)
      if(${list_of_all_versions_length} GREATER 1)
        message(WARNING "[MATLAB] Found several versions, taking the first one (versions found ${list_of_all_versions})")
      endif()
    endif()
  endif()
endif()


if(MATLAB_FIND_DEBUG)
  message(STATUS "[MATLAB] Current version is ${MATLAB_VERSION} located ${MATLAB_ROOT}")
endif()




file(TO_CMAKE_PATH ${MATLAB_ROOT} MATLAB_ROOT)

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(_matlab_64Build FALSE)
else()
  set(_matlab_64Build TRUE)
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
if(_matlab_64Build)
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

unset(_matlab_64Build)


if(NOT DEFINED MATLAB_MEX_EXTENSION)
  set(_matlab_mex_extension "")
  matlab_get_mex_suffix(${MATLAB_ROOT} _matlab_mex_extension)

  # This variable goes to the cache.
  set(MATLAB_MEX_EXTENSION ${_matlab_mex_extension} CACHE STRING "Extensions for the mex targets (automatically given by Matlab)")
  unset(_matlab_mex_extension)
endif()


if(MATLAB_FIND_DEBUG)
  message(STATUS "[MATLAB] [DEBUG]MATLAB_LIB_PREFIX_FOR_LOOKUP = ${MATLAB_LIB_PREFIX_FOR_LOOKUP} | MATLAB_LIB_DIR_FOR_LOOKUP = ${MATLAB_LIB_DIR_FOR_LOOKUP}")
endif()

# WARNING: this thing pollutes the CMAKE_FIND_LIBRARY_PREFIXES global variable. 
# Should it be restored afterwards? Is there a more appropriate way to do that?
set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} ${MATLAB_LIB_PREFIX_FOR_LOOKUP})


set(MATLAB_REQUIRED_VARIABLES)


# the MEX library/header are required
find_path(
  MATLAB_INCLUDE_DIR
  mex.h
  PATHS ${MATLAB_INCLUDE_DIR_TO_LOOK}
  NO_DEFAULT_PATH
  )
list(APPEND MATLAB_REQUIRED_VARIABLES MATLAB_INCLUDE_DIR)

find_library(
  MATLAB_MEX_LIBRARY
  mex
  PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
  NO_DEFAULT_PATH
)
list(APPEND MATLAB_REQUIRED_VARIABLES MATLAB_MEX_LIBRARY)

# the MEX extension is required
list(APPEND MATLAB_REQUIRED_VARIABLES MATLAB_MEX_EXTENSION)


# component Mex Compiler
list(FIND MATLAB_FIND_COMPONENTS MEX_COMPILER _matlab_find_mex_compiler)
if(_matlab_find_mex_compiler GREATER -1)
  find_program(
    MATLAB_MEX_COMPILER
    "mex"
    PATHS ${MATLAB_BIN_DIR_ARCH}
    DOC "Matlab MEX compiler"
    NO_DEFAULT_PATH
  )
  
  if(MATLAB_MEX_COMPILER)
    set(MATLAB_MEX_COMPILER_FOUND TRUE)
  endif()
endif()  

# component Matlab program
list(FIND MATLAB_FIND_COMPONENTS MAIN_PROGRAM _matlab_find_matlab_program)
if(_matlab_find_matlab_program GREATER -1)
  # todo cleanup with code above
  if(NOT DEFINED MATLAB_PROGRAM)
    find_program(
      MATLAB_PROGRAM
      matlab
      PATHS ${MATLAB_ROOT} ${MATLAB_ROOT}/bin
      DOC "Matlab main program"
      NO_DEFAULT_PATH
    )
  endif()
  if(MATLAB_PROGRAM)
    set(MATLAB_MAIN_PROGRAM_FOUND TRUE)
  endif()

endif()  


# Component MX library
list(FIND MATLAB_FIND_COMPONENTS MX_LIBRARY _matlab_find_mx)
if(_matlab_find_mx GREATER -1)
  
  find_library(
    MATLAB_MX_LIBRARY
    mx
    PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
    NO_DEFAULT_PATH
  )
  
  if(MATLAB_MX_LIBRARY)
    set(MATLAB_MX_LIBRARY_FOUND TRUE)
  endif()
endif()

# Component ENG library
list(FIND MATLAB_FIND_COMPONENTS ENG_LIBRARY _matlab_find_eng)
if(_matlab_find_eng GREATER -1)
  find_library(
    MATLAB_ENG_LIBRARY
    eng
    PATHS ${MATLAB_LIB_DIR_FOR_LOOKUP}
    NO_DEFAULT_PATH
  )
  if(MATLAB_ENG_LIBRARY)
    set(MATLAB_ENG_LIBRARY_FOUND TRUE)
  endif()
endif()


unset(_matlab_find_matlab_program)
unset(_matlab_find_mex_compiler)
unset(_matlab_find_mx)
unset(_matlab_find_eng)


set(MATLAB_LIBRARIES ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY})

if(CMAKE_VERSION VERSION_LESS "2.8.11")
  find_package_handle_standard_args(
    MATLAB
    REQUIRED_VARS ${MATLAB_REQUIRED_VARIABLES} MATLAB_MEX_LIBRARY
    VERSION_VAR MATLAB_VERSION)
else()
  find_package_handle_standard_args(
    MATLAB 
    FOUND_VAR MATLAB_FOUND
    REQUIRED_VARS ${MATLAB_REQUIRED_VARIABLES} MATLAB_MEX_LIBRARY #MATLAB_REQUIRED_PROGRAMS MATLAB_REQUIRED_LIBRARIES MATLAB_REQUIRED_INCLUDE_DIRS
    VERSION_VAR MATLAB_VERSION
    HANDLE_COMPONENTS)
endif()




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
    MATLAB_PROGRAM
    MATLAB_MEX_EXTENSION
  )
endif()






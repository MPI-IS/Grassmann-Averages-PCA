# Usage: cmake 
#   -Dtest_timeout=180 
#   -Dworking_directory="." 
#   -Doutput_directory=
#   -Dadditional_paths=""
#   -DMatlab_PROGRAM=matlab_exe_location 
#   -DMatlab_ADDITIONNAL_STARTUP_OPTIONS=""
#   -Dtest_name=name_of_the_test 
#   -Dunittest_file_to_run
#   -P FindMatlab_TestsRedirect.cmake

set(Matlab_UNIT_TESTS_CMD -nosplash -nojvm -nodesktop -nodisplay ${Matlab_ADDITIONNAL_STARTUP_OPTIONS})
if(WIN32)
  set(Matlab_UNIT_TESTS_CMD ${Matlab_UNIT_TESTS_CMD} -wait)
endif()

if(NOT test_timeout)
  set(test_timeout 180)
endif()

get_filename_component(unittest_file_directory   ${unittest_file_to_run} DIRECTORY)
get_filename_component(unittest_file_to_run_name ${unittest_file_to_run} NAME_WE)

set(concat_string '${unittest_file_directory}')
foreach(s IN LISTS additional_paths)
  if(NOT "${s}" STREQUAL "")
    set(concat_string "${concat_string}, '${s}'")
  endif()
endforeach()

set(Matlab_SCRIPT_TO_RUN
    "addpath(${concat_string})\; path, runtests('${unittest_file_to_run_name}'), exit(max([ans(1,:).Failed]))"
   )

set(Matlab_LOG_FILE ${output_directory}/${test_name}.log)


message(STATUS "FULL COMMAND IS ${Matlab_PROGRAM} 
      ${Matlab_UNIT_TESTS_CMD}
      -logfile ${Matlab_LOG_FILE}
      -r ${Matlab_SCRIPT_TO_RUN}")

execute_process(
  COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${Matlab_LOG_FILE} -r ${Matlab_SCRIPT_TO_RUN}
  RESULT_VARIABLE res
  TIMEOUT ${test_timeout}
  OUTPUT_QUIET # we do not want the output twice
  )


if(NOT (res EQUAL 0))
  file(READ ${Matlab_LOG_FILE} matlab_log_content)
  message("Matlab test ${name_of_the_test} FAILED\n${matlab_log_content}") # if we put FATAL_ERROR here, the file is indented.
  message( FATAL_ERROR "TEST FAILED" )
endif()



if(FALSE)
  add_test(NAME ${PROJECT_NAME}_matlabtest-1
           COMMAND ${Matlab_PROGRAM} ${Matlab_UNIT_TESTS_CMD} -logfile ${CMAKE_BINARY_DIR}/Matlab/${PROJECT_NAME}_matlabtest_1.log 
           -r "addpath('$<TARGET_FILE_DIR:${GAPCA_MEX_Project}>', '${CMAKE_SOURCE_DIR}/test'); path, runtests('grassmannpca_matlab_unit_tests'), exit(max([ans(1,:).Failed]))")
  set_tests_properties(${PROJECT_NAME}_matlabtest-1 PROPERTIES TIMEOUT 180)
endif()
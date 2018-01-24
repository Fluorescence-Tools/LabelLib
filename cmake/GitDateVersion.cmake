function (gitDateVersionString _verStr)
  if(NOT GIT_FOUND)
    find_package(Git QUIET)
  endif()
  
  if(NOT GIT_FOUND)
    set(${_verStr} "GIT-NOTFOUND" PARENT_SCOPE)
    return()
  endif()
  
  execute_process(COMMAND
		"${GIT_EXECUTABLE}" log -1 --format=%h
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE res
		OUTPUT_VARIABLE hash
		OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT res EQUAL 0)
    set(hash "nohash")
  endif()
  
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" branch --contains ${hash}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE res
    OUTPUT_VARIABLE branch
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT res EQUAL 0)
    set(branch "nobranch")
  endif()
  STRING(REPLACE "* " "" branch ${branch})
  STRING(REGEX REPLACE "\n.*" "" branch ${branch})
  
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" show -s --format=%cd --date=short ${hash}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE res
    OUTPUT_VARIABLE date
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT res EQUAL 0)
    set(date "nodate")
  endif()
  STRING(REPLACE "-" "." date ${date})
    
  SET(${_verStr} "${date}-${branch}-${hash}" PARENT_SCOPE)
endfunction(gitDateVersionString)
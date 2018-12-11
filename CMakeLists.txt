cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_CXX_STANDARD 11)
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++11")
endif(MSVC)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(Git QUIET)
set(GIT_VERSION_STRING "version-unknown")
if(GIT_FOUND)
    include(GitDateVersion)
    gitDateVersionString(GIT_VERSION_STRING)
    message(STATUS "GIT version: ${GIT_VERSION_STRING}")
endif()
add_definitions("-DGIT_VERSION_STRING=\"${GIT_VERSION_STRING}\"")

enable_testing() # needed on top-level CMakeLists.txt

# Target doc definition
# The sub projects must all define a custom target but with different names to avoid
# conflicts.
# We define as well a function AddDocTarget in cmake/DoxygenTarget.cmake
add_custom_target(doc DEPENDS docFlexLabel)
# ADD_DEPENDENCIES(doc docFlexLabel)

add_subdirectory(FlexLabel)

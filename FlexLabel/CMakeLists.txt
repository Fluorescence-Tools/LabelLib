cmake_minimum_required(VERSION 3.15)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake) # main (top) cmake dir
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

project("FlexLabel")

include(Common)

find_package(Eigen3 3.2 QUIET NO_MODULE)
if (NOT TARGET Eigen3::Eigen)
        find_path(EIGEN3_INCLUDE_DIR Eigen/Core
                PATHS
                "$ENV{CONDA_PREFIX}/include/eigen3"
                "$ENV{PREFIX}/include/eigen3"
                PATH_SUFFIXES eigen3 include/eigen3
        )
        if (EIGEN3_INCLUDE_DIR)
                message(STATUS "Using system Eigen include path: ${EIGEN3_INCLUDE_DIR}")
        elseif (EXISTS "${CMAKE_SOURCE_DIR}/thirdparty/eigen/Eigen/src/Core/util/Macros.h")
                message(STATUS "Using thirdparty/Eigen")
                set(EIGEN3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/eigen)
        else()
                message(FATAL_ERROR
                        "Eigen3 was not found and thirdparty/eigen is unavailable. "
                        "Install Eigen3 (for example: conda install eigen or apt install libeigen3-dev) "
                        "or clone with --recursive.")
        endif()
endif (NOT TARGET Eigen3::Eigen)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(pybind11 2.2 QUIET)
if(NOT pybind11_FOUND)
    if (EXISTS "${CMAKE_SOURCE_DIR}/thirdparty/pybind11/tools/pybind11Tools.cmake")
        message(STATUS "Using thirdparty/pybind11")
        add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/pybind11 pybind11_bin)
    endif()
    if (NOT TARGET pybind11::module)
        message(FATAL_ERROR
                "pybind11 was not found and thirdparty/pybind11 is unavailable. "
                "Install pybind11 (for example: pip install pybind11 or conda install pybind11) "
                "or clone with --recursive.")
    endif()
endif(NOT pybind11_FOUND)

#PCG Random Number Generation
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/pcg-cpp/include)

#Halton quasi-random numbers
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/Halton)

#Spline approximation
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/SplineApprox)

#Flat map
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/pubby-flat/include)

# This sets the include directory for the reference project. This is the -I flag in gcc.
# All the includes should be in this variable To add new directories to the include path,
# just add them in this list (after a space), or add an include_directories statement in a
# subdirectory's CMakeLists.
include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

# Mac needed variables (adapt according to your needs)
set(CMAKE_MACOSX_RPATH ON)
#set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

enable_testing()
add_subdirectory(src)
# add_subdirectory(doc)
if(NOT PYTHON_LIBRARY_OUTPUT_DIRECTORY)
  add_subdirectory(test)
endif()
add_subdirectory(python)

#
# SYSTEM PACKAGING (RPM, TGZ, ...)
# _____________________________________________________________________________

include(CPackConfig)

#
# CMAKE PACKAGING (for other CMake projects to use this one easily)
# _____________________________________________________________________________

include(PackageConfigurator)

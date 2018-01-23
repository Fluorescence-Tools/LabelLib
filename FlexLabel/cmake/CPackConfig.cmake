# General CPack configuration 
# Info: http://www.itk.org/Wiki/CMake:Component_Install_With_CPack

set(CPACK_PACKAGE_NAME "FlexLabel")

if(APPLE)
  set(CPACK_PACKAGE_VENDOR "Fluorescence-Tools") # PackageMaker doesn't like http://
else()
  set(CPACK_PACKAGE_VENDOR "https://github.com/Fluorescence-Tools") # deb lintian insists on URL
endif()

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Library for coarse-grained simulations of probes flexibly coupled to biomolecules")
set(CPACK_PACKAGE_CONTACT "Mykola Dimura <dimura@hhu.de>")
set(CPACK_PACKAGE_VERSION ${GIT_VERSION_STRING})
# set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
# set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
# set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_RPM_PACKAGE_DEBUG 0)

# Select package generator
if(MSVC)
  set(CPACK_GENERATOR "NSIS")
endif(MSVC)

if(APPLE)
  set(CPACK_GENERATOR "PackageMaker")
  set(CPACK_SET_DESTDIR ON)
endif(APPLE)

if (${CMAKE_SYSTEM_NAME} MATCHES Linux)
  find_program(RPM_EXE rpmbuild)
  if(${RPM_EXE} MATCHES RPM_EXE-NOTFOUND)
    set(CPACK_GENERATOR "TGZ;DEB")
  else()
    set(CPACK_GENERATOR "TGZ;DEB;RPM")
  endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES Linux)

# Components
# See www.cmake.org/Wiki/CMake:Component_Install_With_CPack#Controlling_Differents_Ways_of_packaging_components
# _____________________________________________________________________________
set(CPACK_COMPONENT_INSTALL ON)
set(CPACK_COMPONENTS_ALL libs devel doc)

set(CPACK_COMPONENT_LIBS_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_LIBS_DESCRIPTION "Runtime libraries.")
set(CPACK_COMPONENT_LIBS_GROUP "Runtime")
# express component dependencies this way, it will translate into package dependencies where applicable
set(CPACK_COMPONENT_LIBS_DEPENDS doc)

set(CPACK_COMPONENT_DOCS_DISPLAY_NAME "Documents")
set(CPACK_COMPONENT_DOCS_DESCRIPTION "User Documentation")
set(CPACK_COMPONENT_DOCS_GROUP "Documentation")

set(CPACK_COMPONENT_DEV_DISPLAY_NAME "Development files")
set(CPACK_COMPONENT_DEV_DESCRIPTION "Development header files and libraries, as well as cmake config files.")
set(CPACK_COMPONENT_DEV_GROUP "Development")
set(CPACK_COMPONENT_DEV_DEPENDS libs)

# Debian specific configuration (minimum)
# _____________________________________________________________________________

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${CPACK_PACKAGE_CONTACT}")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.3.1-6)")
SET(CPACK_DEBIAN_PACKAGE_CONFLICTS "")

# RPM specific configuration (minimum)
# _____________________________________________________________________________
set(CPACK_RPM_PACKAGE_LICENSE "LGPL")
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
set(CPACK_RPM_PACKAGE_VERSION ${GIT_VERSION_STRING})
set(CPACK_RPM_COMPONENT_INSTALL ON)   # necessary even if CPACK_COMPONENT_INSTALL set to ON.

# OS X PackageMaker
# _____________________________________________________________________________
set(CPACK_OSX_PACKAGE_VERSION "10.5")

include (InstallRequiredSystemLibraries)
include (CPack)

set(LIB_HEADERS
  ../include/FlexLabel/FlexLabel.h # Normal header
)

set(LIB_SOURCES
  FlexLabel.cxx
)

set(APP_HEADERS
)

set(APP_SOURCES
)

# Define two variables in order not to repeat ourselves.
set(LIBRARY_NAME FlexLabel)

# Define the library
add_library(${LIBRARY_NAME} SHARED ${LIB_SOURCES})
add_library(${LIBRARY_NAME}_static STATIC ${LIB_SOURCES})
set_target_properties(${LIBRARY_NAME}_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set the build version. It will be used in the name of the lib, with corresponding
# symlinks created. SOVERSION could also be specified for api version. 
set_target_properties(${LIBRARY_NAME} PROPERTIES
  VERSION ${GIT_VERSION_STRING}
  FRAMEWORK FALSE
  PUBLIC_HEADER "${LIB_HEADERS}" # specify the public headers
)

# Says how and where to install software
# Targets:
#   * <prefix>/lib/<libraries>
#   * header location after install: <prefix>/include/<project>/*.h
#   * headers can be included by C++ code `#<project>/Bar.hpp>`
install(TARGETS ${LIBRARY_NAME} ${APPLICATION_NAME}
  EXPORT ${PROJECT_NAME}Targets            # for downstream dependencies
  ARCHIVE DESTINATION lib COMPONENT libs   # static lib
  RUNTIME DESTINATION bin COMPONENT libs   # binaries
  LIBRARY DESTINATION lib COMPONENT libs   # shared lib
  FRAMEWORK DESTINATION bin COMPONENT libs # for mac
  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME} COMPONENT devel   # headers for mac (note the different component -> different package)
  INCLUDES DESTINATION include             # headers
)


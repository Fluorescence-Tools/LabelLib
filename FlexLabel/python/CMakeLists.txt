if(PYTHON_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYTHON_LIBRARY_OUTPUT_DIRECTORY}")
endif()

pybind11_add_module(LabelLib bindings.cpp)
target_link_libraries(LabelLib PRIVATE FlexLabel_static)

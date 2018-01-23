set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PYTHON_LIBRARY_OUTPUT_DIRECTORY}")

pybind11_add_module(LabelLib bindings.cpp)
target_link_libraries(LabelLib PRIVATE FlexLabel_static)

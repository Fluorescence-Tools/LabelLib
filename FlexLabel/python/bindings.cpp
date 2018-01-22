#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "FlexLabel/FlexLabel.h"

namespace py = pybind11;

PYBIND11_MODULE(LabelLib, m) {
    m.doc() = "LabelLib bindings"; // optional module docstring
    m.def("minLinkerLength", &minLinkerLength, "Returns single source shortest path lengths, given the obstacles.");
    
    py::class_<Grid3D>(m,"Grid3D")
    .def_readonly("discStep", &Grid3D::discStep)
    .def_readonly("originXYZ", &Grid3D::originXYZ)
    .def_readonly("shape", &Grid3D::shape)
    .def_readonly("grid", &Grid3D::grid);
}

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "FlexLabel/FlexLabel.h"

namespace py = pybind11;

Grid3D dyeDensityAV1(const Eigen::Matrix4Xf &atomsXyzr,
		     const Eigen::Vector3f &sourceXyz, const float linkerLength,
		     const float linkerDiameter, const float dyeRadius,
		     const float discStep)
{
	return dyeDensity(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			  dyeRadius, discStep);
}

Grid3D dyeDensityAV3(const Eigen::Matrix4Xf &atomsXyzr,
		     const Eigen::Vector3f &sourceXyz, const float linkerLength,
		     const float linkerDiameter,
		     const Eigen::Vector3f &dyeRadii, const float discStep)
{
	return dyeDensity(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			  dyeRadii, discStep);
}

PYBIND11_MODULE(LabelLib, m)
{
	m.doc() = "LabelLib bindings"; // optional module docstring
	m.def("minLinkerLength", &minLinkerLength,
	      "Returns single source shortest path lengths, given the obstacles.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadii"), py::arg("discStep"));
	m.def("dyeDensityAV1", &dyeDensityAV1,
	      "Returns an Accessible Volume assuming a spherical dye.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadius"), py::arg("discStep"));
	m.def("dyeDensityAV3", &dyeDensityAV3,
	      "Returns an Accessible Volume assuming a quasi-ellipsoidal dye.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadii"), py::arg("discStep"));

	py::class_<Grid3D>(m, "Grid3D")
		.def_readonly("discStep", &Grid3D::discStep)
		.def_readonly("originXYZ", &Grid3D::originXYZ)
		.def_readonly("shape", &Grid3D::shape)
		.def_readonly("grid", &Grid3D::grid);
}

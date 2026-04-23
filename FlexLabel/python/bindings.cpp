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
	m.doc() = "LabelLib C++ core bindings for accessible-volume modeling.";
	m.def("minLinkerLength", &minLinkerLength,
	      "Compute shortest linker-path lengths on a 3D grid with steric obstacles.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadius"), py::arg("discStep"));
	m.def("dyeDensityAV1", &dyeDensityAV1,
	      "Compute accessible-volume density for a spherical probe.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadius"), py::arg("discStep"));
	m.def("dyeDensityAV3", &dyeDensityAV3,
	      "Compute accessible-volume density for a quasi-ellipsoidal probe.",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadii"), py::arg("discStep"));
	m.def("addWeights", &addWeights,
	      "Add local density weights to a Grid3D.", py::arg("grid"),
	      py::arg("xyzRQ"));
	m.def("meanDistance", &meanDistance,
	      "Estimate mean inter-grid distance by stochastic sampling.",
	      py::arg("g1"), py::arg("g2"), py::arg("nsamples"));
	m.def("meanEfficiency", &meanEfficiency,
	      "Estimate mean FRET efficiency by stochastic sampling.", py::arg("g1"), py::arg("g2"),
	      py::arg("R0"), py::arg("nsamples"));
	m.def("sampleDistanceDistInv", &sampleDistanceDistInv,
		"Sample inter-grid distances with inverse transform sampling.", py::arg("g1"), py::arg("g2"),
		py::arg("nsamples"));

	py::class_<Grid3D>(m, "Grid3D",
		      "3D scalar grid used for path-length and density representations.")
	        .def(py::init<const std::array<uint32_t, 3> &,
	                      const std::array<float, 3> &, float>(),
		     py::arg("shape"), py::arg("originXYZ"), py::arg("discStep"))
	        .def_readwrite("discStep", &Grid3D::discStep)
	        .def_readwrite("originXYZ", &Grid3D::originXYZ)
	        .def_readwrite("shape", &Grid3D::shape)
	        .def_readwrite("grid", &Grid3D::grid)
	        .def("points", &Grid3D::points,
		     "Return non-zero grid points as [x, y, z, weight] columns.");
}

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "FlexLabel/FlexLabel.h"

#include <algorithm>
#include <string>
#include <stdexcept>

namespace py = pybind11;
using Matrix5Xf = Eigen::Matrix<float, 5, Eigen::Dynamic>;

Eigen::Matrix4Xf asMatrix4Xf(const py::array_t<float, py::array::c_style | py::array::forcecast> &arr,
                             const char *name)
{
    auto arrView = arr;
    py::buffer_info info = arrView.request();
    if (info.ndim != 2) {
        throw std::invalid_argument(std::string(name) + " must have 2 dimensions");
    }

    const py::ssize_t d0 = info.shape[0];
    const py::ssize_t d1 = info.shape[1];
    const float *ptr = static_cast<const float *>(info.ptr);

    if (d0 == 4) {
        Eigen::Matrix4Xf out(4, d1);
        for (py::ssize_t i = 0; i < d0; ++i) {
            for (py::ssize_t j = 0; j < d1; ++j) {
                out(i, j) = ptr[i * d1 + j];
            }
        }
        return out;
    }

    if (d1 == 4) {
        Eigen::Matrix4Xf out(4, d0);
        for (py::ssize_t i = 0; i < d0; ++i) {
            for (py::ssize_t j = 0; j < d1; ++j) {
                out(j, i) = ptr[i * d1 + j];
            }
        }
        return out;
    }

    throw std::invalid_argument(std::string(name) + " must be shaped (4, N) or (N, 4)");
}

Matrix5Xf asMatrix5Xf(const py::array_t<float, py::array::c_style | py::array::forcecast> &arr,
                      const char *name)
{
    auto arrView = arr;
    py::buffer_info info = arrView.request();
    if (info.ndim != 2) {
        throw std::invalid_argument(std::string(name) + " must have 2 dimensions");
    }

    const py::ssize_t d0 = info.shape[0];
    const py::ssize_t d1 = info.shape[1];
    const float *ptr = static_cast<const float *>(info.ptr);

    if (d0 == 5) {
        Matrix5Xf out(5, d1);
        for (py::ssize_t i = 0; i < d0; ++i) {
            for (py::ssize_t j = 0; j < d1; ++j) {
                out(i, j) = ptr[i * d1 + j];
            }
        }
        return out;
    }

    if (d1 == 5) {
        Matrix5Xf out(5, d0);
        for (py::ssize_t i = 0; i < d0; ++i) {
            for (py::ssize_t j = 0; j < d1; ++j) {
                out(j, i) = ptr[i * d1 + j];
            }
        }
        return out;
    }

    throw std::invalid_argument(std::string(name) + " must be shaped (5, N) or (N, 5)");
}

Eigen::Vector3f asVector3f(const py::array_t<float, py::array::c_style | py::array::forcecast> &arr,
                           const char *name)
{
    auto arrView = arr;
    py::buffer_info info = arrView.request();
    const float *ptr = static_cast<const float *>(info.ptr);

    if (info.ndim == 1 && info.shape[0] == 3) {
        return Eigen::Vector3f(ptr[0], ptr[1], ptr[2]);
    }

    if (info.ndim == 2 && info.shape[0] == 1 && info.shape[1] == 3) {
        return Eigen::Vector3f(ptr[0], ptr[1], ptr[2]);
    }

    if (info.ndim == 2 && info.shape[0] == 3 && info.shape[1] == 1) {
        return Eigen::Vector3f(ptr[0], ptr[1], ptr[2]);
    }

    throw std::invalid_argument(std::string(name) + " must be shaped (3,), (1, 3), or (3, 1)");
}

py::array_t<float> asArray1f(const std::vector<float> &values)
{
    py::array_t<float> out(values.size());
    py::buffer_info info = out.request();
    float *ptr = static_cast<float *>(info.ptr);
    std::copy(values.begin(), values.end(), ptr);
    return out;
}

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

Grid3D minLinkerLengthArr(const py::array_t<float, py::array::c_style | py::array::forcecast> &atomsXyzr,
                          const py::array_t<float, py::array::c_style | py::array::forcecast> &sourceXyz,
                          const float linkerLength, const float linkerDiameter,
                          const float dyeRadius, const float discStep)
{
    return minLinkerLength_arr(asMatrix4Xf(atomsXyzr, "atomsXyzr"),
                               asVector3f(sourceXyz, "sourceXyz"), linkerLength,
                               linkerDiameter, dyeRadius, discStep);
}

Grid3D dyeDensityAV1Arr(const py::array_t<float, py::array::c_style | py::array::forcecast> &atomsXyzr,
                        const py::array_t<float, py::array::c_style | py::array::forcecast> &sourceXyz,
                        const float linkerLength, const float linkerDiameter,
                        const float dyeRadius, const float discStep)
{
    return dyeDensity_arr(asMatrix4Xf(atomsXyzr, "atomsXyzr"),
                          asVector3f(sourceXyz, "sourceXyz"), linkerLength,
                          linkerDiameter, dyeRadius, discStep);
}

Grid3D dyeDensityAV3Arr(const py::array_t<float, py::array::c_style | py::array::forcecast> &atomsXyzr,
                        const py::array_t<float, py::array::c_style | py::array::forcecast> &sourceXyz,
                        const float linkerLength, const float linkerDiameter,
                        const py::array_t<float, py::array::c_style | py::array::forcecast> &dyeRadii,
                        const float discStep)
{
    return dyeDensity_arr(asMatrix4Xf(atomsXyzr, "atomsXyzr"),
                          asVector3f(sourceXyz, "sourceXyz"), linkerLength,
                          linkerDiameter, asVector3f(dyeRadii, "dyeRadii"), discStep);
}

Grid3D addWeightsArr(const Grid3D &grid,
                     const py::array_t<float, py::array::c_style | py::array::forcecast> &xyzRQ)
{
    return addWeights_arr(grid, asMatrix5Xf(xyzRQ, "xyzRQ"));
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
	      py::arg("g1"), py::arg("g2"), py::arg("nsamples") = 100000);
	m.def("meanEfficiency", &meanEfficiency,
	      "Estimate mean FRET efficiency by stochastic sampling.", py::arg("g1"), py::arg("g2"),
	      py::arg("R0"), py::arg("nsamples") = 100000);
	m.def("sampleDistanceDistInv", &sampleDistanceDistInv,
		"Sample inter-grid distances with inverse transform sampling.", py::arg("g1"), py::arg("g2"),
		py::arg("nsamples") = 1000000);
	m.def("minLinkerLength_arr", &minLinkerLengthArr,
	      "Array API: accepts atomsXyzr shape (4, N) or (N, 4).",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadius"), py::arg("discStep"));
	m.def("dyeDensityAV1_arr", &dyeDensityAV1Arr,
	      "Array API: accepts atomsXyzr shape (4, N) or (N, 4).",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadius"), py::arg("discStep"));
	m.def("dyeDensityAV3_arr", &dyeDensityAV3Arr,
	      "Array API: accepts atomsXyzr shape (4, N) or (N, 4) and dyeRadii shape (3,).",
	      py::arg("atomsXyzr"), py::arg("sourceXyz"),
	      py::arg("linkerLength"), py::arg("linkerDiameter"),
	      py::arg("dyeRadii"), py::arg("discStep"));
	m.def("addWeights_arr", &addWeightsArr,
	      "Array API: accepts xyzRQ shape (5, N) or (N, 5).",
	      py::arg("grid"), py::arg("xyzRQ"));
	m.def("meanDistance_arr", &meanDistance_arr,
	      "Array API alias for meanDistance.",
	      py::arg("g1"), py::arg("g2"), py::arg("nsamples") = 100000);
	m.def("meanEfficiency_arr", &meanEfficiency_arr,
	      "Array API alias for meanEfficiency.",
	      py::arg("g1"), py::arg("g2"), py::arg("R0"), py::arg("nsamples") = 100000);
	m.def("sampleDistanceDistInv_arr",
	      [](const Grid3D &g1, const Grid3D &g2, const unsigned nsamples) {
	          return asArray1f(sampleDistanceDistInv_arr(g1, g2, nsamples));
	      },
	      "Array API alias for sampled inter-grid distances, returned as numpy.ndarray.",
	      py::arg("g1"), py::arg("g2"), py::arg("nsamples") = 1000000);

	py::class_<Grid3D>(m, "Grid3D",
		      "3D scalar grid used for path-length and density representations.")
	        .def(py::init<const std::array<uint32_t, 3> &,
	                      const std::array<float, 3> &, float>(),
		     py::arg("shape"), py::arg("originXYZ"), py::arg("discStep"))
	        .def_readwrite("discStep", &Grid3D::discStep)
	        .def_readwrite("originXYZ", &Grid3D::originXYZ)
	        .def_readwrite("shape", &Grid3D::shape)
	        .def_readwrite("grid", &Grid3D::grid)
	        .def("origin_arr", [](const Grid3D &g) {
	            py::array_t<float> out(3);
	            py::buffer_info info = out.request();
	            float *ptr = static_cast<float *>(info.ptr);
	            ptr[0] = g.originXYZ[0];
	            ptr[1] = g.originXYZ[1];
	            ptr[2] = g.originXYZ[2];
	            return out;
	         }, "Return origin as numpy array of shape (3,).")
	        .def("shape_arr", [](const Grid3D &g) {
	            py::array_t<uint32_t> out(3);
	            py::buffer_info info = out.request();
	            auto *ptr = static_cast<uint32_t *>(info.ptr);
	            ptr[0] = g.shape[0];
	            ptr[1] = g.shape[1];
	            ptr[2] = g.shape[2];
	            return out;
	         }, "Return grid shape as numpy array of shape (3,).")
	        .def("grid_arr", [](const Grid3D &g) {
	            return asArray1f(g.grid);
	         }, "Return flattened grid values as numpy array.")
	        .def("points", &Grid3D::points,
		     "Return non-zero grid points as [x, y, z, weight] columns.");
}

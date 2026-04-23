#ifndef LABELLIB_FLEXLABEL_H
#define LABELLIB_FLEXLABEL_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector4f)

#include <vector>
#include <cstdint>
#include <array>
#include <numeric>
#include <ciso646>

/// @brief   3-dimensional grid of float numbers.
/// Can be used for density maps or other maps of scalar variables
class Grid3D
{
public:
	/// @brief   space discretization step
	float discStep;

	/// @brief   cartesian coordinates of the grid origin (grid[0,0,0]) with
	/// respect to the origin of the atoms which were used to calculate the
	/// Grid
	std::array<float, 3> originXYZ;

	/// @brief   Number of grid points in each (X, Y and Z) direction.
	/// {xMax, yMax, zMax}
	std::array<uint32_t, 3> shape;

	/// @brief   1D array of the 3D grid
	std::vector<float> grid;

	/// constructors
	Grid3D(const Grid3D &) = default;
	Grid3D &operator=(const Grid3D &) = default;
	Grid3D(Grid3D &&other) = default;

	Grid3D(const std::array<uint32_t, 3> &shape,
	       const std::array<float, 3> &originXYZ, float discStep)
	    : shape(shape), originXYZ(originXYZ), discStep(discStep)
	{
		grid.resize(shape[0] * shape[1] * shape[2], 0.0f);
	}

	Grid3D(std::array<float, 3> originXYZ, const float edgeL,
	       const float discStep, const float value)
	    : originXYZ(originXYZ), discStep(discStep)
	{
		uint32_t iL = static_cast<uint32_t>(edgeL / discStep + 0.5f);
		shape.fill(iL);
		grid.resize(iL * iL * iL, value);
	}

	/// @brief The scalar value of a certain grid point. The @param ix ,
	/// @param iy , and @param iz refer to the grid point
	float value(int ix, int iy, int iz) const
	{
		return grid[index1D(ix, iy, iz)];
	}

	/// @brief Converts the parameter @param ix @param iy and @param iz to
	/// an the index of the 1D value array. (3D->1D)
	int index1D(int ix, int iy, int iz) const
	{
		/*assert(ix >= 0 && ix < shape[0]);
		assert(iy >= 0 && iy < shape[1]);
		assert(iz >= 0 && iz < shape[2]);*/
		return ix + shape[0] * (iy + iz * shape[1]);
	}

	/// @brief Converts the parameter @param i to an the an 3D array
	/// that accesses the corresponding value of the 3D array (1D->3D).
	std::array<int, 3> index3D(int i) const
	{
		std::array<int, 3> ijk;
		ijk[2] = i / (shape[0] * shape[1]);
		i -= (ijk[2] * shape[0] * shape[1]);
		ijk[1] = i / shape[0];
		ijk[0] = i % shape[0];
		assert(ijk[0] >= 0 && ijk[0] < shape[0]);
		assert(ijk[1] >= 0 && ijk[1] < shape[1]);
		assert(ijk[2] >= 0 && ijk[2] < shape[2]);
		return ijk;
	}

	/// @brief Calculates using the input indeces @param ix, @param iy, and
	/// @param iz an 3D vector pointing to the corresponding grid point.
	std::array<float, 3> xyz(int ix, int iy, int iz) const
	{
		std::array<float, 3> arr = {float(ix), float(iy), float(iz)};
		for (int i = 0; i < 3; ++i) {
			arr[i] = arr[i] * discStep + originXYZ[i];
		}
		return arr;
	}

	/// @brief Converts the index inputs @param ix, @param iy, and
	/// @param iz, to a 3D vector in a real space pointing to the
	/// corresponding grid point.
	Eigen::Vector4f xyz(const Eigen::Vector4i &ijk) const
	{
		return ijk.cast<float>() * discStep
		       + Eigen::Vector4f(originXYZ[0], originXYZ[1],
					 originXYZ[2], 0.0f);
	}

	Eigen::Vector4f xyz(int i) const
	{
		std::array<int, 3> ijk = index3D(i);
		Eigen::Vector4i tmp(ijk[0], ijk[1], ijk[2], 0);
		return xyz(tmp);
	}

	/// @brief returns a list of cartesian points and corresponding
	/// weights for cells, that have density values > 0.0
	Eigen::Matrix4Xf points() const
	{
		const int gridSize = grid.size();
		Eigen::Matrix4Xf p(4, grid.size());
		int iPoint = 0;
		for (int i = 0; i < gridSize; ++i) {
			if (grid[i] > 0.0f) {
				p.col(iPoint) = xyz(i);
				p(3, iPoint) = grid[i];
				++iPoint;
			}
		}
		p.conservativeResize(4, iPoint);
		return p;
	}

	/// @brief returns a coordinate list with all coordinates having
	/// values > 0.0
	std::vector<Eigen::Vector4f> pointsVec() const
	{
		const size_t gridSize = grid.size();
		std::vector<Eigen::Vector4f> p;
		p.reserve(gridSize);
		for (size_t i = 0; i < gridSize; ++i) {
			if (grid[i] > 0.0f) {
				p.push_back(xyz(i));
				p.back()[3] = grid[i];
			}
		}
		p.shrink_to_fit();
		return p;
	}
};

/// @brief Compute shortest linker-path lengths around steric obstacles.
///
/// The linker diffusion problem is discretized on a 3D grid and solved with
/// Dijkstra's algorithm from @p sourceXyz. Cells blocked by atoms plus linker
/// thickness are excluded. The returned grid stores shortest path length from
/// the source for each accessible cell.
///
/// @param atomsXyzr Atom coordinates and radii (xyzr) in columns.
/// @param sourceXyz Attachment/source point in Cartesian coordinates.
/// @param linkerLength Maximum linker extension from source to probe center.
/// @param linkerDiameter Linker diameter used for steric exclusion.
/// @param dyeRadius Probe radius used to define effective obstacle size.
/// @param discStep Isotropic grid spacing in x, y, and z.
/// @return Grid3D with shortest path lengths for accessible cells.
Grid3D minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
		       const Eigen::Vector3f &sourceXyz,
		       const float linkerLength, const float linkerDiameter,
		       const float dyeRadius, const float discStep);

/// @brief Compute accessible-volume density for a spherical probe model.
///
/// Returns a grid where positive values indicate accessible cells. For this
/// overload, the probe is described by a single radius.
///
/// @param atomsXyzr Atom coordinates and radii (xyzr) in columns.
/// @param sourceXyz Attachment/source point in Cartesian coordinates.
/// @param linkerLength Maximum linker extension from source to probe center.
/// @param linkerDiameter Linker diameter used for steric exclusion.
/// @param dyeRadius Spherical probe radius.
/// @param discStep Isotropic grid spacing in x, y, and z.
/// @return Grid3D containing accessible-volume density values.
Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const float dyeRadius,
		  const float discStep);

/// @brief Compute accessible-volume density for a quasi-ellipsoidal probe.
///
/// This overload models probe geometry with directional radii.
///
/// @param atomsXyzr Atom coordinates and radii (xyzr) in columns.
/// @param sourceXyz Attachment/source point in Cartesian coordinates.
/// @param linkerLength Maximum linker extension from source to probe center.
/// @param linkerDiameter Linker diameter used for steric exclusion.
/// @param dyeRadii Probe radii along principal axes.
/// @param discStep Isotropic grid spacing in x, y, and z.
/// @return Grid3D containing accessible-volume density values.
Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const Eigen::Vector3f &dyeRadii,
		  const float discStep);

/// @brief Adds extra weights to the cells close to the specified points
/// @param grid (In) Original grid object
/// @param xyzRQ (In) Coordinates(xyz), "active" radius(R) and extra density(Q)
/// for each point
/// @return Returns a Grid3D with modified weights.
Grid3D addWeights(const Grid3D &grid,
		  const Eigen::Matrix<float, 5, Eigen::Dynamic> &xyzRQ);

/// @brief Estimate mean inter-grid distance from stochastic samples.
/// @param g1 First grid object (typically donor AV).
/// @param g2 Second grid object (typically acceptor AV).
/// @param nsamples Number of random samples; higher values improve precision.
/// @return Mean inter-grid distance.
double meanDistance(const Grid3D &g1, const Grid3D &g2,
		    const unsigned nsamples = 100000);

/// @brief Estimate mean FRET efficiency from stochastic distance samples.
/// @param g1 First grid object (typically donor AV).
/// @param g2 Second grid object (typically acceptor AV).
/// @param R0 Forster radius used in efficiency computation.
/// @param nsamples Number of random samples; higher values improve precision.
/// @return Mean FRET efficiency.
double meanEfficiency(const Grid3D &g1, const Grid3D &g2, const float R0,
		      const unsigned nsamples = 100000);

/// @brief Sample inter-grid distances using inverse transform sampling.
/// @param g1 First grid object (typically donor AV).
/// @param g2 Second grid object (typically acceptor AV).
/// @param nsamples Number of distances to sample.
/// @return Array of sampled distances.
std::vector<float> sampleDistanceDistInv(const Grid3D &g1, const Grid3D &g2,
					 const unsigned nsamples = 1000000);

#endif // LABELLIB_FLEXLABEL_H

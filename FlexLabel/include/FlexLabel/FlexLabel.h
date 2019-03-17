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

	/// @brief returns a coordinate list with all coordinates having
	/// values > 0.0
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

/// \todo Describe all the parameters. Add a more detailed general explanation.
/// @brief Returns single source shortest path lengths, given the obstacles.
/// Implementated using Dijkstra algorithm on a 3D grid.
/// @param atomsXyzr (In) Coordinates and collision radius for each atom in the
/// array. Number of atoms is the number of columns.
/// @param sourceXyz (In) Coordinates of the source point
/// @param discStep (In) Space discretization step.
/// The Same step is used for X, Y, and Z.
/// @return Returns a Grid3D of shortest path lengths.
Grid3D minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
		       const Eigen::Vector3f &sourceXyz,
		       const float linkerLength, const float linkerDiameter,
		       const float dyeRadius, const float discStep);

Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const float dyeRadius,
		  const float discStep);

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

/// @brief Calculate mean inter-dye distance between two accessible volumes
/// @param g1 (In) First grid object (donor)
/// @param g2 (In) Second grid object (acceptor)
/// @param nsamples (In) Number of samples to draw. Mean distance is determined
/// stochastically by choosing random pairs of points. Higher nsamples results
/// in more precise <Rda>.
/// @return Returns mean inter-dye distance (<Rda>).
double meanDistance(const Grid3D &g1, const Grid3D &g2,
		    const unsigned nsamples = 100000);

/// @brief Calculate mean FRET efficiency
/// @param g1 (In) First grid object (e.g. donor)
/// @param g2 (In) Second grid object (e.g. acceptor)
/// @param nsamples (In) Number of samples to draw. Mean distance is determined
/// stochastically by choosing random pairs of points. Higher nsamples results
/// in more precise <E>.
/// @return Returns mean FRET efficiency (<E>).
double meanEfficiency(const Grid3D &g1, const Grid3D &g2, const float R0,
		      const unsigned nsamples = 100000);

/// @brief Sample inter-grid distances. Inverse transform sampling is used
/// intrinsically to generate samples, that accurately approximate underlying
/// distribution.
/// @param g1 (In) First grid object (e.g. donor)
/// @param g2 (In) Second grid object (e.g. acceptor)
/// @param nsamples (In) Number of samples to draw. Distances will be drawn
/// stochastically by choosing random pairs of points. Higher \p nsamples
/// results in more accurate distribution.
/// @return Returns an array of sampled distances.
std::vector<float> sampleDistanceDistInv(const Grid3D &g1, const Grid3D &g2,
					 const unsigned nsamples = 1000000);
/// \todo Enable doxygen documentation
#endif // LABELLIB_FLEXLABEL_H

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

	/*!
	 * @brief returns the scalar value of the grid at a point specified by the
	 * integer grid point numbers @param ix, @param iy, and @param iz.
	 *
	 * @param ix (In) integer specifying a Grid point
	 * @param iy (In) integer specifying a Grid point
	 * @param iz (In) integer specifying a Grid point
	 * @return scalar value a the grid point (ix,iy,iz)
	 */
	float value(int ix, int iy, int iz) const
	{
		return grid[index1D(ix, iy, iz)];
	}

	/// @brief Converts the parameters @param ix @param iy and @param iz to
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

	/// @brief Calculates using the input indeces @param ix, @param iy, and iz an
	/// 3D vector pointing to the corresponding grid point.
	/*!
	 * @param ix (In) integer specifying a Grid point
	 * @param iy (In) integer specifying a Grid point
	 * @param iz (In) integer specifying a Grid point
	 * @return
	 */
	std::array<float, 3> xyz(int ix, int iy, int iz) const
	{
		std::array<float, 3> arr = {float(ix), float(iy), float(iz)};
		for (int i = 0; i < 3; ++i) {
			arr[i] = arr[i] * discStep + originXYZ[i];
		}
		return arr;
	}

    /// @brief an vector of grid indices to a cartesian 3D vector pointing to corresponding coordinates
	/*!
	 * Converts a 3D grid index vector to a cartesian vector pointing to the coordinates of the
	 * corresponding grid point.
	 *
	 * @param ijk (In) a 3D index vector
	 * @return 3D vector pointing the the coordinates of the grid point specified by @param ijk
	 */
	Eigen::Vector4f xyz(const Eigen::Vector4i &ijk) const
	{
		return ijk.cast<float>() * discStep
		       + Eigen::Vector4f(originXYZ[0], originXYZ[1],
					 originXYZ[2], 0.0f);
	}

    /// @brief Converts a grid index to a cartesian 3D vector pointing to corresponding coordinates
	/*!
	 * Converts grid index to a 3D grid index vector. The grid index vector is converted  to a cartesian
	 * vector pointing to the coordinates of the corresponding grid point.
	 *
	 * @param i (In) grid index
	 * @return 3D vector pointing the the coordinates of the grid point specified by @param ijk
	 */
	Eigen::Vector4f xyz(int i) const
	{
		std::array<int, 3> ijk = index3D(i);
		Eigen::Vector4i tmp(ijk[0], ijk[1], ijk[2], 0);
		return xyz(tmp);
	}

	/// @brief returns the coordinates of all grid points with values > 0.0
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

	/// @brief returns the coordinates of all grid points with values > 0.0
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

/*!
 * @brief Returns a @class Grid3D object containing the shortest path lengths from
 * a specified starting point considering a set of obstacles.
 *
 * Calculates for a probe with a radius of @param dyeRadius attached by a linker
 * of the length @param linkerLength at the position @param sourceXyz the shortest
 * path length from the attachment point considering a set of obstacles described
 * by the parameter @param atomsXyzr. The path length are calculated using Dijkstra
 * algorithm on a 3D grid. The shortest path length are stored an returned as
 * a @class Grid3D object with a grid spacing defined by the parameter @param discStep.
 *
 * The linker length @param linkerLength is the maximum extension of the linker
 * to the center of the probe. The attachment point @param sourceXyz is
 * a starting point specified by the by a vector of cartesian coordinates @param sourceXyz.
 * The obstacles are spheres. The cartesian coordinates and the radii of the obstacles
 * are provided by the parameter @param atomsXyzr.
 *
 * @param atomsXyzr (In) Coordinates and radii for all obstacles. The number of
 * obstacles is the number of columns.
 * @param sourceXyz (In) Coordinates of the source point
 * @param linkerLength (In) The maximal extension from the attachment point to the center of the probe.
 * @param linkerDiameter (In) The linker diameter
 * @param dyeRadius (In) The radius of the dye
 * @param discStep (In) Grid spacing in X, Y, and Z.
 * @return Returns a Grid3D with a grid spacing of @param discStep containing shortest path lengths.
 */
Grid3D minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
                       const Eigen::Vector3f &sourceXyz,
                       float linkerLength, float linkerDiameter,
                       float dyeRadius, float discStep);


/*!
 * @brief Computes the accessible volume (AV) for a probe with a single radius and returns the AV as
 * a @class Grid3D object. Accessible grid points contain positive values.
 *
 * @param atomsXyzr (In) Coordinates and radii for all obstacles. The number of
 * obstacles is the number of columns.
 * @param sourceXyz (In) Coordinates of the source point
 * @param linkerLength (In) The maximal extension from the attachment point to the center of the probe.
 * @param linkerDiameter (In) The linker diameter
 * @param dyeRadius (In) The radius of the dye
 * @param discStep (In) Grid spacing in X, Y, and Z.
 * @return Returns a Grid3D with a grid spacing of @param discStep containing the density of the probe.
 */
Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, float linkerLength,
		  float linkerDiameter, float dyeRadius,
		  float discStep);


/*!
 * @brief Adds extra weights to the cells close to the specified points
 *
 * @param grid (In) Original grid object
 * @param xyzRQ (In) Coordinates(xyz), "active" radius(R) and extra density(Q) for each point
 * @return Returns a Grid3D with modified weights.
 */
Grid3D addWeights(const Grid3D &grid,
		  const Eigen::Matrix<float, 5, Eigen::Dynamic> &xyzRQ);

/// @brief Calculate mean distance between two accessible volumes
/*!
 *
 * Calculates the mean distance, <R>, for randomly drawn samples on the @class Grid3D specified by the parameters
 * @param g1 and @param g2. The number of random distances is specified by the parameter @param nsamples. For a
 * large number of random samples the mean distance between two accessible volumes represented by the @param g1
 * and @param g2 can be approximated by <R>. The distance samples are drawn at random. Hence, by choosing
 * random pairs of points, specified by @param nsamples, the precision of <R> can be increased.
 *
 * In a FRET experiment, the @class Grid3D objects @param g1 and @param g2 typically represent donor and acceptor
 * fluorophores.
 *
 * @param g1 (In) First grid object
 * @param g2 (In) Second grid object
 * @param nsamples (In) Number random drawn samples to estimate the mean distance between @param g1 and @param g2
 * @return Returns the mean inter-grid distance (<R>).
 */
double meanDistance(const Grid3D &g1, const Grid3D &g2,
		    unsigned nsamples = 100000);

/// @brief Calculate mean FRET efficiency for two fluorophores that are represented by @class Grid3D objects.
/*!
 *
 * The mean FRET efficiency, <E>, between two fluorophores represented by @class Grid3D objects for @param nsanmples
 * randomly drawn distances between the two @class Grid3D objects specified by @param g1 and @param g2 is calculated.
 * For a large number of samples @param nsamples this corresponds to <E> of the inter-Grid distance distribution.
 * For every randomly drawn distance the FRET efficiency for that distance given the Förster radius specified
 * by @param R0 is calculated. These FRET efficiencies are averaged and returned.
 *
 * @param g1 (In) First grid object (e.g. donor)
 * @param g2 (In) Second grid object (e.g. acceptor)
 * @param R0 (In) Förster radius of the fluorophore pair represented by the @class Grid3D objects @param g1 and @param g1.
 * @param nsamples (In) Number of random samples to draw. Mean distance is determined
 * @return Returns mean FRET efficiency (<E>).
 */
double meanEfficiency(const Grid3D &g1, const Grid3D &g2, float R0,
        unsigned nsamples = 100000);

///@brief Sample the distances between two grids.
/*!
 * This function returns a vector of random distances between two @class Grid3D
 * objects specified by the parameters @param g1 and @param g2. To accurately approximate the distance distribution,
 * the distance distribution is sampled by the inverse transform sampling method. Here, the parameter @param nsamples
 * determines the number of random drawn distances. The accuraty of the sampling can be increased by
 * increasing the number of samples specified by @param nsamples.

 * @param g1 (In) First grid object (e.g. donor)
 * @param g2 (In) Second grid object (e.g. acceptor)
 * @param nsamples (In) Number of random samples.
 * @return Returns an array of sampled distances.
 * @see https://en.wikipedia.org/wiki/Inverse_transform_sampling
 */
std::vector<float> sampleDistanceDistInv(const Grid3D &g1, const Grid3D &g2,
					 unsigned nsamples = 1000000);

#endif // LABELLIB_FLEXLABEL_H

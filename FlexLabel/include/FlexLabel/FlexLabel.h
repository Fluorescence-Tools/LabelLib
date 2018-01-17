///
/// @file    FlexLabel.h
/// @author  Mykola Dimura
///

#ifndef LABELLIB_FLEXLABEL_H
#define LABELLIB_FLEXLABEL_H

#include <Eigen/Dense>

#include <vector>
#include <cstdint>
#include <array>

/// @brief   3-dimensional grid of float numbers.
/// Can be used for density maps or other maps of scalar variables
struct Grid3D {
	/// @brief   space discretization step
	float discStep;

	/// @brief   true coordinates of the grid[0,0,0] cell
	std::array<float, 3> originXYZ;

	/// @brief   Number of cels in each (X, Y and Z) direction.
	/// {xMax, yMax, zMax}
	std::array<uint32_t, 3> shape;

	/// @brief   3D array mapped to 1D
	std::vector<float> grid;

	Grid3D(const Grid3D &) = default;
	Grid3D &operator=(const Grid3D &) = default;
	Grid3D(Grid3D &&other) = default;

	Grid3D(std::array<float, 3> originXYZ, const float edgeL,
	       const float discStep, const float value)
	    : originXYZ(originXYZ), discStep(discStep)
	{
		uint32_t iL = static_cast<uint32_t>(edgeL / discStep + 0.5f);
		shape.fill(iL);
		grid.resize(iL * iL * iL, value);
	}

	int index1D(int ix, int iy, int iz) const
	{
		return ix + shape[0] * (iy + iz * shape[1]);
	}

	std::array<int, 3> index3D(int i) const
	{
		std::array<int, 3> index;
		index[2] = i / (shape[0] * shape[1]);
		i -= (index[2] * shape[0] * shape[1]);
		index[1] = i / shape[0];
		index[0] = i % shape[0];
		return index;
	}

	std::array<float, 3> xyz(int ix, int iy, int iz) const
	{
		std::array<float, 3> xyz = {ix, iy, iz};
		for (int i = 0; i < xyz.size(); ++i) {
			xyz[i] = xyz[i] * discStep + originXYZ[i];
		}
		return xyz;
	}

	float &value(int ix, int iy, int iz)
	{
		return grid[index1D(ix, iy, iz)];
	}
	const float &value(int ix, int iy, int iz) const
	{
		return grid[index1D(ix, iy, iz)];
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

/// \todo implement
Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const float dyeRadius,
		  const float discStep);
// const float contactR=0.0f, const float trappedFrac=-1.0f);

/// \todo implement
Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const Eigen::Vector3f &dyeRadii,
		  const float discStep);

/// \todo Implement python bindings
#endif // LABELLIB_FLEXLABEL_H

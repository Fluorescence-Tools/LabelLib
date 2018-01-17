///
/// @file    FlexLabel.cxx
/// @author  Mykola Dimura
///

#include "FlexLabel/FlexLabel.h"

#include <algorithm>
#include <set>
#include <queue>

struct Grid3DExt : public Grid3D {
	using edge_t = std::pair<float, int>;
	/// List of distances to the nearest neighbours (float) and
	/// corresponding offsets (int). The list is sorted by distance (small
	/// to large).
	std::vector<edge_t> neighIdx;

	// Returns a list of 1-D edge_index_offsets which matter
	// for the path length determination (Dijkstra) algorithm.
	// The shorter this list is, the faster Dijstra algorithm will run
	// at the cost of lower path length precision.
	// For example, including edges to only 6 nearest neighbours will result
	// in isopath surfaces that are cubic instead of spherical.
	std::vector<edge_t> essentialNeighbours();

	Grid3DExt(std::array<float, 3> originXYZ, const float edgeL,
		  const float discStep, const float value)
	    : Grid3D(originXYZ, edgeL, discStep, value)
	{
	}

	using Grid3D::index1D;
	int index1D(const Eigen::Vector4f &xyzR)
	{
		return index1D(int((xyzR[0] - originXYZ[0]) / discStep + 0.5f),
			       int((xyzR[1] - originXYZ[1]) / discStep + 0.5f),
			       int((xyzR[2] - originXYZ[2]) / discStep + 0.5f));
	}
	int index1D(const Eigen::Vector3f &xyzR)
	{
		return index1D(
			Eigen::Vector4f(xyzR[0], xyzR[1], xyzR[2], 0.0f));
	}
	inline bool isInvalidIndex(const int &i) const
	{
		if (i < 0)
			return true;
		if (i > grid.size())
			return true;
		return false;
	}
	std::vector<edge_t> neighbourEdges(const float maxR) const;
	void setMaxNeighbourDistance(const float maxR);

	inline void fillSphere(const Eigen::Vector4f &xyzR, const float value);
	inline void fillSpheres(const Eigen::Matrix4Xf &xyzR,
				const float extraClash, const float value);
	void fillOutsideSphere(const Eigen::Vector4f &xyzR, const float value);
	void setPathLengths(const Eigen::Vector3f &source);

	/// @brief Converts an array of spherical obstacles to a 3D clash grid.
	/// Cells inside obstacles are indicated by negative "density" values.
	static Grid3DExt fromObstacles(const Eigen::Matrix4Xf &xyzR,
				       const Eigen::Vector3f &sourceXyz,
				       const float discStep,
				       const float maxLength,
				       const float extraClash = 0.0f);
};

std::vector<Grid3DExt::edge_t> Grid3DExt::essentialNeighbours()
{
	std::vector<edge_t> fullList = neighbourEdges(3.0f * discStep);
	Eigen::VectorXf allowedDist(5);
	// good compromise, will result in 74 edges:
	allowedDist << 1.0f, 2.0f, 3.0f, 5.0f, 6.0f;
	allowedDist = allowedDist.cwiseSqrt() * discStep;
	std::vector<edge_t> edges;
	float minDiff;
	for (const edge_t &e : fullList) {
		minDiff = (allowedDist.array() - e.first).cwiseAbs().minCoeff();
		if (minDiff < discStep * 0.01f) {
			edges.push_back(e);
		}
	}
	edges.shrink_to_fit();
	return edges;
}

std::vector<Grid3DExt::edge_t> Grid3DExt::neighbourEdges(const float maxR) const
{
	std::vector<edge_t> idxs;
	const float maxRSq = maxR * maxR;
	const int iMaxR = static_cast<int>(maxR / discStep + 0.5f);
	const float discSq = discStep * discStep;
	idxs.clear();
	for (int dz = -iMaxR; dz <= iMaxR; ++dz) {
		for (int dy = -iMaxR; dy <= iMaxR; ++dy) {
			for (int dx = -iMaxR; dx <= iMaxR; ++dx) {
				float dSq =
					discSq * (dx * dx + dy * dy + dz * dz);
				if (dSq <= maxRSq) {
					int di = -index1D(dx, dy, dz);
					idxs.emplace_back(sqrt(dSq), di);
				}
			}
		}
	}
	idxs.shrink_to_fit();
	std::sort(idxs.begin(), idxs.end(),
		  [](const auto &lhs, const auto &rhs) {
			  return lhs.first < rhs.first;
		  });
	return idxs;
}

void Grid3DExt::setMaxNeighbourDistance(const float maxR)
{
	if (!neighIdx.empty()) {
		if (neighIdx.back().first >= maxR) {
			return;
		}
	}
	neighIdx = neighbourEdges(maxR);
}

inline void Grid3DExt::fillSphere(const Eigen::Vector4f &xyzR,
				  const float value)
{
	// Vector4f is faster than Vector3f, therefore we "ignore" Vector4f[3]
	Eigen::Vector4f dr = xyzR - Eigen::Vector4f(originXYZ[0], originXYZ[1],
						    originXYZ[2], 0.0f);
	const Eigen::Vector4f boxSize = {shape[0] * discStep,
					 shape[1] * discStep,
					 shape[2] * discStep, 0.0f};
	if ((dr.array() < -xyzR[3]).any()) {
		return;
	}
	if (((dr - boxSize).array() > xyzR[3]).any()) {
		return;
	}
	const int i0 = index1D(xyzR);
	int i;
	for (const auto &pair : neighIdx) {
		if (pair.first > xyzR[3]) {
			break;
		}
		i = i0 + pair.second;
		if (isInvalidIndex(i)) {
			continue;
		}
		grid[i] = value;
	}
}

void Grid3DExt::fillSpheres(const Eigen::Matrix4Xf &xyzR,
			    const float extraClash, const float value)
{
	const float maxRclash = xyzR.row(3).maxCoeff() + extraClash;
	setMaxNeighbourDistance(maxRclash);
	Eigen::Vector4f at;
	for (int iAtom = 0; iAtom < xyzR.cols(); ++iAtom) {
		at = xyzR.col(iAtom);
		at[3] += extraClash;
		fillSphere(at, value);
	}
}

void Grid3DExt::fillOutsideSphere(const Eigen::Vector4f &xyzR,
				  const float value)
{
	const int x0 =
		static_cast<int>((xyzR[0] - originXYZ[0] + 0.5f) / discStep);
	const int y0 =
		static_cast<int>((xyzR[1] - originXYZ[1] + 0.5f) / discStep);
	const int z0 =
		static_cast<int>((xyzR[2] - originXYZ[2] + 0.5f) / discStep);

	const int maxRSq =
		static_cast<int>(xyzR[3] * xyzR[3] / (discStep * discStep));

	int i = 0;
	for (int z = 0; z < shape[2]; ++z) {
		for (int y = 0; y < shape[1]; ++y) {
			for (int x = 0; x < shape[0]; ++x) {
				assert(i == index1D(x, y, z));
				int dSq = (x - x0) * (x - x0)
					  + (y - y0) * (y - y0)
					  + (z - z0) * (z - z0);
				if (dSq > maxRSq) {
					grid[i] = value;
				}
				++i;
			}
		}
	}
}

void Grid3DExt::setPathLengths(const Eigen::Vector3f &source)
{
	// perform dijkstra algorithm
	using Eigen::Vector4f;
	using std::vector;
	const vector<edge_t> &edges = essentialNeighbours();
	const int sourceVertex = index1D(source);
	grid[sourceVertex] = 0.0f;

	// it is possible to improve performance by 10-20%
	// using more efficient heap/queue implementation like radix heap
	std::priority_queue<edge_t, vector<edge_t>, std::greater<edge_t>> que;
	que.emplace(0.0f, sourceVertex);
	while (!que.empty()) {
		const edge_t qt = que.top();
		que.pop();
		if (qt.first > grid[qt.second])
			continue;
		edge_t t;
		for (const edge_t &e : edges) {
			t.second = qt.second + e.second;
			int &idx = t.second;
			if (isInvalidIndex(idx))
				continue;
			if (grid[idx] < 0.0f) {
				continue;
			}
			t.first = qt.first + e.first;
			if (t.first < grid[idx]) {
				grid[idx] = t.first;
				que.push(t);
			}
		}
	}
}

Grid3DExt Grid3DExt::fromObstacles(const Eigen::Matrix4Xf &xyzR,
				   const Eigen::Vector3f &sourceXyz,
				   const float discStep, const float maxLength,
				   const float extraClash)
{
	const float adjL = ceil(maxLength / discStep) * discStep;
	std::array<float, 3> origin = {sourceXyz[0] - adjL, sourceXyz[1] - adjL,
				       sourceXyz[2] - adjL};
	using Eigen::Vector4f;
	Grid3DExt grid(origin, adjL * 2.0f, discStep,
		       std::numeric_limits<float>::max());
	grid.fillSpheres(xyzR, extraClash, -1.0f);
	return grid;
}

Grid3D minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
		       const Eigen::Vector3f &sourceXyz,
		       const float linkerLength, const float linkerDiameter,
		       const float dyeRadius, const float discStep)
{
	Grid3DExt grid =
		Grid3DExt::fromObstacles(atomsXyzr, sourceXyz, discStep,
					 linkerLength, linkerDiameter * 0.5f);
	Eigen::Vector4f exclusionArea = {sourceXyz[0], sourceXyz[1],
					 sourceXyz[2], linkerLength};
	grid.fillOutsideSphere(exclusionArea, -2.0f);
	grid.setPathLengths(sourceXyz);
	grid.fillSpheres(atomsXyzr, dyeRadius, -3.0f);
	return grid;
}

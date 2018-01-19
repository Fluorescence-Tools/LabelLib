///
/// @file    FlexLabel.cxx
/// @author  Mykola Dimura
///

#include "FlexLabel/FlexLabel.h"

#include <algorithm>
#include <set>
#include <queue>

Eigen::Vector4f Vec4f(const std::array<float, 3> &arr)
{
	return Eigen::Vector4f(arr[0], arr[1], arr[2], 0.0f);
}
Eigen::Vector4f Vec4f(const Eigen::Vector3f &vec)
{
	return Eigen::Vector4f(vec[0], vec[1], vec[2], 0.0f);
}
Eigen::Vector4f Vec4f(const std::array<uint32_t, 3> &arr)
{
	return Eigen::Vector4f(arr[0], arr[1], arr[2], 0.0f);
}

struct Grid3DExt final : public Grid3D {
      public:
        struct edge_t {
                Eigen::Array4i ijk;
                float r;
                bool operator>(const edge_t &other) const
                {
                        return r > other.r;
                }
        };

	Grid3DExt(std::array<float, 3> originXYZ, const float edgeL,
		  const float discStep, const float value)
	    : Grid3D(originXYZ, edgeL, discStep, value)
	{
		originAdj = Vec4f(originXYZ).array() - 0.5f * discStep;
		originAdj[3] = 0.0f;
		shape4i << shape[0], shape[1], shape[2],
			std::numeric_limits<int>::max();
	}

	/// @brief Converts an array of spherical obstacles to a 3D clash grid.
	/// Cells inside obstacles are indicated by negative "density" values.
	static Grid3DExt fromObstacles(const Eigen::Matrix4Xf &xyzR,
				       const Eigen::Vector3f &sourceXyz,
				       const float discStep,
				       const float maxLength,
				       const float extraClash = 0.0f);

	inline void fillSpheres(const Eigen::Matrix4Xf &xyzR,
				const float extraClash, const float gridRef);

	void fillOutsideSphere(const Eigen::Vector4f &xyzR,
			       const float gridRef);

	void setPathLengths(const Eigen::Vector3f &source);

	void setAboveThreshold(const float threshold, const float rho);

	void excludeConcentricSpheres(const Eigen::Matrix4Xf &xyzR,
				      const Eigen::VectorXf extraClashes,
				      const float maxRho);

      private:
        Eigen::Vector4f originAdj; // originXYZ - discStep*0.5
        Eigen::Array4i shape4i;
        float devNull;
        Eigen::Array4i getIjk(const Eigen::Vector4f &xyz);
        float &gridRef(const Eigen::Vector4f &xyz);
        float &gridRef(const Eigen::Array4i &ijk);
        float gridVal(const Eigen::Array4i &ijk);

	/// List of distances to the nearest neighbours and
	/// corresponding offsets. The list is sorted by distance (small
	/// to large).
	std::vector<edge_t> neighbours;

	// The shorter this list is, the faster Dijstra algorithm will run
	// at the cost of lower path length precision.
	// For example, including edges to only 6 nearest neighbours will result
	// in isopath surfaces that are cubic instead of spherical.
	std::vector<edge_t> essentialNeighbours();


	using Grid3D::index1D;
	std::vector<edge_t> neighbourEdges(const float maxR) const;
	void setMaxNeighbourDistance(const float maxR);

	inline void fillSphere(const Eigen::Vector4f &xyzR, const float val);
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
		minDiff = (allowedDist.array() - e.r).cwiseAbs().minCoeff();
		if (minDiff < discStep * 0.01f) {
			edges.push_back(e);
		}
	}
	edges.shrink_to_fit();
	return edges;
}

float &Grid3DExt::gridRef(const Eigen::Vector4f &xyz)
{
	return gridRef(getIjk(xyz));
}

float &Grid3DExt::gridRef(const Eigen::Array4i &ijk)
{
	assert(ijk[3] == 0);
	if ((ijk < 0).any()) {
		return devNull;
	}
	if ((shape4i - ijk <= 0).any()) {
		return devNull;
	}
	return grid[index1D(ijk[0], ijk[1], ijk[2])];
}

float Grid3DExt::gridVal(const Eigen::Array4i &ijk)
{
	assert(ijk[3] == 0);
	if ((ijk < 0).any()) {
		return -1.0f;
	}
	if ((shape4i - ijk <= 0).any()) {
		return -1.0f;
	}
	return grid[index1D(ijk[0], ijk[1], ijk[2])];
}

std::vector<Grid3DExt::edge_t> Grid3DExt::neighbourEdges(const float maxR) const
{
	std::vector<edge_t> idxs;
	const float maxRSq = maxR * maxR;
	const int iMaxR = static_cast<int>(maxR / discStep + 0.5f);
	const float discSq = discStep * discStep;
	idxs.clear();
	Eigen::Array4i ijk(0, 0, 0, 0);
	int &dx = ijk.coeffRef(0);
	int &dy = ijk.coeffRef(1);
	int &dz = ijk.coeffRef(2);
	for (dz = -iMaxR; dz <= iMaxR; ++dz) {
		for (dy = -iMaxR; dy <= iMaxR; ++dy) {
			for (dx = -iMaxR; dx <= iMaxR; ++dx) {
				float dSq = discSq * ijk.matrix().squaredNorm();
				if (dSq <= maxRSq) {
					idxs.push_back({ijk, sqrtf(dSq)});
				}
			}
		}
	}
	idxs.shrink_to_fit();
	std::sort(
		idxs.begin(), idxs.end(),
		[](const auto &lhs, const auto &rhs) { return lhs.r < rhs.r; });
	return idxs;
}

void Grid3DExt::setMaxNeighbourDistance(const float maxR)
{
	if (!neighbours.empty()) {
		if (neighbours.back().r >= maxR) {
			return;
		}
	}
	neighbours = neighbourEdges(maxR);
}

void Grid3DExt::excludeConcentricSpheres(const Eigen::Matrix4Xf &xyzR,
					 const Eigen::VectorXf extraClashes,
					 const float maxRho)
{
	using Eigen::Vector4f;
	using Eigen::VectorXf;
	using Eigen::Array4i;
	using std::min;
	const float maxRclash =
		xyzR.row(3).maxCoeff() + extraClashes.maxCoeff();
	setMaxNeighbourDistance(maxRclash);

	const int numClashes = extraClashes.size();
	VectorXf effR = extraClashes;
	std::sort(effR.data(), effR.data() + numClashes);
	const VectorXf sortedR = effR;
	const VectorXf rhos = VectorXf::LinSpaced(numClashes + 1, 0.0f, maxRho);

	for (int iAtom = 0; iAtom < xyzR.cols(); ++iAtom) {
		const Vector4f &xyzRc = xyzR.col(iAtom);
		Array4i ijk0 = getIjk(xyzRc);
		effR = sortedR.array() + xyzRc.coeff(3);
		int iNei = 0;
		for (int iClash = 0; iClash < numClashes; ++iClash) {
			const float &curR = effR[iClash];
			for (; neighbours[iNei].r <= curR; ++iNei) {
				if (iNei >= neighbours.size()) {
					break;
				}
				Array4i ijk = ijk0 + neighbours[iNei].ijk;
				float &ref = gridRef(ijk);
				ref = min(ref, rhos[iClash]);
			}
		}
	}
}

Eigen::Array4i Grid3DExt::getIjk(const Eigen::Vector4f &xyz)
{
	Eigen::Array4i ijk = ((xyz - originAdj) / discStep).cast<int>();
	ijk.coeffRef(3) = 0;
	return ijk;
}

inline void Grid3DExt::fillSphere(const Eigen::Vector4f &xyzR, const float val)
{
	using Eigen::Array4i;
	const Array4i ijk0 = getIjk(xyzR);
	for (const auto &e : neighbours) {
		if (e.r > xyzR[3]) {
			break;
		}
		gridRef(Array4i(ijk0 + e.ijk)) = val;
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
	const Eigen::Array4i ijk0 = getIjk(xyzR);
	const int maxRSq = static_cast<int>(
		xyzR[3] * xyzR[3] / (discStep * discStep) + 0.5f);

	int i = 0, dSq;
	Eigen::Array4i ijk(0, 0, 0, 0);
	int &x = ijk.coeffRef(0);
	int &y = ijk.coeffRef(1);
	int &z = ijk.coeffRef(2);
	for (z = 0; z < shape[2]; ++z) {
		for (y = 0; y < shape[1]; ++y) {
			for (x = 0; x < shape[0]; ++x) {
				assert(i == index1D(x, y, z));
				dSq = (ijk - ijk0).matrix().squaredNorm();
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
	using Eigen::Array4i;
	using std::vector;
	const vector<edge_t> &edges = essentialNeighbours();
	const Array4i ijk0 = getIjk(Vec4f(source));
	gridRef(ijk0) = 0.0f;

	// it is possible to improve performance by 10-20%
	// using more efficient heap/queue implementation like radix heap
	std::priority_queue<edge_t, vector<edge_t>, std::greater<edge_t>> que;
	que.push({ijk0, 0.0f});
	while (!que.empty()) {
		const edge_t qt = que.top();
		que.pop();
		if (qt.r > gridRef(qt.ijk))
			continue;
		edge_t t;
		for (const edge_t &e : edges) {
			t.ijk = qt.ijk + e.ijk;
			const float &val = gridVal(t.ijk);
			if (val < 0.0f) {
				continue;
			}
			t.r = qt.r + e.r;
			if (t.r < val) {
				gridRef(t.ijk) = t.r;
				que.push(t);
			}
		}
	}
}

void Grid3DExt::setAboveThreshold(const float threshold, const float rho)
{
	for (float &val : grid) {
		val = val > threshold ? rho : val;
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

Grid3DExt minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
			  const Eigen::Vector3f &sourceXyz,
			  const float linkerLength, const float linkerDiameter,
			  const float discStep)
{
	Grid3DExt grid =
		Grid3DExt::fromObstacles(atomsXyzr, sourceXyz, discStep,
					 linkerLength, linkerDiameter * 0.5f);
	Eigen::Vector4f exclusionArea = {sourceXyz[0], sourceXyz[1],
					 sourceXyz[2], linkerLength};
	grid.fillOutsideSphere(exclusionArea, -2.0f);
	grid.setPathLengths(sourceXyz);
	return grid;
}

Grid3D minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
		       const Eigen::Vector3f &sourceXyz,
		       const float linkerLength, const float linkerDiameter,
		       const float dyeRadius, const float discStep)
{
	Grid3DExt grid = minLinkerLength(atomsXyzr, sourceXyz, linkerLength,
					 linkerDiameter, discStep);
	grid.fillSpheres(atomsXyzr, dyeRadius, -3.0f);
	return grid;
}

Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const Eigen::VectorXf &dyeRadii,
		  const float discStep)
{
	Grid3DExt grid = minLinkerLength(atomsXyzr, sourceXyz, linkerLength,
					 linkerDiameter, discStep);
	grid.setAboveThreshold(linkerLength, -4.0f);
	grid.setAboveThreshold(0.0f, 1.0f);
	grid.excludeConcentricSpheres(atomsXyzr, dyeRadii, 1.0f);
	return grid;
}

Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const float dyeRadius,
		  const float discStep)
{
	Eigen::VectorXf dyeRadii(1);
	dyeRadii << dyeRadius;
	return dyeDensity(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			  dyeRadii, discStep);
}

Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const Eigen::Vector3f &dyeRadii,
		  const float discStep)
{
	return dyeDensity(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			  Eigen::VectorXf(dyeRadii), discStep);
}

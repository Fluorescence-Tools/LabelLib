#include "FlexLabel/FlexLabel.h"

#include "pcg_random.hpp"

#include <map>
#include <random>
#include <queue>

Eigen::Vector4f Vec4f(const std::array<float, 3> &arr)
{
	return Eigen::Vector4f(arr[0], arr[1], arr[2], 0.0f);
}

Eigen::Vector4f Vec4f(const Eigen::Vector3f &vec, const float v3 = 0.0f)
{
	return Eigen::Vector4f(vec[0], vec[1], vec[2], v3);
}

Eigen::Vector4f Vec4f(const std::array<uint32_t, 3> &arr)
{
	return Eigen::Vector4f(arr[0], arr[1], arr[2], 0.0f);
}

std::array<float, 3> StdArr3f(const Eigen::Vector4f &xyz)
{
	return std::array<float, 3>() = {xyz[0], xyz[1], xyz[2]};
}

using Matrix5Xf = Eigen::Matrix<float, 5, Eigen::Dynamic>;

class Grid3DExt final : public Grid3D
{
public:
	struct edge_t {
		Eigen::Array4i ijk;
		float r;

		bool operator>(const edge_t &other) const
		{
			return r > other.r;
		}
	};

	struct edge1D_t {
		float r;
		int idx;

		bool operator>(const edge1D_t &other) const
		{
			return r > other.r;
		}
	};

	Grid3DExt(const Eigen::Vector4f &originXYZ, const float edgeL,
		  const float discStep, const float value)
	    : Grid3D(StdArr3f(originXYZ), edgeL, discStep, value)
	{
		originAdj = originXYZ.array() - 0.5f * discStep;
		originAdj[3] = 0.0f;
		shape4i << shape[0], shape[1], shape[2],
			std::numeric_limits<int>::max();
	}

	Grid3DExt(const Grid3D &grid) : Grid3D(grid)
	{
		originAdj = Vec4f(originXYZ).array() - 0.5f * discStep;
		originAdj[3] = 0.0f;
		shape4i << shape[0], shape[1], shape[2],
			std::numeric_limits<int>::max();
	}

	static Grid3DExt fromSource(const Eigen::Vector3f &sourceXyz,
				    const float discStep, const float maxLength,
				    const float maxClashR)
	{
		using std::max;
		const float effL = maxLength + max(maxClashR, discStep * 3.0f);
		const float adjL = ceil(effL / discStep) * discStep;
		using Eigen::Vector4f;
		const Vector4f origin = Vec4f(sourceXyz).array() - adjL;
		const float maxFloat = std::numeric_limits<float>::max();
		Grid3DExt grid(origin, adjL * 2.0f, discStep, maxFloat);
		const Vector4f exclusionArea = Vec4f(sourceXyz, maxLength);
		grid.fillOutsideSphere(exclusionArea, -1.0f);
		return grid;
	}

	inline void fillSpheres(const Eigen::Matrix4Xf &xyzR,
				const float extraClash, const float gridRef);

	void fillOutsideSphere(const Eigen::Vector4f &xyzR,
			       const float gridRef);

	void setPathLengths(const Eigen::Vector3f &source);

	void setAboveThreshold(const float threshold, const float rho);

	void excludeConcentricSpheres(const Eigen::Matrix4Xf &xyzR,
				      const Eigen::VectorXf extraClashes,
				      const float maxRho);

	void addDensity(const Matrix5Xf &xyzRQ);

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

	std::vector<edge1D_t> toEdges1D(const std::vector<edge_t> &in) const
	{
		std::vector<edge1D_t> out;
		for (const edge_t &e : in) {
			edge1D_t t;
			t.r = e.r;
			t.idx = index1D(e.ijk[0], e.ijk[1], e.ijk[2]);
			out.push_back(t);
		}
		return out;
	}

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
	std::sort(idxs.begin(), idxs.end(),
		  [](const edge_t &lhs, const edge_t &rhs) {
			  return lhs.r < rhs.r;
		  });
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
	using Eigen::Array4i;
	using Eigen::Vector4f;
	using Eigen::VectorXf;
	using std::min;
	const float maxRclash =
		xyzR.row(3).maxCoeff() + extraClashes.maxCoeff();
	setMaxNeighbourDistance(maxRclash);

	const int numClashes = extraClashes.size();
	VectorXf effR = extraClashes;
	std::sort(effR.data(), effR.data() + numClashes);
	const VectorXf sortedR = effR;
	const VectorXf rhos = VectorXf::LinSpaced(numClashes + 1, 0.0f, maxRho);

	const Vector4f center =
		Vec4f(shape) * discStep * 0.5f + Vec4f(originXYZ);
	const float maxL = shape4i.minCoeff() * discStep * 0.5f;
	const float l2out = std::pow(maxL + maxRclash, 2.0f);

	Vector4f at;
	float rAt;
	for (int iAtom = 0; iAtom < xyzR.cols(); ++iAtom) {
		at = xyzR.col(iAtom);
		rAt = at[3];
		at[3] = 0.0f;
		const float dist2 = (center - at).squaredNorm();
		if (dist2 > l2out) {
			continue;
		}
		Array4i ijk0 = getIjk(at);
		effR = sortedR.array() + rAt;
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

void Grid3DExt::addDensity(const Matrix5Xf &xyzRQ)
{
	using Eigen::Array4i;
	using Eigen::Vector4f;
	using Eigen::VectorXf;
	using std::min;
	const float maxRclash = xyzRQ.row(3).maxCoeff();
	setMaxNeighbourDistance(maxRclash + discStep);

	const Vector4f center =
		Vec4f(shape) * discStep * 0.5f + Vec4f(originXYZ);
	const float maxL = shape4i.minCoeff() * discStep * 0.5f;
	const float loutSq = std::pow(maxL + maxRclash, 2.0f);

	Vector4f xyz;
	float rAt;
	for (int iAtom = 0; iAtom < xyzRQ.cols(); ++iAtom) {
		xyz = xyzRQ.col(iAtom).head<4>();
		const float newDensity = xyzRQ(4, iAtom);
		const float curR = xyz[3];
		xyz[3] = 0.0f;
		const float dist2 = (center - xyz).squaredNorm();
		if (dist2 > loutSq) {
			continue;
		}

		const Array4i ijkAt = getIjk(xyz);
		int iNei = 0;
		for (; neighbours[iNei].r <= curR; ++iNei) {
			assert(iNei < neighbours.size());
			Array4i ijk = ijkAt + neighbours[iNei].ijk;
			float &ref = gridRef(ijk);
			if (ref > 0) {
				ref += newDensity;
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

	using Eigen::Vector4f;
	Vector4f center = Vec4f(shape) * discStep * 0.5f + Vec4f(originXYZ);
	const float maxL = shape4i.minCoeff() * discStep * 0.5f;
	const float l2out = std::pow(maxL + maxRclash, 2.0f);

	Vector4f at;
	for (int iAtom = 0; iAtom < xyzR.cols(); ++iAtom) {
		at = xyzR.col(iAtom);
		center[3] = at[3];
		const float dist2 = (center - at).squaredNorm();
		if (dist2 > l2out) {
			continue;
		}
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
	// the grid must be padded, i.e. N outermost layers of edges must be set
	// to a value < 0. Where N = the maximum depth of essentialNeighbours().
	// This can be done by fillOutsideSphere()
	using Eigen::Array4i;
	using Eigen::Vector4f;
	using std::vector;
	const vector<edge1D_t> &edges = toEdges1D(essentialNeighbours());
	const Array4i ijk0 = getIjk(Vec4f(source));
	const int i0 = index1D(ijk0[0], ijk0[1], ijk0[2]);
	grid[i0] = 0.0f;

	// it is possible to improve performance by 10-20%
	// using more efficient heap/queue implementation like radix heap
	std::priority_queue<edge1D_t, vector<edge1D_t>, std::greater<edge1D_t>>
		que;
	que.push({0.0f, i0});
	while (!que.empty()) {
		const edge1D_t qt = que.top();
		que.pop();
		if (qt.r > grid[qt.idx])
			continue;
		edge1D_t t;
		for (const edge1D_t &e : edges) {
			t.idx = qt.idx + e.idx;
			const float &val = grid[t.idx];
			if (val < 0.0f) {
				continue;
			}
			t.r = qt.r + e.r;
			if (t.r < val) {
				grid[t.idx] = t.r;
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

Grid3DExt minLinkerLength(const Eigen::Matrix4Xf &atomsXyzr,
			  const Eigen::Vector3f &sourceXyz,
			  const float linkerLength, const float linkerDiameter,
			  const float discStep)
{
	Grid3DExt grid =
		Grid3DExt::fromSource(sourceXyz, discStep, linkerLength, 0.0f);
	grid.fillSpheres(atomsXyzr, linkerDiameter * 0.5f, -1.0f);
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

Grid3DExt dyeDensityExt(const Eigen::Matrix4Xf &atomsXyzr,
			const Eigen::Vector3f &sourceXyz,
			const float linkerLength, const float linkerDiameter,
			const Eigen::VectorXf &dyeRadii, const float discStep)
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
	return dyeDensityExt(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			     dyeRadii, discStep);
}

Grid3D dyeDensity(const Eigen::Matrix4Xf &atomsXyzr,
		  const Eigen::Vector3f &sourceXyz, const float linkerLength,
		  const float linkerDiameter, const Eigen::Vector3f &dyeRadii,
		  const float discStep)
{
	return dyeDensityExt(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
			     Eigen::VectorXf(dyeRadii), discStep);
}

Grid3D addWeights(const Grid3D &grid, const Matrix5Xf &xyzRQ)
{
	Grid3DExt ext(grid);
	ext.addDensity(xyzRQ);
	return ext;
}

template <typename T> class InverseSampler
{
	// Draw elements using Inverse transform sampling
private:
	// probability is the key and index in the points array
	// is the value. Probabilities may have arbitrary normalization, i.e
	// last element in the map has maximum possible potability as the key.
	std::map<uint32_t, unsigned> map;
	T vec;
	using val_t = typename T::value_type;
	mutable pcg32_fast rng{pcg_extras::seed_seq_from<std::random_device>{}};

public:
	typename T::value_type get_random() const
	{
		// return element index, assuming that elements between the keys
		// in the map are uniformly sampled

		unsigned rnd = rng();
		// find the relevant index range (idxmin..idxmax]
		auto it = map.upper_bound(rnd);
		if (it == map.end()) {
			return vec[map.rbegin()->second];
		}
		uint32_t rmax = it->first;
		unsigned idxmax = it->second;
		--it;
		uint32_t rmin = it->first;
		unsigned idxmin = it->second;

		float relRnd = float(rnd - rmin) / (rmax - rmin);
		unsigned i = idxmin + unsigned(relRnd * (idxmax - idxmin));
		return vec[i];
	}

	template <typename KeyAccessor>
	InverseSampler(T &&_vec, KeyAccessor accessor) : vec(_vec)
	{
		// Since entries in the flat_map are only added for unique
		// weights(densities), the map stays small and fast.
		using Eigen::Vector4f;

		auto adder = [accessor](double sum, const val_t &p) {
			return sum + accessor(p);
		};
		const float sumWeights =
			std::accumulate(vec.begin(), vec.end(), 0.0, adder);
		const float pRatio = rng.max() / sumWeights;

		auto comparator = [accessor](const val_t &v1,
					     const val_t &v2) -> bool {
			return accessor(v1) > accessor(v2);
		};
		if (not std::is_sorted(vec.begin(), vec.end(), comparator)) {
			std::sort(vec.begin(), vec.end(), comparator);
		}

		map.emplace(0, 0);
		double cumSum = accessor(vec[0]) * pRatio;
		float prevSlope = cumSum;
		for (unsigned i = 1; i < vec.size(); ++i) {
			float slope = accessor(vec[i]) * pRatio;
			if (slope < prevSlope * 0.999) {
				map.emplace(cumSum, i - 1);
				prevSlope = slope;
			}
			cumSum += slope;
		}
		map.emplace(rng.max(), vec.size() - 1);
	}
};

double meanDistanceInv(const Grid3D &g1, const Grid3D &g2,
		       const unsigned nsamples)
{
	// Draw points using Inverse transform sampling

	auto el3getter = [](const Eigen::Vector4f &p) { return p[3]; };
	using points_type = std::vector<Eigen::Vector4f>;
	InverseSampler<points_type> sampler1(g1.pointsVec(), el3getter);
	InverseSampler<points_type> sampler2(g2.pointsVec(), el3getter);

	double r = 0.;
	Eigen::Vector4f tmp;
	for (unsigned s = 0; s < nsamples; s++) {
		tmp = sampler1.get_random() - sampler2.get_random();
		tmp[3] = 0.0f;
		r += tmp.norm();
	}
	return r / nsamples;
}

std::vector<float> sampleDistanceDistInv(const Grid3D &g1, const Grid3D &g2,
	const unsigned nsamples)
{
	// Draw distances using Inverse transform sampling

	std::vector<float> distances;
	distances.reserve(nsamples);

	auto el3getter = [](const Eigen::Vector4f &p) { return p[3]; };
	using points_type = std::vector<Eigen::Vector4f>;
	InverseSampler<points_type> sampler1(g1.pointsVec(), el3getter);
	InverseSampler<points_type> sampler2(g2.pointsVec(), el3getter);

	Eigen::Vector4f tmp;
	for (unsigned s = 0; s < nsamples; ++s) {
		tmp = sampler1.get_random() - sampler2.get_random();
		tmp[3] = 0.0f;
		distances.emplace_back(tmp.norm());
	}
	return distances;
}

double FRETefficiency(double R, double R0)
{
	// return 1.0 / (1.0 + std::pow(R / R0, 6));
	double ratSq = R / R0;
	ratSq *= ratSq;
	return 1.0 / (1.0 + ratSq * ratSq * ratSq);
}

double meanEfficiencyUniform(const Grid3D &g1, const Grid3D &g2, const float R0,
			     const unsigned nsamples)
{
	// Draw point indexes from uniform distribution
	pcg32 rng(pcg_extras::seed_seq_from<std::random_device>{});

	Eigen::Matrix4Xf p1 = g1.points();
	Eigen::Matrix4Xf p2 = g2.points();

	const unsigned size1 = p1.cols();
	const unsigned size2 = p2.cols();

	double totalW = 0.0;
	double e = 0., w;
	unsigned i, j;
	Eigen::Vector4f tmp;
	for (unsigned s = 0; s < nsamples; s++) {
		i = rng(size1);
		j = rng(size2);
		w = p1(3, i) * p2(3, j);
		totalW += w;
		tmp = p1.col(i) - p2.col(j);
		tmp[3] = 0.0f;
		e += w * FRETefficiency(tmp.norm(), R0);
	}
	return e / totalW;
}

double meanDistanceUniform(const Grid3D &g1, const Grid3D &g2,
			   const unsigned nsamples)
{
	// Draw point indexes from uniform distribution
	pcg32 rng(pcg_extras::seed_seq_from<std::random_device>{});

	Eigen::Matrix4Xf p1 = g1.points();
	Eigen::Matrix4Xf p2 = g2.points();

	const unsigned size1 = p1.cols();
	const unsigned size2 = p2.cols();

	double totalW = 0.0;
	double r = 0., w;
	unsigned i, j;
	Eigen::Vector4f tmp;
	for (unsigned s = 0; s < nsamples; s++) {
		i = rng(size1);
		j = rng(size2);
		w = p1(3, i) * p2(3, j);
		totalW += w;
		tmp = p1.col(i) - p2.col(j);
		tmp[3] = 0.0f;
		r += w * tmp.norm();
	}
	return r / totalW;
}

double meanEfficiencyInv(const Grid3D &g1, const Grid3D &g2, const float R0,
			 const unsigned nsamples)
{
	// Draw points using Inverse transform sampling
	auto el3getter = [](const Eigen::Vector4f &p) { return p[3]; };
	using points_type = std::vector<Eigen::Vector4f>;
	InverseSampler<points_type> sampler1(g1.pointsVec(), el3getter);
	InverseSampler<points_type> sampler2(g2.pointsVec(), el3getter);

	double e = 0.;
	Eigen::Vector4f tmp;
	for (unsigned s = 0; s < nsamples; s++) {
		tmp = sampler1.get_random() - sampler2.get_random();
		tmp[3] = 0.0f;
		e += FRETefficiency(tmp.norm(), R0);
	}
	return e / nsamples;
}

double meanDistance(const Grid3D &g1, const Grid3D &g2, const unsigned nsamples)
{
	return meanDistanceInv(g1, g2, nsamples);
}

double meanEfficiency(const Grid3D &g1, const Grid3D &g2, const float R0,
		      const unsigned nsamples)
{
	return meanEfficiencyInv(g1, g2, R0, nsamples);
}

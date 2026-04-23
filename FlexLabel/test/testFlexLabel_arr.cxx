#include <FlexLabel/FlexLabel.h>

#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

float checksum(const Grid3D &grid)
{
    return grid.discStep * grid.grid.size()
           + std::accumulate(grid.originXYZ.begin(), grid.originXYZ.end(), 0.0f)
           + std::accumulate(grid.shape.begin(), grid.shape.end(), 0)
           + std::accumulate(grid.grid.begin(), grid.grid.end(), 0.0f);
}

int main()
{
    Eigen::Vector3f source;
    source << 0.0f, -4.0f, 0.0f;

    const int nAtoms = 128;
    Eigen::Matrix4Xf atoms(4, nAtoms);
    atoms.col(0) << 0.0f, -4.0f, 22.0f, 1.5f;
    atoms.col(1) << 9.0f, 0.0f, 0.0f, 3.0f;
    atoms.col(2) << 9.0f, 8.0f, 0.0f, 3.0f;
    atoms.col(3) << 0.0f, -10.5f, 0.0f, 1.5f;
    atoms.col(4) << 5.0f, 0.0f, 0.0f, 2.0f;
    atoms.col(5) << 0.0f, -4.0f, -10.5f, 1.3f;
    for (int i = 6; i < nAtoms; ++i) {
        atoms.col(i) = atoms.col(0);
    }

    const float discStep = 0.9f;
    const float linkerL = 20.0f;
    const float linkerD = 2.0f;
    const float dyeR = 3.5f;
    const Eigen::Vector3f dyeRadii(3.5f, 5.0f, 1.5f);

    const Grid3D ml = minLinkerLength(atoms, source, linkerL, linkerD, dyeR, discStep);
    const Grid3D mlArr = minLinkerLength_arr(atoms, source, linkerL, linkerD, dyeR, discStep);
    assert(std::fabs(checksum(ml) - checksum(mlArr)) < 1e-5f);

    const Grid3D av1 = dyeDensity(atoms, source, linkerL, linkerD, dyeR, discStep);
    const Grid3D av1Arr = dyeDensity_arr(atoms, source, linkerL, linkerD, dyeR, discStep);
    assert(std::fabs(checksum(av1) - checksum(av1Arr)) < 1e-5f);

    const Grid3D av3 = dyeDensity(atoms, source, linkerL, linkerD, dyeRadii, discStep);
    const Grid3D av3Arr = dyeDensity_arr(atoms, source, linkerL, linkerD, dyeRadii, discStep);
    assert(std::fabs(checksum(av3) - checksum(av3Arr)) < 1e-5f);

    Eigen::Matrix<float, 5, Eigen::Dynamic> xyzRQ(5, 2);
    xyzRQ.col(0) << 0.0f, -4.0f, 0.0f, 2.0f, 2.5f;
    xyzRQ.col(1) << 6.0f, -2.0f, 4.0f, 2.0f, 1.5f;
    const Grid3D w = addWeights(av1, xyzRQ);
    const Grid3D wArr = addWeights_arr(av1, xyzRQ);
    assert(std::fabs(checksum(w) - checksum(wArr)) < 1e-5f);

    const double md = meanDistance(av1, av1, 5000);
    const double mdArr = meanDistance_arr(av1, av1, 5000);
    assert(std::fabs(md - mdArr) < 1e-6);

    const double me = meanEfficiency(av1, av1, 52.0f, 5000);
    const double meArr = meanEfficiency_arr(av1, av1, 52.0f, 5000);
    assert(std::fabs(me - meArr) < 1e-6);

    const std::vector<float> d = sampleDistanceDistInv(av1, av1, 2000);
    const std::vector<float> dArr = sampleDistanceDistInv_arr(av1, av1, 2000);
    assert(d.size() == dArr.size());

    return 0;
}

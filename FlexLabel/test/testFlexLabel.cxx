///
/// @file    testFlexLabel.cxx
/// @author  Mykola Dimura
///

#include <FlexLabel/FlexLabel.h>

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

void savePqr(const std::string &fName, const Grid3D &grid)
{
	std::fstream out(fName, std::fstream::out);
	out << std::setprecision(3) << std::setfill(' ');

	using std::setw;
	int idx = 0;
	for (int k = 0; k < grid.shape[2]; ++k) {
		for (int j = 0; j < grid.shape[1]; ++j) {
			for (int i = 0; i < grid.shape[0]; ++i) {
				float val = grid.value(i, j, k);
				if (val <= 0.0f) {
					continue;
				}
				val = std::min(val, 99999.9f);
				auto xyz = grid.xyz(i, j, k);
				out << "ATOM" << setw(7) << idx << "   AV  AV"
				    << setw(6) << idx << setw(12) << xyz[0]
				    << setw(8) << xyz[1] << setw(8) << xyz[2]
				    << setw(8) << val << setw(7)
				    << grid.discStep * 0.5f << std::endl;
				++idx;
			}
		}
	}
}

Eigen::Matrix4Xf readPqr(const std::string &fName)
{
	std::fstream in(fName, std::fstream::in);
	std::string buf((std::istreambuf_iterator<char>(in)),
			std::istreambuf_iterator<char>());
	int nAtoms = 0;
	size_t pos = 0;
	while ((pos = buf.find("ATOM ", pos)) != std::string::npos) {
		++nAtoms;
		++pos;
	}
	Eigen::Matrix4Xf m(4, nAtoms);
	std::istringstream stream(buf);
	std::string tmp;
	int iat = 0;
	while (stream >> tmp) {
		if (tmp != "ATOM") {
			continue;
		}
		stream >> tmp; // id
		stream >> tmp; // name
		stream >> tmp; // res
		stream >> tmp; // resi
		float x, y, z, r;
		stream >> x >> y >> z;
		stream >> tmp; // charge
		stream >> r;
		m.col(iat) << x, y, z, r;
		++iat;
	}
	assert(iat == nAtoms);
	return m;
}

float checksum(const Grid3D &grid)
{
	float sum = grid.discStep * grid.grid.size();
	for (int i = 0; i < 3; ++i) {
		sum += grid.originXYZ[i] + grid.shape[i];
	}
	for (const float &v : grid.grid) {
		sum += v;
	}
	return sum;
}

int main()
{
	using std::cout;
	using std::endl;

	Eigen::Vector3f source;
	Eigen::Matrix4Xf atoms;

	const int nAtoms = 2000;
	atoms.resize(4, nAtoms);
	atoms.col(0) << 0.0f, -4.0f, 22.0f, 1.5f;
	atoms.col(1) << 9.0f, 0.0f, 0.0f, 3.0f;
	atoms.col(2) << 9.0f, 8.0f, 0.0f, 3.0f;
	atoms.col(3) << 0.0f, -10.5f, 0.0f, 1.5f;
	atoms.col(4) << 5.0f, 0.0f, 0.0f, 2.0f;
	atoms.col(5) << 0.0f, -4.0f, -10.5f, 1.3f;
	atoms.col(6) << 0.0f, -4.0f, -11.5f, 1.0f;
	atoms.col(7) << 0.0f, -4.0f, -12.5f, 2.5f;
	atoms.col(8) << 0.0f, -4.0f, -13.5f, 1.7f;
	atoms.col(9) << 0.0f, -4.0f, -14.5f, 1.8f;
	atoms.col(10) << 0.0f, -4.0f, 95.0f, 1.5f;
	for (int iat = 11; iat < nAtoms; ++iat) {
		atoms.col(iat) = atoms.col(10);
	}

	source << 0.0f, -4.0f, 0.0f;

	/*atoms= readPqr("short_helix.pqr");
	atoms(3, 56) = 0.0f; //remove the attachment atom itself
	source = atoms.col(56).head<3>();*/

	const float discStep = 0.9f;
	const float linkerL = 20.0f;
	const float linkerD = 2.0f;
	const Eigen::Vector3f dyeRadii(3.5f, 5.0f, 1.5f);
	const float &dyeR = dyeRadii.coeffRef(0);

	Grid3D grid = minLinkerLength(atoms, source, linkerL, linkerD, dyeR,
				      discStep);
	// savePqr("testMinL.pqr", grid);
	if (fabs(checksum(grid) - 667073.62500f) > 0.00001f) {
		cout << "minLinkerLength() produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		// return 1;
	}

	auto start = std::chrono::steady_clock::now();
	grid = dyeDensity(atoms, source, linkerL, linkerD, dyeR, discStep);
	auto diff = std::chrono::steady_clock::now() - start;
	double dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 19 ms on a laptop with Core i5-4200U CPU
	cout << "FlexLabel version " << GIT_VERSION_STRING << std::endl;
	cout << "AV calculation took: " << dtMs << " ms" << std::endl;

	// savePqr("testDensityAV1.pqr", grid);
	if (fabs(checksum(grid) - 50183.00000f) > 0.00001f) {
		cout << "dyeDensity() AV1 produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		// return 2;
	}

	grid = dyeDensity(atoms, source, linkerL, linkerD, dyeRadii, discStep);
	// savePqr("testDensityAV3.pqr", grid);
	if (fabs(checksum(grid) - 49989.75391f) > 0.00001f) {
		cout << "dyeDensity() AV3 produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		return 3;
	}

	return 0;
}

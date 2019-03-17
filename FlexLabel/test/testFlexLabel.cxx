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


#include "pcg_random.hpp"
#include <random>

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

void printDistanceHist(const std::vector<float> &distances)
{
	const float minD =
		*std::min_element(distances.begin(), distances.end());
	const float maxD =
		*std::max_element(distances.begin(), distances.end());

	const unsigned nBins = 20;
	const float step = (maxD - minD) / (nBins - 1);
	std::vector<int> hist(nBins, 0);
	for (float dist : distances) {
		int iBin = (dist - minD) / step;
		hist[iBin]++;
	}
	std::cout << "Distance\tFreq.\n";
	for (int iBin = 0; iBin < nBins; ++iBin) {
		std::cout << step * (0.5f + iBin) + minD << "\t" << hist[iBin]
			  << "\n";
	}
	std::cout << std::endl;
}

int main()
{
	using std::cout;
	using std::endl;

	cout << "FlexLabel version " << GIT_VERSION_STRING << endl;

	Eigen::Vector3f source;
	source << 0.0f, -4.0f, 0.0f;
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
		return 1;
	}

	auto start = std::chrono::steady_clock::now();
	grid = dyeDensity(atoms, source, linkerL, linkerD, dyeR, discStep);
	auto diff = std::chrono::steady_clock::now() - start;
	double dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 19 ms on a laptop with Core i5-4200U CPU
	cout << "AV1 calculation took: " << dtMs << " ms" << endl;
	// savePqr("testDensityAV1.pqr", grid);

	if (fabs(checksum(grid) - 50183.00000f) > 0.00001f) {
		cout << "dyeDensity() AV1 produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		return 2;
	}

	start = std::chrono::steady_clock::now();
	float Rda = meanDistance(grid, grid, 30000);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 2.7 ms on a Core i7-4930K CPU
	cout << "meanDistance took: " << dtMs << " ms " << endl;
	if (fabs(Rda - 20.35f) > 0.2f) {
		cout << "meanDistance() produced an unexpected result\n";
		cout << "<Rda> = " << std::setprecision(2) << std::fixed << Rda
		     << endl;
		return 3;
	}

	start = std::chrono::steady_clock::now();
	grid = dyeDensity(atoms, source, linkerL, linkerD, dyeRadii, discStep);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 16 ms on a Core i7-4930K CPU
	cout << "AV3 calculation took: " << dtMs << " ms" << endl;
	// savePqr("testDensityAV3.pqr", grid);
	if (fabs(checksum(grid) - 49989.75391f) > 0.00001f) {
		cout << "dyeDensity() AV3 produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		return 4;
	}

	// Contact surface weighting
	Eigen::Matrix<float, 5, Eigen::Dynamic> xyzRQ;
	xyzRQ.resize(5, nAtoms);
	xyzRQ.col(0) << 0.0f, -4.0f, 0.0f, 3.0f, 3.0f;
	xyzRQ.col(1) << 13.0f, -4.0f, 13.0f, 4.0f, 0.12f;
	for (int iat = 2; iat < nAtoms; ++iat) {
		xyzRQ.col(iat) = xyzRQ.col(1);
	}
	start = std::chrono::steady_clock::now();
	grid = addWeights(grid, xyzRQ);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 4.5 ms on a Core i7-4930K CPU
	cout << "addWeights took: " << dtMs << " ms " << endl;
	// savePqr("testContactDensityAV3.pqr", grid);
	if (fabs(checksum(grid) - 117905.71875f) > 0.00001f) {
		cout << "addWeights() produced an unexpected result\n";
		cout << "checksum = " << std::setprecision(5) << std::fixed
		     << checksum(grid) << endl;
		return 5;
	}

	start = std::chrono::steady_clock::now();
	Rda = meanDistance(grid, grid, 30000);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 3.8 ms on a Core i7-4930K CPU
	cout << "meanDistance took: " << dtMs << " ms " << endl;
	const float RdaMeanRef = 14.2f;
	if (fabs(Rda - RdaMeanRef) > 0.2f) {
		cout << "meanDistance() produced an unexpected result\n";
		cout << "<Rda> = " << std::setprecision(2) << std::fixed << Rda
		     << endl;
		return 6;
	}

	auto shiftedGrid = grid;
	shiftedGrid.originXYZ[0] += 52.0;
	start = std::chrono::steady_clock::now();
	auto E = meanEfficiency(grid, shiftedGrid, 52.0, 30000);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 4.0 ms on a Core i7-4930K CPU
	cout << "meanEfficiency took: " << dtMs << " ms " << endl;
	if (fabs(E - 0.479f) > 0.003f) {
		cout << "meanEfficiency() produced an unexpected result\n";
		cout << "<E> = " << std::setprecision(3) << std::fixed << E
		     << endl;
		return 7;
	}

	start = std::chrono::steady_clock::now();
	std::vector<float> distances = sampleDistanceDistInv(grid, grid);
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 50 ms on a Core i5-4200U CPU
	cout << "sampleDistanceDistInv() took: " << dtMs << " ms " << endl;
	Rda = std::accumulate(distances.begin(), distances.end(), 0.0f)
	      / distances.size();
	if (fabs(Rda - RdaMeanRef) > 0.03f) {
		cout << "sampleDistanceDistInv() produced an unexpected result\n";
		cout << "<Rda> = " << std::setprecision(2) << std::fixed << Rda
		     << endl;
		printDistanceHist(distances);
		return 8;
	}

	return 0;
}

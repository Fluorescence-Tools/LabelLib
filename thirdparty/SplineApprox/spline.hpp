//
// Created by dimura on 09/02/18.
//

#ifndef SPLINE_H
#define SPLINE_H

#include "polynomial.hpp"
#include <Eigen/Dense>
#include <array>
#include <Eigen/StdVector>
#include <iostream>

template <unsigned Deg> class Spline
{
      private:
        std::vector<Polynomial<Deg>> polynomials;
        float xmin, xmax;
        float invdx;
        Spline(const unsigned NPieces)
	{
        	polynomials.resize(NPieces);
        }

      public:
	//Spline(const Eigen::MatrixX2f &unsortedXY, const unsigned NPieces);
        float of(const float &x) const
        {
                int i = (x - xmin) * invdx;
                assert(i >= 0);
                assert(i < polynomials.size());
                return polynomials[i].of(x);
        }
	static Spline<Deg> fromSorted(const Eigen::Matrix2Xf& sortedXY, const unsigned NPieces);

        template <unsigned Nd>
        friend std::ostream &operator<<(std::ostream &os,
                                        const Spline<Nd> &spline);
};
template <unsigned Deg>
std::ostream &operator<<(std::ostream &os,
                         const Spline<Deg> &spline)
{
	for (int i = 0; i < spline.polynomials.size(); ++i) {
		float dx = 1.0 / spline.invdx;
		float from = spline.xmin + dx * i;
		float to = from + dx;
		os << '[' << from << ".." << to << "] " << spline.polynomials[i]
		   << '\n';
	}
	return os;
}
/*
template <unsigned Deg>
Spline<Deg>::Spline(const Eigen::MatrixX2f &unsortedXY, const unsigned NPieces)
{
	using Eigen::Matrix;

	polynomials.resize(NPieces);

	const int nRows = unsortedXY.rows();
	using pairf = std::pair<float, float>;
	std::vector<pairf> xy;
	xy.reserve(nRows);
	for (int i = 0; i < nRows; ++i) {
		xy.emplace_back(unsortedXY(i, 0), unsortedXY(i, 1));
	}
	std::sort(xy.begin(), xy.end(),
	          [](pairf a, pairf b) { return a.first < b.first; });
	xmin = xy.front().first;
	xmax = nextafterf(xy.back().first, xy.back().first + 1.0f);
	const float dx = (xmax - xmin) / NPieces;
	invdx = 1.0f / dx;
	using Matrix2XfRM = Matrix<float, Eigen::Dynamic, 2, Eigen::RowMajor>;
	Eigen::Map<Matrix2XfRM> mXY(&(xy.front().first), nRows, 2);

	int rStart;
	for (rStart = 0; rStart < nRows;) {
		float xStart = xy[rStart].first;
		int piece = (xStart - xmin) / dx;
		float xEnd = (piece + 1) * dx + xmin;
		if (xEnd >= xmax) {
			break;
		}

		assert(piece < NPieces - 1);

		int rEnd;
		for (rEnd = rStart; xy[rEnd].first < xEnd; ++rEnd) {
		} // do nothing

		Polynomial<Deg> p(mXY.middleRows(rStart, rEnd - rStart));
		polynomials[piece] = std::move(p);
		rStart = rEnd;
	}
	Polynomial<Deg> p(mXY.middleRows(rStart, nRows - rStart));
	polynomials.back() = std::move(p);
}*/

template<unsigned int Deg>
Spline<Deg> Spline<Deg>::fromSorted(const Eigen::Matrix2Xf& sortedXY, const unsigned NPieces)
{
	Spline<Deg> sp(NPieces);

	sp.xmin=sortedXY(0,0);
	sp.xmax=sortedXY(0,sortedXY.cols()-1);
	sp.xmax=nextafterf(sp.xmax,sp.xmax*2.0f);
	const float dx = (sp.xmax - sp.xmin) / NPieces;
	sp.invdx = 1.0f / dx;

	const int nPoints = sortedXY.cols();
	int iStart;
	for (iStart = 0; iStart < nPoints;) {
		float xStart = sortedXY(0,iStart);
		int piece = (xStart - sp.xmin) * sp.invdx;
		float xEnd = (piece + 1) * dx + sp.xmin;
		if (xEnd >= sp.xmax) {
			break;
		}

		assert(piece < NPieces - 1);

		int iEnd = iStart;
		while (sortedXY(0,iEnd) < xEnd) {
			++iEnd;
		}

		Polynomial<Deg> p(sortedXY.middleCols(iStart, iEnd - iStart));
		sp.polynomials[piece] = std::move(p);
		iStart = iEnd;
	}
	Polynomial<Deg> p(sortedXY.middleCols(iStart, nPoints - iStart));
	sp.polynomials.back() = std::move(p);

	return sp;
}

#endif // SPLINE_H

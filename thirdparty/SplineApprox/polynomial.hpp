#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <Eigen/Dense>
#include <iosfwd>

template <unsigned D> class Polynomial
{
private:
	Eigen::Matrix<float, D + 1, 1> c;
	static Eigen::Matrix<float, D + 1, 1> powers(const float &x)
	{
		Eigen::Matrix<float, D + 1, 1> vec;
		vec[0] = 1.0f;
		for (int i = 1; i < D + 1; ++i) {
			vec[i] = vec[i - 1] * x;
		}
		return vec;
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	inline float of(const float &x) const
	{
		return c.cwiseProduct(powers(x)).sum();
	}
	//Polynomial(const Eigen::MatrixX2f &xy);
	Polynomial(const Eigen::Matrix2Xf &xy);
	Polynomial() = default;
	template <unsigned Deg>
	friend std::ostream &operator<<(std::ostream &os,
	                                const Polynomial<Deg> &p);
};
/*
template <unsigned D> Polynomial<D>::Polynomial(const Eigen::MatrixX2f &xy)
{
	const int N = xy.rows();
	Eigen::Matrix<float, Eigen::Dynamic, D + 1> X(N, D + 1);
	for (int i = 0; i < N; ++i) {
		X.row(i) = powers(xy(i, 0));
	}
	c = X.householderQr().solve(xy.col(1));
}*/


template <unsigned D> Polynomial<D>::Polynomial(const Eigen::Matrix2Xf &xy)
{
	const int N = xy.cols();
	Eigen::Matrix<float, Eigen::Dynamic, D + 1> X(N, D + 1);
	for (int i = 0; i < N; ++i) {
		X.row(i) = powers(xy(0, i));
	}
	c = X.householderQr().solve(xy.row(1).transpose());
}


template <unsigned D>
std::ostream &operator<<(std::ostream &os, const Polynomial<D> &p)
{
	os << p.c[0];
	if (D > 0) {
		os << " + " << p.c[1] << " * x";
	}
	for (int i = 2; i < D + 1; ++i) {
		os << " + " << p.c[i] << " * x^" << i;
	}
	return os;
}

#endif // POLYNOMIAL_H

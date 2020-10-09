#ifndef INTERPOLATION
#define INTERPOLATION

#include "Vect_Utility.hpp"
#include "Constants.hpp"
#include <vector>
#include <complex>
#include <algorithm>
#include <numeric> 
#include <cmath>
#include <limits>

void InterpPattern
(
	double theta, 
	double phi, 
	const std::vector<std::vector<double>>& Gxy, 
	const std::vector<std::vector<double>>& Gxz, 
	const std::vector<std::vector<double>>& Gyz,
	double & GBS,
	double & GBP
);

void InterpLine
(
	double theta, 
	double phi, 
	const std::vector<std::vector<double>>& Gxy, 
	const std::vector<std::vector<double>>& Gxz, 
	const std::vector<std::vector<double>>& Gyz,
	std::complex<double>& gt1,
	std::complex<double>& gt2,
	std::complex<double>& gp1,
	std::complex<double>& gp2
);

std::vector<size_t> sort_indices(const std::vector<double> &v);
void reorder_vector(std::vector<std::complex<double>>& v, const std::vector<size_t>& order);
void reorder_vector(std::vector<double>& v, const std::vector<size_t>& order);

int findNearestNeighbourIndex(double value, std::vector< double > &x );
std::complex<double> interp1(std::vector<double> &x, std::vector<std::complex<double>> &y, double &x_new );

#endif

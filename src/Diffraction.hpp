#ifndef DIFFRACTION
#define DIFFRACTION

#include "Vect_Utility.hpp"
#include "Constants.hpp"
#include <vector>
#include <complex>

void wedge_diff_coeff
(
	double r,
	double ph, 
	double php, 
	double bo, 
	double fn, 
	std::complex<double>& ds, 
	std::complex<double>& dh, 
	std::complex<double>& dps, 
	std::complex<double>& dph
);

void dielec_wedge_diff_coeff
(
	double r,
	double ph, 
	double php, 
	double bo, 
	double fn,
	std::complex<double> gammah,
	std::complex<double> gammas,
	std::complex<double>& ds, 
	std::complex<double>& dh, 
	std::complex<double>& dps, 
	std::complex<double>& dph
);

std::complex<double> di (double r, double bet, double bo, double fn);
std::complex<double> dpi (double r, double bet, double bo, double fn);
std::pair<double,double> frnels_int(double xs);

#endif
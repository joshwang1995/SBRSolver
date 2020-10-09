// Header file used to store constants
#ifndef CONSTANTS
#define CONSTANTS

#include <cmath>
#include <complex>
#define _USE_MATH_DEFINES

const double SPEED_OF_LIGHT = 2.99792458e8;
const double U0 = 4*M_PI*1e-7;
const double E0 = 1.0/(U0*SPEED_OF_LIGHT*SPEED_OF_LIGHT);

// Threshold value used to compare a double or float number to zero
const double SMALL_DOUBLE = 1e-5;

// Value used to initialize vectors
const double LARGE_DOUBLE = 1e300;

// Easy way to reference 2PI
const double TWOPI = 2*M_PI;
const double PI = M_PI;

// Complex number
const std::complex<double> j (0.0,1.0);

#endif
#pragma once

#include <complex>
#include <limits>

const double SPEED_OF_LIGHT = 2.99792458e8; // Speed of light in vaccum
const double U0 = 4*EIGEN_PI*1e-7; // permeability of free space 
const double E0 = 1.0/(U0*SPEED_OF_LIGHT*SPEED_OF_LIGHT); // Permittivity of free space
const double ETA0 = U0 * SPEED_OF_LIGHT; // Free space impedance

// Threshold value used to compare a double or float number to zero
const double SMALL_DOUBLE = 1e-5;
const double EPSILON = std::numeric_limits<float>::epsilon();

// Value used to initialize vectors
const double INF = std::numeric_limits<double>::infinity();

// Easy way to reference 2PI
const double TWOPI = 2* EIGEN_PI;
const double PI = EIGEN_PI;

// Complex number
const std::complex<double> j (0.0,1.0);
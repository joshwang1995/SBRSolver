#pragma once
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

struct Vect3f
{
	float x, y, z;
};

struct Vect3d
{
	double x, y, z;
};

struct Vect3u
{
	unsigned v0, v1, v2;
};

#pragma region Constructors
Vect3f make_Vect3f(float x, float y, float z);
Vect3f make_Vect3f(double x, double y, double z);
Vect3f make_Vect3f(const float s);
Vect3d make_Vect3d(double x, double y, double z);
Vect3d make_Vect3d(const double s);
#pragma endregion

#pragma region Operators
Vect3f operator+(const Vect3f& a, const Vect3f& b);
Vect3d operator+(const Vect3d& a, const Vect3d& b);
Vect3f operator+(const Vect3f& a, const float b);
Vect3d operator+(const Vect3d& a, const double b);
Vect3f operator+(const float a, const Vect3f& b);
Vect3d operator+(const double a, const Vect3d& b);
void operator+=(Vect3f& a, const Vect3f& b);
void operator+=(Vect3d& a, const Vect3d& b);
Vect3f operator-(const Vect3f& a, const Vect3f& b);
Vect3d operator-(const Vect3d& a, const Vect3d& b);
Vect3f operator-(const Vect3f& a, const float b);
Vect3d operator-(const Vect3d& a, const double b);
Vect3f operator-(const float a, const Vect3f& b);
Vect3d operator-(const double a, const Vect3d& b);
void operator-=(Vect3f& a, const Vect3f& b);
void operator-=(Vect3d& a, const Vect3d& b);
Vect3f operator-(const Vect3f& a);
Vect3d operator-(const Vect3d& a);
Vect3f operator*(const Vect3f& a, const Vect3f& b);
Vect3d operator*(const Vect3d& a, const Vect3d& b);
Vect3f operator*(const Vect3f& a, const float s);
Vect3f operator*(const Vect3f& a, const double s);
Vect3d operator*(const Vect3d& a, const double s);
Vect3f operator*(const float s, const Vect3f& a);
Vect3d operator*(const double s, const Vect3d& a);
void operator*=(Vect3f& a, const Vect3f& s);
void operator*=(Vect3d& a, const Vect3d& s);
void operator*=(Vect3f& a, const float s);
void operator*=(Vect3d& a, const double s);
Vect3f operator/(const Vect3f& a, const Vect3f& b);
Vect3d operator/(const Vect3d& a, const Vect3d& b);
Vect3f operator/(const Vect3f& a, const float s);
Vect3d operator/(const Vect3d& a, const double s);
Vect3f operator/(const float s, const Vect3f& a);
Vect3d operator/(const double s, const Vect3d& a);
void operator/=(Vect3f& a, const float s);
void operator/=(Vect3d& a, const double s);
#pragma endregion

#pragma region Vector Operations
float dot(const Vect3f& a, const Vect3f& b);
double dot(const Vect3d& a, const Vect3d& b);
float length(const Vect3f& v);
double length(const Vect3d& v);
Vect3f cross(const Vect3f& a, const Vect3f& b);
Vect3d cross(const Vect3d& a, const Vect3d& b);
Vect3f normalize(const Vect3f& v);
Vect3d normalize(const Vect3d& v);
Vect3f reflect(const Vect3f& i, const Vect3f& n);
Vect3d reflect(const Vect3d& i, const Vect3d& n);
float DistanceBetween(Vect3f pt1, Vect3f pt2);
float AngleBetween(Vect3f pt1, Vect3f pt2, bool alreadyNormalized = false, bool resulInRadian = false);
double AngleBetween(Vect3d pt1, Vect3d pt2, bool alreadyNormalized = false, bool resulInRadian = false);
#pragma endregion

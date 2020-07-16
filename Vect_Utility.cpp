#include "Vect_Utility.hpp"

#pragma region Constructors
Vect3f make_Vect3f(float x, float y, float z)
{
	Vect3f t_Vect{ x, y, z };
	return t_Vect;
}

Vect3f make_Vect3f(double x, double y, double z)
{
    Vect3f t_Vect{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) };
    return t_Vect;
}

Vect3f make_Vect3f(const float s)
{
	return make_Vect3f(s, s, s);
}

Vect3d make_Vect3d(double x, double y, double z)
{
	Vect3d t_Vect{ x, y, z };
	return t_Vect;
}

Vect3d make_Vect3d(const double s)
{
	return make_Vect3d(s, s, s);
}
#pragma endregion

#pragma region Operators
Vect3f operator+(const Vect3f& a, const Vect3f& b)
{
	return make_Vect3f(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vect3d operator+(const Vect3d& a, const Vect3d& b)
{
	return make_Vect3d(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vect3f operator+(const Vect3f& a, const float b)
{
	return make_Vect3f(a.x + b, a.y + b, a.z + b);
}

Vect3d operator+(const Vect3d& a, const double b)
{
	return make_Vect3d(a.x + b, a.y + b, a.z + b);
}

Vect3f operator+(const float a, const Vect3f& b)
{
	return make_Vect3f(a + b.x, a + b.y, a + b.z);
}

Vect3d operator+(const double a, const Vect3d& b)
{
	return make_Vect3d(a + b.x, a + b.y, a + b.z);
}

void operator+=(Vect3f& a, const Vect3f& b)
{
	a.x += b.x; a.y += b.y; a.z += b.z;
}

void operator+=(Vect3d& a, const Vect3d& b)
{
	a.x += b.x; a.y += b.y; a.z += b.z;
}

Vect3f operator-(const Vect3f& a, const Vect3f& b)
{
	return make_Vect3f(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vect3d operator-(const Vect3d& a, const Vect3d& b)
{
	return make_Vect3d(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vect3f operator-(const Vect3f& a, const float b)
{
	return make_Vect3f(a.x - b, a.y - b, a.z - b);
}

Vect3d operator-(const Vect3d& a, const double b)
{
	return make_Vect3d(a.x - b, a.y - b, a.z - b);
}

Vect3f operator-(const float a, const Vect3f& b)
{
	return make_Vect3f(a - b.x, a - b.y, a - b.z);
}

Vect3d operator-(const double a, const Vect3d& b)
{
	return make_Vect3d(a - b.x, a - b.y, a - b.z);
}

void operator-=(Vect3f& a, const Vect3f& b)
{
	a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

void operator-=(Vect3d& a, const Vect3d& b)
{
	a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

Vect3f operator-(const Vect3f& a)
{
	return make_Vect3f(-a.x, -a.y, -a.z);
}

Vect3d operator-(const Vect3d& a)
{
	return make_Vect3d(-a.x, -a.y, -a.z);
}

Vect3f operator*(const Vect3f& a, const Vect3f& b)
{
	return make_Vect3f(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vect3d operator*(const Vect3d& a, const Vect3d& b)
{
	return make_Vect3d(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vect3f operator*(const Vect3f& a, const float s)
{
	return make_Vect3f(a.x * s, a.y * s, a.z * s);
}

Vect3f operator*(const Vect3f& a, const double s)
{    
    return make_Vect3f(static_cast<float>(a.x * s), static_cast<float>(a.y * s), static_cast<float>(a.z * s));
}

Vect3d operator*(const Vect3d& a, const double s)
{
	return make_Vect3d(a.x * s, a.y * s, a.z * s);
}

Vect3f operator*(const float s, const Vect3f& a)
{
	return make_Vect3f(a.x * s, a.y * s, a.z * s);
}

Vect3d operator*(const double s, const Vect3d& a)
{
	return make_Vect3d(a.x * s, a.y * s, a.z * s);
}

void operator*=(Vect3f& a, const Vect3f& s)
{
	a.x *= s.x; a.y *= s.y; a.z *= s.z;
}

void operator*=(Vect3d& a, const Vect3d& s)
{
	a.x *= s.x; a.y *= s.y; a.z *= s.z;
}

void operator*=(Vect3f& a, const float s)
{
	a.x *= s; a.y *= s; a.z *= s;
}

void operator*=(Vect3d& a, const double s)
{
	a.x *= s; a.y *= s; a.z *= s;
}

Vect3f operator/(const Vect3f& a, const Vect3f& b)
{
	return make_Vect3f(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vect3d operator/(const Vect3d& a, const Vect3d& b)
{
	return make_Vect3d(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vect3f operator/(const Vect3f& a, const float s)
{
	float inv = 1.0f / s;
	return a * inv;
}

Vect3d operator/(const Vect3d& a, const double s)
{
	double inv = 1.0 / s;
	return a * inv;
}

Vect3f operator/(const float s, const Vect3f& a)
{
	return make_Vect3f(s / a.x, s / a.y, s / a.z);
}

Vect3d operator/(const double s, const Vect3d& a)
{
	return make_Vect3d(s / a.x, s / a.y, s / a.z);
}

void operator/=(Vect3f& a, const float s)
{
	float inv = 1.0f / s;
	a *= inv;
}

void operator/=(Vect3d& a, const double s)
{
	double inv = 1.0 / s;
	a *= inv;
}
#pragma endregion

#pragma region Vector Operations
float dot(const Vect3f& a, const Vect3f& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double dot(const Vect3d& a, const Vect3d& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

float length(const Vect3f& v)
{
	return sqrt(dot(v, v));
}

double length(const Vect3d& v)
{
	return sqrt(dot(v, v));
}

Vect3f cross(const Vect3f& a, const Vect3f& b)
{
	return make_Vect3f(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

Vect3d cross(const Vect3d& a, const Vect3d& b)
{
	return make_Vect3d(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

Vect3f normalize(const Vect3f& v)
{
	float invLen = 1.0f / sqrt(dot(v, v));
	return v * invLen;
}

Vect3d normalize(const Vect3d& v)
{
	double invLen = 1.0 / sqrt(dot(v, v));
	return v * invLen;
}

//Reflect ray is calculated as --> inRay - 2 * wallNormal * Dot(wallNormal, inRay)
Vect3f reflect(const Vect3f& i, const Vect3f& n)
{
	return i - 2.0f * n * dot(n, i);
}

Vect3d reflect(const Vect3d& i, const Vect3d& n)
{
	return i - 2.0 * n * dot(n, i);
}

float DistanceBetween(Vect3f pt1, Vect3f pt2)
{
	return length(pt2 - pt1);
}

float AngleBetween(Vect3f pt1, Vect3f pt2, bool alreadyNormalized, bool resulInRadian)
{
	if (!alreadyNormalized)
	{
		pt1 = normalize(pt1);
		pt2 = normalize(pt2);
	}

	float ratio = dot(pt1, pt2);
	float theta;

	if (ratio < 0)
	{
		theta = static_cast<float>(M_PI - 2.0f * std::asin(length(-pt1 - pt2) / 2.0f));
	}
	else
	{
		theta = 2.0f * std::asin(length(pt1 - pt2) / 2.0f);
	}

	return resulInRadian ? theta : static_cast<float>(theta*(180.0 / M_PI));
}

double AngleBetween(Vect3d pt1, Vect3d pt2, bool alreadyNormalized, bool resulInRadian)
{
	if (!alreadyNormalized)
	{
		pt1 = normalize(pt1);
		pt2 = normalize(pt2);
	}

	double ratio = dot(pt1, pt2);
	double theta;

	if (ratio < 0)
	{
		theta = M_PI - 2.0 * std::asin(length(-pt1 - pt2) / 2.0);
	}
	else
	{
		theta = 2.0 * std::asin(length(pt1 - pt2) / 2.0);
	}

	return resulInRadian ? theta : theta * (180.0 / M_PI);
}
#pragma endregion
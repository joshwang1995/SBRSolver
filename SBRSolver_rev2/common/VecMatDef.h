//
// (c) 2017-2018 Takahiro Hashimoto
//
#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>

// #define Dynamic Eigen::Dynamic

typedef	Eigen::Vector2d		Vec2;	// (x, y) value
typedef Eigen::Vector2cd    Vec2c;
typedef	Eigen::Vector2i		Idx2;	// (x, y) location/idx
typedef	Eigen::Vector3d		Vec3;	// (x, y, z) value
typedef Eigen::Vector3cd    Vec3c;
typedef	Eigen::Vector3i		Idx3;	// (x, y, z) location/idx
typedef Eigen::Matrix<bool, 3, 1> Bool3;
typedef	Eigen::Vector4d		Vec4;
typedef	Eigen::Vector4i		Idx4;
typedef	Eigen::VectorXd		Vec;
typedef Eigen::Matrix<double, 6, 1> Vec6; //for (angle-axis, translation)
typedef Eigen::Matrix<double, 2, 3> Mat23;

typedef	Eigen::Matrix2d		Mat2;
typedef	Eigen::Matrix3d		Mat3;
typedef	Eigen::Matrix4d		Mat4;
typedef	Eigen::MatrixXd		Mat;
typedef	Eigen::MatrixXi		Mati;

typedef Eigen::Matrix2cd	Mat2c;
typedef Eigen::Matrix<Vec2c, Eigen::Dynamic, Eigen::Dynamic> MatXc2;

typedef	std::vector<int>									VecIdx;
typedef	std::vector<VecIdx>									VecVecIdx;
typedef std::vector<bool>									VecBool;

typedef	std::vector<Vec2, Eigen::aligned_allocator<Vec2>>	VecVec2;
typedef	std::vector<Idx2, Eigen::aligned_allocator<Idx2>>	VecIdx2;
typedef	std::vector<Vec3, Eigen::aligned_allocator<Vec3>>	VecVec3;
typedef	std::vector<Idx3, Eigen::aligned_allocator<Idx3>>	VecIdx3;
typedef	std::vector<Vec4, Eigen::aligned_allocator<Vec4>>	VecVec4;
typedef	std::vector<Idx4, Eigen::aligned_allocator<Idx4>>	VecIdx4;

typedef std::vector<Mat2, Eigen::aligned_allocator<Mat2>>   VecMat2;
typedef std::vector<Mat2c, Eigen::aligned_allocator<Mat2c>> VecMat2c;

typedef std::vector<std::complex<double>>    Vecc;
typedef std::vector<double>    Vecd;
typedef Eigen::dcomplex cdouble;

inline Vec3 Reflect(const Vec3& i, const Vec3& n)
{
	return i - 2.0 * n * n.dot(i);
}

inline double AngleBetween(Vec3 v1, Vec3 v2, bool alreadyNormalized = false, bool resulInRadian = false)
{
	if (!alreadyNormalized)
	{
		v1.normalize();
		v2.normalize();
	}

	double ratio = v1.dot(v2);
	double theta;

	if (ratio < 0)
	{
		theta = static_cast<double>(EIGEN_PI - 2.0 * std::asin((-v1 - v2).norm() / 2.0));
	}
	else
	{
		theta = static_cast<double>(2.0 * std::asin((v1 - v2).norm() / 2.0));
	}

	return resulInRadian ? theta : (static_cast<double>(theta * (180.0 / EIGEN_PI)));
}

inline Mat3 RotationMatrix(const Mat3& oldSys, const Mat3& newSys)
{
	Mat3 result = newSys * oldSys.transpose();
	return result;
}

inline Vec3 PointInNewCoordSys(const Vec3& point, const Mat3& oldSys, const Mat3& newSys)
{
	return RotationMatrix(oldSys, newSys) * point;
}

inline Vec3 CartesianToSpherical(const Vec3& v)
{
	Vec3 result;
	result(0) = v.norm();
	result(1) = atan(sqrt(v.x() * v.x() + v.y() * v.y()) / v.z()); // Note: atan restricts the result to be between -pi/2 to pi/2
	result(2) = atan(v.y() / v.x() );
	return result;
}

inline Vec3 SphericalToCartesian(const Vec3& v)
{
	Vec3 result;
	result(0) = v(0) * sin(v(1)) * cos(v(2));
	result(1) = v(0) * sin(v(1)) * sin(v(2));
	result(2) = v(0) * sin(v(1));
	return result;
}
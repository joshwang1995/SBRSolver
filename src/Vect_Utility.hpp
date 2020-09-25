#ifndef VECT_UTIL
#define VECT_UTIL

#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <complex>
#include "Constants.hpp"

template <typename T>
class Vect3
{
    public:
        T x,y,z;
        
        //Constructors
        Vect3<T>(T x_ = 0.0, T y_=0.0, T z_=0.0)
            :x(x_), y(y_), z(z_) { };
        Vect3<T>(const Vect3<T>& v)
            :x(v.x), y(v.y), z(v.z) { };
        
        //Operators
        Vect3<T> operator -() const;    
        Vect3<T>& operator +=(const Vect3<T>&);
        Vect3<T>& operator +=(const T&);
        Vect3<T>& operator -=(const Vect3<T>&);
        Vect3<T>& operator -=(const T&);
        Vect3<T>& operator *=(const Vect3<T>&);
        Vect3<T>& operator *=(const T&);
        Vect3<T>& operator /=(const Vect3<T>&);
        Vect3<T>& operator /=(const T&);
		
        Vect3<T> operator +(const Vect3<T>&) const;
        Vect3<T> operator +(const T&) const;
        Vect3<T> operator -(const Vect3<T>&) const;
        Vect3<T> operator -(const T&) const;
        Vect3<T> operator *(const Vect3<T>&) const;
        Vect3<T> operator *(const T&) const;
        Vect3<T> operator /(const Vect3<T>&) const;
        Vect3<T> operator /(const T&) const;
		
		bool operator==(const Vect3<T>&);
		bool operator!=(const Vect3<T>&);
		
		// Overloading for LHS of the operator of type U
		//	Example: (double a) / (Vect3<double> b)
		
		template<typename U>
		friend Vect3<U> operator + (const U&, const Vect3<U>&);
		
		template<typename U>
		friend Vect3<U> operator - (const U&, const Vect3<U>&);
		
		template<typename U>
		friend Vect3<U> operator * (const U&, const Vect3<U>&);
		
		template<typename U>
		friend Vect3<U> operator / (const U&, const Vect3<U>&);
        
        template<typename U>
        friend std::ostream& operator <<(std::ostream&, const Vect3<U>&);
};

template<typename T>
inline Vect3<T> Vect3<T>::operator -() const 
{
    return Vect3(-x, -y, -z);
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator +=(const Vect3<T>& v2) 
{
    x += v2.x;
    y += v2.y;
    z += v2.z;
    return *this;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator +=(const T& c) 
{
    x += c;
    y += c;
    z += c;
    return *this;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator -=(const Vect3<T>& v2)
{
    return *this += -v2;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator -=(const T& c)
{
    return *this += -c;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator *=(const Vect3<T>& v2)
{
    x *= v2.x;
    y *= v2.y;
    z *= v2.z;
    return *this;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator *=(const T& c)
{
    x *= c;
    y *= c;
    z *= c;
    return *this;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator /=(const Vect3<T>& v2)
{
    x *= (1.0/v2.x);
    y *= (1.0/v2.y);
    z *= (1.0/v2.z);
    return *this;
}

template<typename T>
inline Vect3<T>& Vect3<T>::operator /=(const T& c)
{
	T inv = 1.0 / c;
    x *= inv;
    y *= inv;
    z *= inv;
    return *this;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator +(const Vect3<T>& v2) const 
{
    return Vect3<T>(*this) += v2;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator +(const T& c) const 
{
    return Vect3<T>(*this) += c;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator -(const Vect3<T>& v2) const 
{
    return Vect3<T>(*this) -= v2;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator -(const T& c) const 
{
    return Vect3<T>(*this) -= c;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator *(const Vect3<T>& v2) const 
{
    return Vect3<T>(*this) *= v2;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator *(const T& c) const 
{
    return Vect3<T>(*this) *= c;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator /(const Vect3<T>& v2) const 
{
    return Vect3<T>(*this) /= v2;
}

template<typename T>
inline Vect3<T> Vect3<T>::operator /(const T& c) const 
{
    return Vect3<T>(*this) /= c;
}

template<typename T>
inline bool Vect3<T>::operator==(const Vect3<T>& b)
{
	T diff_x = fabs(x-b.x);
	T diff_y = fabs(y-b.y);
	T diff_z = fabs(z-b.z);
	T diff_tot = diff_x + diff_y + diff_z;
	return diff_tot < SMALL_DOUBLE? true:false;
}

template<typename T>
inline bool Vect3<T>::operator!=(const Vect3<T>& b)
{
	T diff_x = fabs(x-b.x);
	T diff_y = fabs(y-b.y);
	T diff_z = fabs(z-b.z);
	T diff_tot = diff_x + diff_y + diff_z;
	return diff_tot < SMALL_DOUBLE? false:true;
}

template<typename T>
inline Vect3<T> operator + (const T& c, const Vect3<T>& v)
{
	return Vect3<T>(v+c);
}

template<typename T>
inline Vect3<T> operator - (const T& c, const Vect3<T>& v)
{
	return Vect3<T>(-v+c);
}

template<typename T>
inline Vect3<T> operator * (const T& c, const Vect3<T>& v)
{
	return Vect3<T>(v*c);
}

template<typename T>
inline Vect3<T> operator / (const T& c, const Vect3<T>& v)
{
	return Vect3<T>(c/v.x,c/v.y,c/v.z);
}

template<typename T>
inline std::ostream& operator <<(std::ostream& os, const Vect3<T>& v) 
{
    os << "<" << v.x << "," << v.y << "," << v.z << ">";
    return os;
}

template<typename T, typename U>
inline U dot(const Vect3<T>& v1, const Vect3<U>& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename T>
inline T length(const Vect3<T>& v)
{
	return sqrt(dot(v,v));
}

template<typename T>
inline Vect3<T> cross(const Vect3<T>& v1, const Vect3<T>& v2)
{
	return Vect3<T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template<typename T>
inline Vect3<T> normalize(const Vect3<T>& v)
{
	T invLen = 1.0/length(v);
	return v * invLen;
}

template<typename T>
inline Vect3<T> reflect(const Vect3<T>& i,const Vect3<T>& n)
{
	return i - 2.0 * n * dot(n,i);
}

template<typename T>
inline T DistanceBetween(const Vect3<T>& v1, const Vect3<T>& v2)
{
	return length(v2-v1);
}

template<typename T>
T AngleBetween(Vect3<T> v1, Vect3<T> v2, bool alreadyNormalized = false, bool resulInRadian = false)
{	
	if (!alreadyNormalized)
	{
		v1 = normalize(v1);
		v2 = normalize(v2);
	}

	T ratio = dot(v1,v2);
	T theta;
	
	if (ratio < 0)
	{
		theta = static_cast<T>(M_PI - 2.0 * std::asin(length(-v1 - v2) / 2.0));
	}
	else
	{
		theta = static_cast<T>(2.0 * std::asin(length(v1-v2)/ 2.0));
	}

	return resulInRadian ? theta : (static_cast<T>(theta*(180.0 / M_PI)));
}

template <>
class Vect3 <unsigned>
{
    public:
        unsigned v0,v1,v2;
		
        friend std::ostream& operator <<(std::ostream& os, const Vect3<unsigned>& v)
		{
			os << "(" << v.v0 << "," << v.v1 << "," << v.v2 << ")";
			return os;
		}
};

template <typename T>
class Vect3Sph: public Vect3<T>
{
    public:
        T r,theta,phi;
		
		//Constructors
        Vect3Sph(T r_ = 0.0, T theta_=0.0, T phi_=0.0) : Vect3<T>(r_,theta_,phi_)
		{
			r = r_;
			theta = theta_;
			phi = phi_;
		};
		
        Vect3Sph(const Vect3Sph& v) : Vect3<T>(v)
		{
			r = v.x;
			theta = v.y;
			phi = v.z;
		};
		 
		
        friend std::ostream& operator <<(std::ostream& os, const Vect3Sph<T>& v)
		{
			os << "(" << v.r << "," << v.theta << "," << v.phi << ")";
			return os;
		}
};

template <typename T>
struct Ray_T
{
	Vect3<T> orig;
	Vect3<T> dir;
};

template <typename T>
struct RectCoord_T
{
	Vect3<T> ax;
	Vect3<T> ay;
	Vect3<T> az;
};

template <typename T>
struct SphCoord_T
{
	Vect3Sph<T> ar;
	Vect3Sph<T> at;
	Vect3Sph<T> ap;
};

template <typename T, typename U>
Vect3<U> VectRectToRect(const RectCoord_T<T>& old_sys, const RectCoord_T<T>& new_sys, const Vect3<U>& p)
{
	RectCoord_T<T> rotation_matrix;

    rotation_matrix.ax = Vect3<T>(dot(old_sys.ax,new_sys.ax), dot(old_sys.ay, new_sys.ax), dot(old_sys.az,new_sys.ax));
    rotation_matrix.ay = Vect3<T>(dot(old_sys.ax,new_sys.ay), dot(old_sys.ay, new_sys.ay), dot(old_sys.az,new_sys.ay));
    rotation_matrix.az = Vect3<T>(dot(old_sys.ax,new_sys.az), dot(old_sys.ay, new_sys.az), dot(old_sys.az,new_sys.az));
	
	Vect3<U> result;
	
    result.x = dot(rotation_matrix.ax,p);
    result.y = dot(rotation_matrix.ay,p);
    result.z = dot(rotation_matrix.az,p);
	return result;
}

template <typename T, typename U>
Vect3<T> VectSphToRect(const Vect3Sph<T>& vect_sph, U theta, U phi)
{
	U st = sin(theta); 
	U sp = sin(phi) ; 
	U ct = cos(theta); 
	U cp = cos(phi) ; 
	
	RectCoord_T<U> Rsc {Vect3<U>(st*cp,ct*cp,-1*sp),Vect3<U>(st*sp,ct*sp,cp),Vect3<U>(ct,-1*st,0)};

	Vect3<T> result;
	result.x = Rsc.ax.x*vect_sph.r + Rsc.ax.y*vect_sph.theta + Rsc.ax.z*vect_sph.phi;
	result.y = Rsc.ay.x*vect_sph.r + Rsc.ay.y*vect_sph.theta + Rsc.ay.z*vect_sph.phi;
	result.z = Rsc.az.x*vect_sph.r + Rsc.az.y*vect_sph.theta + Rsc.az.z*vect_sph.phi;
	
	return result;
}

template <typename T, typename U>
Vect3Sph<T> VectRectToSph(const Vect3<T>& vect, U theta, U phi)
{
	U st = sin(theta); 
	U sp = sin(phi) ; 
	U ct = cos(theta); 
	U cp = cos(phi) ; 
	
	RectCoord_T<U> Rsc {Vect3<U>(st*cp,st*sp,ct),Vect3<U>(ct*cp,ct*sp,-st),Vect3<U>(-sp,cp,0)};

	Vect3Sph<T> result;
	result.r = Rsc.ax.x*vect.x + Rsc.ax.y*vect.y + Rsc.ax.z*vect.z;
	result.theta = Rsc.ay.x*vect.x + Rsc.ay.y*vect.y + Rsc.ay.z*vect.z;
	result.phi = Rsc.az.x*vect.x + Rsc.az.y*vect.y + Rsc.az.z*vect.z;
	
	return result;
}


template <typename T>
Vect3Sph<T> PointRectToSph(const Vect3<T>& a)
{
    Vect3Sph<T> result;
	
    // In this case, x,y,z coresponds to r,theta,phi
	result.r = sqrt(a.x*a.x + a.y*a.y +a.z*a.z);
	result.theta = atan(sqrt(a.x*a.x + a.y*a.y)/a.z); // Note: atan restricts the result to be between -pi/2 to pi/2
	result.phi = atan(a.y/a.x);
	return result;
}

template <typename T>
Vect3<T> PointSphToRect(const Vect3Sph<T>& a)
{
    Vect3<T> result;
	result.x = a.r*sin(a.theta)*cos(a.phi);
	result.y = a.r*sin(a.theta)*sin(a.phi);
	result.z = a.r*sin(a.theta);
	return result;
}

typedef Vect3<unsigned> Vect3u;
typedef Vect3<float> Vect3f;

typedef Vect3<double> Vect3d;
typedef Vect3<std::complex<double>> Cvect3d;
typedef Vect3Sph<double> Vect3dSph;
typedef Vect3Sph<std::complex<double>> Cvect3dsph;

typedef Ray_T<double> Ray;
typedef RectCoord_T<double> RectCoord;

#endif
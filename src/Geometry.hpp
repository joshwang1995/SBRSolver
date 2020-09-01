#pragma once
#include "Vect_Utility.hpp"
#include "Constants.hpp"
#include <vector>

struct Ray
{
	Vect3d orig;
	Vect3d dir;
};

struct RectCoord
{
	Vect3d ax;
	Vect3d ay;
	Vect3d az;
};

class Geometry
{
	public:
		//Virtual functions
		virtual bool Intersects(const Ray& inc_ray, Ray& ref_ray) = 0;
		
		//Static functions
		static Vect3d CalcUnitNorm(const std::vector<Vect3d>& p);
		static bool IsCoplanar(const std::vector<Vect3d>& p);
		static bool IsValidPolygon(const std::vector<Vect3d>& p);
};


class FinitePlane: public Geometry
{
	protected:
		std::vector<Vect3d> vertices_;
		Vect3d unit_norm_;
		double d_;
		
	public:
		FinitePlane(std::vector<Vect3d> vertices);
		FinitePlane(std::vector<Vect3d> vertices, Vect3d unit_norm, double d);
		~FinitePlane();
		
		bool Intersects(const Ray& inc_ray, Ray& ref_ray);
		Ray Compute_Reflect_Ray(const Ray& inc_ray);
};
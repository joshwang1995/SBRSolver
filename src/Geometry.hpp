#pragma once
#include "Vect_Utility.hpp"
#include <vector>

struct Ray
{
	Vect3d orig;
	Vect3d dir;
};

class Geometry
{
	protected:
		std::vector<Vect3d> vertices;
		Vect3d unit_norm;
		
	public:
		//Constructor and destructor
		Geometry();
		virtual ~Geometry();
		
		//Member functions
		virtual bool Intersects(Ray ray, Ray& ref_ray) = 0;
		
		//Static functions are class object independent
		static Vect3d CalcUnitNorm(std::vector<Vect3d> p);
		static bool IsCoplanar(std::vector<Vect3d> p);
		static bool IsValidPolygon(std::vector<Vect3d> p);
};


class FinitePlane: public Geometry
{
	protected:
		double d;
		
	public:
		FinitePlane(std::vector<Vect3d> vertices_);
		FinitePlane(Vect3d unit_norm_, double d);
		~FinitePlane();
		
		bool Intersects(Ray ray, Ray& ref_ray);
};
#pragma once
#include "Vect_Utility.hpp"
#include "Constants.hpp"
#include <vector>

class Geometry
{
	public:
		//Virtual functions
		virtual bool Intersects(const Ray&, Vect3d&) = 0;
		
		//Static functions
		static Vect3d CalcUnitNorm(const std::vector<Vect3d>&);
		static bool IsCoplanar(const std::vector<Vect3d>&);
		static bool IsValidPolygon(const std::vector<Vect3d>&);
        static bool IsPointInPolygon(const std::vector<Vect3d>& v, const Vect3d&);
};


class FinitePlane: public Geometry
{
	public:
		std::vector<Vect3d> vertices;
		Vect3d unit_norm;
		double d;
		FinitePlane(std::vector<Vect3d>);
		FinitePlane(std::vector<Vect3d>, Vect3d, double);
		
		bool Intersects(const Ray&, Vect3d&);
};

/*
class Box: public Geometry
{
	public:
		std::vector<Vect3d> vertices;
		std::vector<Vect3u> faces;
		
		Box(double length, double width, double height, Vect3d center); 
}
*/
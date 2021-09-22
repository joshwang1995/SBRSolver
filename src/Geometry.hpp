#ifndef GEOMETRY
#define GEOMETRY

#include "Vect_Utility.hpp"
#include "Constants.hpp"
#include <vector>
#include <memory>

struct Interaction
{
    Point3d p_int;
	Vect3d normal;
	double ray_length;
};

class Geometry
{
	public:
		//Virtual functions
		virtual bool Intersects(const Ray&, Interaction&) const = 0;
		
		//Static functions
		static Vect3d CalcUnitNorm(const std::vector<Vect3d>&);
		static bool IsCoplanar(const std::vector<Vect3d>&);
		static bool IsValidPolygon(const std::vector<Vect3d>&);
        static bool IsPointInPolygon(const std::vector<Vect3d>& v, const Vect3d&);
};

class Geometry_List: Geometry
{
    public:
        // Instance variable
        std::vector<std::shared_ptr<Geometry>> objects;
        
        Geometry_List() {}
        Geometry_List(std::shared_ptr<Geometry> object) { add(object); }

        void clear() { objects.clear(); }
        void add(std::shared_ptr<Geometry> object) { objects.push_back(object); }
        virtual bool Intersects(const Ray&, Interaction&) const override; 
};

class RectWall: public Geometry
{
	public:
        //Constructors
		RectWall(std::vector<Vect3d>, double, double, double);
		RectWall(std::vector<Vect3d>, Vect3d, double, double, double, double);
        
        //Instance functions
		virtual bool Intersects(const Ray&, Interaction&) const override;
	
	public:
		// Wall geometry properties
		std::vector<Vect3d> vertices;
		Vect3d unit_norm;
		double d; // d for plane (ax + by + cz + d = 0), radius of the sphere
		double depth;
		
		// Wall electrical properties
		double sigma;
		double rel_perm;
};

class ReceptionSphere: public Geometry
{
	public:
		//Constructors
		ReceptionSphere() {};
		ReceptionSphere(Point3d cen, double r): center(cen), radius(r) {};

		//Instance functions
		virtual bool Intersects(const Ray&, Interaction&) const override;
		inline Vect3d get_sphere_norm(Point3d);

	public:
		Point3d center;
		double radius;
};

#endif
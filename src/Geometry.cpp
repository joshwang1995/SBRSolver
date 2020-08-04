#include "Geometry.hpp"
//#include <iostream>

// *******************************************
// Geometry Class

// Constructor and destructor
Geometry::Geometry()
{
	unit_norm = make_Vect3d(0,0,0);
}

Geometry::~Geometry()
{
	vertices.clear();
	unit_norm = make_Vect3d(0,0,0);
}

// Static member function to calculate unit normal of a geometry
Vect3d Geometry::CalcUnitNorm(std::vector<Vect3d> p)
{
	Vect3d n1 = p[1] - p[0];
	Vect3d n2 = p[2] - p[1];
	
	Vect3d n = cross(n1,n2);
	
	n = normalize(n);
	
	return n;
}

bool Geometry::IsCoplanar(std::vector<Vect3d> p)
{
	// (p2 - p0) DOT [(p1-p0) x (p3-p2)] = 0
	// a DOT (b x c)
	
	Vect3d a = p[2] - p[0];
	Vect3d b = p[1] - p[0];
	Vect3d c = p[3] - p[2];

	double ans = dot(a,cross(b,c));

	if (fabs(ans) < 0.0001)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Geometry::IsValidPolygon(std::vector<Vect3d> p)
{
	// Check if point cloud form a self-intersecting polygon
	// Check if shape is convex
	// Check orientation of points cw or ccw
	
	
}

// *******************************************
// End of Geometry Class


// *********************************************
// FinitePlane Class Derived from Geometry Class
FinitePlane::FinitePlane(std::vector<Vect3d> vertices_)
{
	if(IsCoplanar(vertices_))
	{
		vertices = vertices_;
		unit_norm = CalcUnitNorm(vertices);
		d = -1*(dot(unit_norm,vertices[0]));
	}
	else
	{
		//Throw an error
	}
}

FinitePlane::FinitePlane(Vect3d unit_norm_, double d_)
{
	unit_norm = unit_norm_;
	d = d_;
}

FinitePlane::~FinitePlane()
{
	d = 0;
}

bool FinitePlane::Intersects(Ray ray, Ray& ref_ray)
{
	// ALGORITHM: To find whether the ray interesects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray interesects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	//const double TH = 0.0001;
	
	const double TH = 1e-7;
	const double SMALL_DOUBLE = 1e-5;
	const double TWOPI = 2*M_PI;

	
	double t = 0;
	Vect3d p_int;
	double denom = dot(ray.dir,unit_norm);
	
	if (denom != 0)
	{
		t = (-1*(dot(unit_norm,ray.orig)+d))/denom;
	}
	if (t > 0)
	{			
		p_int = ray.orig + t*ray.dir;
		double sum_angles = 0;
		double sum_angles2 = 0;
		double mag1, mag2;
		Vect3d p1, p2;
		bool p_is_vertex = false;

		for (int i = 0; (i < vertices.size()) && !p_is_vertex; i++)
		{
			p1 = vertices[i] - p_int;
			p2 = vertices[(i+1)%vertices.size()] - p_int;

			mag1 = length(p1);
			mag2 = length(p2);
			if (mag1*mag2 <= SMALL_DOUBLE)
			{
				// This is true when the intersection point is one of the vertices
				// Consider it inside
				p_is_vertex = true;
			}
			else
			{
				bool result_radian = true;
				bool normalized = false;
				sum_angles += AngleBetween(p1,p2,normalized,result_radian);
			}
		}
		if (p_is_vertex || ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH)))
		{
			ref_ray.orig = p_int;
			ref_ray.dir = reflect(ray.dir, unit_norm);
			return true;
		}
	}
	ref_ray.orig = make_Vect3d(0,0,0);
	ref_ray.dir = make_Vect3d(0,0,0);
	return false;
}
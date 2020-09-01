#include "Geometry.hpp"

// *********************************************
// Geometry Base Class

Vect3d Geometry::CalcUnitNorm(const std::vector<Vect3d>& p)
{
	// Need to check if the vertices are clockwise or counter-clockwise
	// The direction of the unit normal (inward or outward pointing) depends on the orientation of the vertices
	
	Vect3d n1 = p[1] - p[0];
	Vect3d n2 = p[2] - p[1];
	
	Vect3d n = cross(n1,n2);
	
	n = normalize(n);
	
	return n;
}

bool Geometry::IsCoplanar(const std::vector<Vect3d>& p)
{
	// Source: https://mathworld.wolfram.com/Coplanar.html
	// Four points are coplanar if
	// (p2 - p0) DOT [(p1-p0) x (p3-p2)] = 0
	
	Vect3d a = p[2] - p[0];
	Vect3d b = p[1] - p[0];
	Vect3d c = p[3] - p[2];

	double ans = dot(a,cross(b,c));

	if (fabs(ans) < SMALL_DOUBLE)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Geometry::IsValidPolygon(const std::vector<Vect3d>& p)
{
	// Check if point cloud forms a self-intersecting polygon
	// Check if shape is convex
	// Check orientation of points cw or ccw
}

// *********************************************

// *********************************************
// FinitePlane Class Derived from Geometry Class
FinitePlane::FinitePlane(std::vector<Vect3d> vertices)
{
	if(IsCoplanar(vertices))
	{
		// Source: https://mathworld.wolfram.com/Plane.html
		vertices_ = vertices;
		unit_norm_ = CalcUnitNorm(vertices_);
		d_ = -1*(dot(unit_norm_,vertices_[0]));
	}
	else
	{
		//Throw an error
	}
}

FinitePlane::FinitePlane(std::vector<Vect3d> vertices, Vect3d unit_norm, double d)
{
	vertices_ = vertices;
	unit_norm_ = unit_norm;
	d_ = d;
}

FinitePlane::~FinitePlane()
{
	d_ = LARGE_DOUBLE;
}

bool FinitePlane::Intersects(const Ray& inc_ray, Ray& ref_ray)
{
	// ALGORITHM: To find whether the ray intersects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray intersects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	
	const double TWOPI = 2*M_PI;

	double t = 0;
	Vect3d p_int;
	double denom = dot(inc_ray.dir,unit_norm_);
	
	if (denom != 0)
	{
		t = (-1*(dot(unit_norm_,inc_ray.orig)+d_))/denom;
	}
	if (t > 0)
	{			
		p_int = inc_ray.orig + t*inc_ray.dir;
		double sum_angles = 0;
		double mag1, mag2;
		Vect3d p1, p2;
		bool p_is_vertex = false;

		for (int i = 0; (i < vertices_.size()) && !p_is_vertex; i++)
		{
			p1 = vertices_[i] - p_int;
			p2 = vertices_[(i+1)%vertices_.size()] - p_int;

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
		if (p_is_vertex || ((sum_angles < TWOPI + SMALL_DOUBLE) && (sum_angles > TWOPI - SMALL_DOUBLE)))
		{
			ref_ray.orig = p_int;
			// Reflect function could have an error if the direction of the unit normal (inward or outward) is not defined
			ref_ray.dir = reflect(inc_ray.dir, unit_norm_);
			return true;
		}
	}
	ref_ray.orig = make_Vect3d(LARGE_DOUBLE,LARGE_DOUBLE,LARGE_DOUBLE);
	ref_ray.dir = make_Vect3d(LARGE_DOUBLE,LARGE_DOUBLE,LARGE_DOUBLE);
	return false;
}


Ray FinitePlane::Compute_Reflect_Ray(const Ray& incident_ray)
{
	// Future implementation
	// Change Intersects to only return the point of intersection
	// Compute_Reflect_Ray will calculate the reflected ray
	//if(Intersects(incident_ray) == true)
	//return;
}


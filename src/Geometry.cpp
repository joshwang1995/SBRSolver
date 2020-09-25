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

/**
  * Checks if test point "p_test" lies within the boundaries of the polygon formed by the vertices "v" 
  *     This test assumes the vertices in vector "v" forms a non-self-intersecting convex polygon
  *     
  *
  * @param v - vertices of the polygon
  *        p_test - test point
  *     
  * @return true if test point is within the polygon, false if it is not or the test point is the one of the vertices
*/
bool Geometry::IsPointInPolygon(const std::vector<Vect3d>& v, const Vect3d& p_test)
{
    double sum_angles = 0;
    double mag1, mag2;
    Vect3d p1, p2;
    bool p_is_vertex = false;

    for (int i = 0; (i < v.size()) && !p_is_vertex; i++)
    {
        p1 = v[i] - p_test;
        p2 = v[(i+1)%v.size()] - p_test;

        mag1 = length(p1);
        mag2 = length(p2);
        if (mag1*mag2 <= SMALL_DOUBLE)
        {
            // If the magnitude of either p1 or p2 = 0, that means the point is the vertex
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
        return true;
    }
    else
    {
        return false;
    }
    
}
// *********************************************

// *********************************************
// FinitePlane Class Derived from Geometry Class
FinitePlane::FinitePlane(std::vector<Vect3d> v_vect)
{
	if(IsCoplanar(v_vect))
	{
		// Source: https://mathworld.wolfram.com/Plane.html
		vertices = v_vect;
		unit_norm = CalcUnitNorm(v_vect);
		d = -1*(dot(unit_norm,v_vect[0]));
	}
	else
	{
		//Throw an error
	}
}

FinitePlane::FinitePlane(std::vector<Vect3d> v_vect, Vect3d n_hat, double d_)
{
	vertices = v_vect;
	unit_norm = n_hat;
	d = d_;
}

/**
  * Checks if test point "p_test" lies within the boundaries of the polygon specified by the vertices "v" 
  *     This test assumes the vertices in vector "v" forms a non-self-intersecting convex polygon 
  *
  * @param v - vertices of the polygon
  *        p_test - test point
  *     
  * @return true if test point is within the polygon, false if it is not or the test point is the one of the vertices
*/
bool FinitePlane::Intersects(const Ray& inc_ray, Vect3d& p_int)
{
	// ALGORITHM: To find whether the ray intersects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray intersects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.

	double t = 0;
	Vect3d p_i;
	double denom = dot(inc_ray.dir,unit_norm);
	
	if (denom != 0)
	{
		t = (-1*(dot(unit_norm,inc_ray.orig)+d))/denom;
	}
    
	if (t > 0)
	{
		p_i = inc_ray.orig + inc_ray.dir*t;
        if(IsPointInPolygon(vertices,p_i) == true)
        {
			p_int = p_i;
            return true;
        }
	}
	p_int = Vect3d(LARGE_DOUBLE,LARGE_DOUBLE,LARGE_DOUBLE);
	return false;
}

#pragma once
#include "IcosahedronMesh.hpp"
#include "Vect_Utility.hpp"
#include <vector>
#include <iostream>

using namespace std;

struct SpherePoint
{
	double r;
	double theta;
	double phi;
};

struct RectCoord
{
	Vect3d ax;
	Vect3d ay;
	Vect3d az;
};

struct SphereCoord
{
	Vect3d ar;
	Vect3d at;
	Vect3d ap;
};

//Output and debug functions
void print_rays(vector<pair<double,double>>& ray_dir);
void print_rays(vector<pair<double,double>>& ray_dir, vector<vector<int>>& ray_tubes);

//Generate rays and ray tubes
vector<pair<double,double>> GetRayDir(vector<Vect3f> rays);
vector<vector<int>> ConstructRayTubes(vector<Vect3u>& faces, int num_vertices);
vector<Vect3f> GetRayVectOnIcosahedron(Vect3f originPoint, vector<Vect3f>& vertices, vector<Vect3u>& faces);
vector<pair<double,double>> GetLaunchAngle(vector<Vect3f>& ray_vect);

//Coordinate transformation
SpherePoint RectPoint2Sphere(Vect3f n);
Vect3d SpherePoint2Rect(SpherePoint sph_point);
double keep_between_zero_and_2pi(double angle);
pair<double,double> calc_launch_angle (float x, float y, float z);

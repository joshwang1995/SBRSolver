#include "IcosahedronMesh.hpp"
#include "Vect_Utility.hpp"
#include "Geometry.hpp"
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

//Generate rays and ray tubes
vector<vector<int>> ConstructRayTubes(vector<Vect3u>& faces, int num_vertices);
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, vector<Vect3f>& vertices, vector<Vect3u>& faces);
vector<pair<double,double>> GetLaunchAngle(vector<Vect3d>& ray_vect);

//Output and debug functions
void print_ray(Ray& r);
void print_ray_vect(vector<Ray> &r_vect);
void print_rays(vector<pair<double,double>>& ray_dir, vector<vector<int>>& ray_tubes);

//Coordinate transformation
SpherePoint RectPoint2Sphere(Vect3d n);
Vect3d SpherePoint2Rect(SpherePoint sph_point);
double keep_between_zero_and_2pi(double angle);
pair<double,double> calc_launch_angle (float x, float y, float z);

int main(int argc, char** argv)
{
	float x,y,z;
	float theta,phi;
	int subdivisions=0;
	bool useRayTubes = false;
	
	Vect3d origin = {0,0,0};
	
	//Construct regular icosahedron (no tessellation), hence subdivision = 0
	auto icosahedron = ConstructIcosahedron(0);
	vector<Vect3f> vertices = icosahedron.first;
	vector<Vect3u> faces = icosahedron.second;
	
	vector<Ray> ray_vects = GetRaysOnIcosahedron(origin, vertices, faces);
	//vector<vector<int>> ray_tubes = ConstructRayTubes(faces, vertices.size());
	
	//Create a square plane(4*4) that is 2.5m away from center of icosahedron
	Vect3d p0 = make_Vect3d(2,2.5,2);
	Vect3d p1 = make_Vect3d(-2,2.5,2);
	Vect3d p2 = make_Vect3d(-2,2.5,-2);
	Vect3d p3 = make_Vect3d(2,2.5,-2);
	vector<Vect3d> p {p0,p1,p2,p3};
	
	FinitePlane fp = FinitePlane(p);
	
	vector<Ray> ref_rays;
	vector<Vect3u> subdivide_faces;
	
	//Find all rays from regular icosahedron that intersects the line/plane
	int j = 0;
	while(j<4)
	{
		cout<< "Ray vects size: " << ray_vects.size() << "\n";
		for(int i=0; i<ray_vects.size(); i++)
		{
			Ray ref_ray;
			if(fp.Intersects(ray_vects[i], ref_ray))
			{
				ref_rays.push_back(ref_ray);
				subdivide_faces.push_back(faces[i]);
			}
		}
		subdivide_faces = Subdivide(vertices, subdivide_faces);
		ray_vects = GetRaysOnIcosahedron(origin, vertices, subdivide_faces);
		j++;
	}
}

//*************************************************************************************************************
// Exposed Functions
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, vector<Vect3f>& vertices, vector<Vect3u>& faces)
{
	vector<Ray> outbound_rays;
	
	for (int i = 0; i < faces.size(); i++)
	{
		//Compute the cente of the triangle
        double center_x = (vertices[faces[i].v0].x + vertices[faces[i].v1].x + vertices[faces[i].v2].x) / 3.0;
        double center_y = (vertices[faces[i].v0].y + vertices[faces[i].v1].y + vertices[faces[i].v2].y) / 3.0;
        double center_z = (vertices[faces[i].v0].z + vertices[faces[i].v1].z + vertices[faces[i].v2].z) / 3.0;
		Vect3d center_Triangle = make_Vect3d(center_x, center_y, center_z);
		Ray ray {originPoint, normalize(center_Triangle - originPoint)};
		outbound_rays.push_back(ray);
	}
	
	return outbound_rays;
}

vector<vector<int>> ConstructRayTubes(vector<Vect3u>& faces, int num_vertices)
{
	vector<vector<int>> ray_tubes;
	
	for(int v = 0; v<num_vertices; v++)
	{
		vector<int> ray_group;
		for(vector<Vect3u>::iterator it = faces.begin() ; it != faces.end(); ++it)
		{
			//Each Vect3u in Faces vector contains 3 indices to the vertices vector
			// Example: <0,4,1> means vertices[0], vertices[4] and vertices[1] form a triangle
			if(it->v0 == v || it->v1 == v || it->v2 == v)
			{
				ray_group.push_back(it - faces.begin());
			}
			
			if(v<=11 && ray_group.size() == 5)
			{
				//Vertices 0-11 are Original Vertex, they have pentagonal wavefronts
				//When all 5 rays are found, break out of for loop
				ray_tubes.push_back(ray_group);
				break;
			}
			else if(v>11 && ray_group.size() == 6)
			{
				//Vertices 12 and above are either edge or interior vertex, they have hexagonal wavefronts
				//When all 6 rays are found, break out of for loop
				ray_tubes.push_back(ray_group);
				break;	
			}
		}
	}
	return ray_tubes;
}

vector<pair<double,double>> GetLaunchAngle(vector<Vect3d>& ray_vect)
{
	vector<pair<double,double>> result;
	for(Vect3d ray: ray_vect)
	{
		pair<double,double> launch_angle = calc_launch_angle(ray.x,ray.y,ray.z);
		result.push_back(launch_angle);
	}
	return result;
}
// End of Exposed Functions
//*************************************************************************************************************


//*********************************************************
//--------------- HELPER FUNCTIONS BEGIN ------------------
SpherePoint RectPoint2Sphere(Vect3d n)
{
	SpherePoint result;
	double r, theta, phi;
	
	theta = acos(n.z/sqrt(n.x*n.x + n.y*n.y + n.z*n.z));
	phi = atan2(n.y, n.x);
	phi = keep_between_zero_and_2pi(phi);	
	
	//Convert answer to degrees
	theta = theta * (180.0/M_PI);
	phi = phi * (180.0/M_PI);
	
	result.r = sqrt(n.x*n.x + n.y*n.y +n.z*n.z);
	result.theta = theta;
	result.phi = phi;
	return result;
}

Vect3d SpherePoint2Rect(SpherePoint n)
{
	Vect3d result;
	result.x = n.r*sin(n.theta)*cos(n.phi);
	result.y = n.r*sin(n.theta)*sin(n.phi);
	result.z = n.r*sin(n.theta);
	return result;
}

pair<double,double> calc_launch_angle (float x, float y, float z)
{
	pair<double,double> launch_angle; 
	
	double theta, phi;
	theta = acos(z/sqrt(x*x + y*y + z*z));
	phi = atan2(y, x);
	phi = keep_between_zero_and_2pi(phi);	
	
	//Convert answer to degrees
	theta = theta * (180.0/M_PI);
	phi = phi * (180.0/M_PI);
	
	launch_angle.first = theta;
	launch_angle.second = phi;
	
	return launch_angle;
}

double keep_between_zero_and_2pi(double angle)
{
	return fmod(angle,2*M_PI) + (fmod(angle,2*M_PI) >= 0 ? 0 : 2*M_PI);			// This makes sure that the value is in [0:2*PI)
}

//--------------- HELPER FUNCTIONS END --------------------
//*********************************************************

//**************************************************
//--------------- DEBUGGING BEGIN ------------------
void print_ray(Ray& r)
{
	cout << "Ray Equation: <" << r.orig.x << "," << r.orig.y << "," << r.orig.z << "> + t";
	cout << "<" << r.dir.x << "," << r.dir.y  << "," << r.dir.z << ">\n";
}

void print_ray_vect(vector<Ray> &r_vect)
{
	for(Ray r: r_vect)
	{
		print_ray(r);
	}
}

void print_rays(vector<pair<double,double>>& ray_dir, vector<vector<int>>& ray_tubes)
{	
	int i = 0;
	
	for(vector<int> ray_tube: ray_tubes)
	{
		bool orig_vertex = (ray_tube.size() == 5)? true:false;
		cout << "Ray tube " << i << ": " ;
		
		if(orig_vertex)
		{
			cout << "Original Vertex" << "\n";
		}
		else
		{
			cout << "Edge or interior vertex" << "\n";
		}
		
		for(int j: ray_tube)
		{
			cout << "\t(Theta,Phi): (" << ray_dir[j].first << "," << ray_dir[j].second << ")\n";		
		}
		i++;
	}
}

//----------------- DEBUGGING END ------------------
//**************************************************
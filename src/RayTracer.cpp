#include "IcosahedronMesh.hpp"
#include "Vect_Utility.hpp"
#include "Geometry.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>

using namespace std;

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

/*
********************
Function Prototypes
*******************
*/

//Generate rays and ray tubes
vector<vector<int>> ConstructRayTubes(const vector<Vect3u>& faces, int num_vertices);
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, const vector<Vect3f>& vertices, const vector<Vect3u>& faces);

//Output and debug functions

//Coordinate transformation
RectCoord axes_rotation_matrix(const RectCoord& src, const RectCoord& target);
Vect3d rotate_point(const RectCoord& coord_sys, const Vect3d& p);

//Others
double keep_between_zero_and_2pi(double angle);

/*
****************************
End of Function Prototypes
****************************
*/

//**********************************
//* Global Variables and Constants
const double SPEED_OF_LIGHT = 2.99792458e8;
double freq = 2.4e9;
double lamda = SPEED_OF_LIGHT/freq;
double k = (2*M_PI)/lamda;
RectCoord global_coord {make_Vect3d(1,0,0), make_Vect3d(0,1,0), make_Vect3d(0,0,1)};
Vect3d tx_origin = {4,4,3}; //In global coordinate system
int tessellation=1;
//****************************

int main(int argc, char** argv)
{
	//Construct regular icosahedron (no tessellation), hence subdivision = 0
	auto icosahedron = ConstructIcosahedron(tessellation);
	vector<Vect3f> vertices = icosahedron.first;
	vector<Vect3u> faces = icosahedron.second;
	
	vector<Ray> ray_vects = GetRaysOnIcosahedron(tx_origin, vertices, faces);
	//vector<vector<int>> ray_tubes = ConstructRayTubes(faces, vertices.size());
	
	//Create a square plane(4*4) that is 2.5m away from center of icosahedron
	Vect3d p0 = make_Vect3d(6,2,1);
	Vect3d p1 = make_Vect3d(6,6,1);
	Vect3d p2 = make_Vect3d(2,6,1);
	Vect3d p3 = make_Vect3d(2,2,1);
	vector<Vect3d> p {p0,p1,p2,p3};
	FinitePlane fp = FinitePlane(p);
	
	// Create transmitter coordinate system
	RectCoord tx_coord {make_Vect3d(0,0,-1), make_Vect3d(1,0,0), make_Vect3d(0,-1,0)};
	RectCoord Rgt = axes_rotation_matrix(global_coord,tx_coord);
	RectCoord Rtg = axes_rotation_matrix(tx_coord,global_coord);
	
	complex<double> propagation_term;
	double rt, theta_t, phi_t;
	double st,sp,ct,cp;
	Ray ref_ray;
	
	for (Ray ray: ray_vects)
	{
		if(fp.Intersects(ray, ref_ray))
		{
			
            Vect3d coord_global = ref_ray.orig - ray.orig;
			Vect3d coord_tx = rotate_point(Rgt,coord_global);
            Vect3d coord_tx_sph = toSphereVect(coord_tx);
            
            rt = coord_tx_sph.x;
            theta_t = coord_tx_sph.y;
            phi_t = coord_tx_sph.z;
			
			
			propagation_term.real(cos(-1*k*rt)/rt);
			propagation_term.imag(sin(-1*k*rt)/rt);
			
			// Field of Hertzian Dipole in TX Spherical Coordinate System
			//complex<double> E_r = 0;
			//complex<double> E_theta = sin(theta_t)*propagation_term;
			//complex<double> E_phi = 0;
			
			// Field of half wave dipole in TX SphereCoord
			double pt = 1.0;
			complex<double> E_r = 0;
			complex<double> E_theta = sqrt(60*pt)*(cos(M_PI*cos(theta_t)/2.0)/sin(theta_t))*propagation_term;
			complex<double> E_phi = 0;
			
			// Field in TX Cartesian Coordinate System
			st = sin(theta_t); 
			sp = sin(phi_t) ; 
			ct = cos(theta_t); 
			cp = cos(phi_t) ; 
			
			RectCoord Rsc {make_Vect3d(st*cp,ct*cp,-1*sp),make_Vect3d(st*sp,ct*sp,cp),make_Vect3d(ct,-1*st,0)};
			
			complex<double> Ext = Rsc.ax.x*E_r + Rsc.ax.y*E_theta + Rsc.ax.z*E_phi;
			complex<double> Eyt = Rsc.ay.x*E_r + Rsc.ay.y*E_theta + Rsc.ay.z*E_phi;
			complex<double> Ezt = Rsc.az.x*E_r + Rsc.az.y*E_theta + Rsc.az.z*E_phi;
			
			complex<double> Ex = Rtg.ax.x*Ext + Rtg.ax.y*Eyt + Rtg.ax.z*Ezt;
			complex<double> Ey = Rtg.ay.x*Ext + Rtg.ay.y*Eyt + Rtg.ay.z*Ezt;
			complex<double> Ez = Rtg.az.x*Ext + Rtg.az.y*Eyt + Rtg.az.z*Ezt;
			
			
			cout << "*********************************************************\n";
			cout << "POI: (" << ref_ray.orig.x << "," << ref_ray.orig.y << "," << ref_ray.orig.z << ")\n";
			cout << "coordT: " << coord_tx.x << "," << coord_tx.y << "," << coord_tx.z << "\n";
			cout << "R, Theta, Phi: " << rt << "," << theta_t << "," << phi_t << "\n";
			cout << "E_theta: " << real(E_theta) << " + " << imag(E_theta) << "j\t" << "Mag: " << abs(E_theta) << "\t" << "Phase: " << arg(E_theta) << "\n";
			cout << "Rsc: \n(" << Rsc.ax.x << "," << Rsc.ax.y << "," << Rsc.ax.z << ")\n(" << Rsc.ay.x << "," << Rsc.ay.y << "," << Rsc.ay.z << ")\n(" << Rsc.az.x << "," << Rsc.az.y << "," << Rsc.az.z << ")\n";
			cout << "k0: " << k << "\n";
			cout << "Propagation Term: " << real(propagation_term) << " + " << imag(propagation_term) << "j\t" << "Mag: " << abs(propagation_term) << "\t" << "Phase: " << arg(propagation_term) << "\n";
			cout << "Ex: " << real(Ex) << " + " << imag(Ex) << "j\t" << "Mag: " << abs(Ex) << "\t" << "Phase: " << arg(Ex) << "\n";
			cout << "Ey: " << real(Ey) << " + " << imag(Ey) << "j\t" << "Mag: " << abs(Ey) << "\t" << "Phase: " << arg(Ey) << "\n";
			cout << "Ez: " << real(Ez) << " + " << imag(Ez) << "j\t" << "Mag: " << abs(Ez) << "\t" << "Phase: " << arg(Ez) << "\n";
            

            
			//complex<double> e_theta = 3*sin(theta)*propagation_term;
			//cout << "Efield Theta: " << e_theta << "\n";
		}
	}
}


/*
********************
Exposed Functions
********************
*/
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, const vector<Vect3f>& vertices, const vector<Vect3u>& faces)
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

vector<vector<int>> ConstructRayTubes(const vector<Vect3u>& faces, int num_vertices)
{
	vector<vector<int>> ray_tubes;
	
	for(int v = 0; v<num_vertices; v++)
	{
		vector<int> ray_group;
		for(vector<Vect3u>::const_iterator it = faces.cbegin() ; it != faces.cend(); ++it)
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
/*
********************
End of Exposed Functions
********************
*/

/*
********************
Helper Functions
********************
*/

double keep_between_zero_and_2pi(double angle)
{
	return fmod(angle,2*M_PI) + (fmod(angle,2*M_PI) >= 0 ? 0 : 2*M_PI);	// This makes sure that the value is in [0:2*PI)
}

//Returns the rotation matrix used in transforming local to global coordinate system and vice versa
RectCoord axes_rotation_matrix(const RectCoord& src, const RectCoord& target)
{
    RectCoord result_matrix;
	
    // Source coordinate system
    // ux = [xu yu zu]
    // uy = [xu yu zu]
    // uz = [xu yu zu]
    
    // Target coordinate system
    // wx = [xw yw zw]
    // wy = [xw yw zw]
    // wz = [xw yw zw]
    
    // Source to target coord_transform_matrix
    // [dot(ux,wx) dot(uy,wx) dot(uz,wx)]
    // [dot(ux,wy) dot(uy,wy) dot(uz,wy)]
    // [dot(ux,wz) dot(uy,wz) dot(uz,wz)]
    
    result_matrix.ax = make_Vect3d(dot(src.ax,target.ax), dot(src.ay,target.ax), dot(src.az,target.ax));
    result_matrix.ay = make_Vect3d(dot(src.ax,target.ay), dot(src.ay,target.ay), dot(src.az,target.ay));
    result_matrix.az = make_Vect3d(dot(src.ax,target.az), dot(src.ay,target.az), dot(src.az,target.az));
    return result_matrix;
}

// Gets the point in specified coordinate system
Vect3d rotate_point(const RectCoord& coord_sys, const Vect3d& p)
{
	Vect3d result;
	
    result.x = dot(coord_sys.ax, p);
    result.y = dot(coord_sys.ay, p);
    result.z = dot(coord_sys.az, p);
	return result;
}

/*
********************
End of Helper Functions
********************
*/

//**************************************************
//--------------- DEBUGGING BEGIN ------------------

//----------------- DEBUGGING END ------------------
//**************************************************
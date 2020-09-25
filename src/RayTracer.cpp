#include "IcosahedronMesh.hpp"
#include "Vect_Utility.hpp"
#include "Geometry.hpp"
#include <vector>
#include <iostream>
#include <complex>
#include <fstream>
#include "Constants.hpp"

using namespace std;

//Generate rays and ray tubes
vector<vector<Ray>> ConstructRayTubes(const vector<Ray>& ico_rays,const vector<Vect3u>& faces, int num_vertices);
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, const vector<Vect3f>& vertices, const vector<Vect3u>& faces);

Cvect3dsph ComputeAnalyticalEfield
(
	int antenna_type,
	double r,
	double theta,
	double phi,
	const RectCoord& global_coord_sys, 
	const RectCoord& local_coord_sys,
	double pt
);

// Functions to compute Fresnel reflection and transmission coefficients
void te_coeff
(
    complex<double> theta_i, complex<double> theta_t, 
    complex<double> rel_perm_i, complex<double> rel_perm_t,
    double lamda, double width, bool inf_wall,
    complex<double> &ref_coeff, complex<double> &tran_coeff
);

void tm_coeff
(
    complex<double> theta_i, complex<double> theta_t, 
    complex<double> rel_perm_i, complex<double> rel_perm_t,
    double lamda, double width, bool inf_wall,
    complex<double> &ref_coeff, complex<double> &tran_coeff
);

RectCoord GetSurfCoordSys(FinitePlane& wall, const Ray& inc_ray);
complex<double> trans_angle (double theta_i, complex<double> e_i, complex<double> e_t);
complex<double> invsin (const complex<double> &z);

//Output and debug functions
void csv_writeln(ofstream& file, Vect3d point_int, Vect3<complex<double>> E);
void log_rays(const vector<Ray>& rays);
void log_ray_tubes(const vector<vector<Ray>>& ray_tube);

/*
****************************
End of Function Prototypes
****************************
*/

//**********************************
//* Global Variables and Constants
double freq = 2.4e9;
double lamda = SPEED_OF_LIGHT/freq;
double k = (2*PI)/lamda;
RectCoord global_coord {Vect3d(1,0,0), Vect3d(0,1,0), Vect3d(0,0,1)};
Vect3d tx_origin = {4,4,3}; //In global coordinate system
//****************************

int main(int argc, char** argv)
{
    int antenna_type = 0;
    int tessellation= -1;
    double pt = 0;
	
    // User prompts
    cout << "Available antenna type:\n";
    cout << "\t1 - Hertzian Dipole\n";
    cout << "\t2 - Half wave Dipole\n";
    
    do
    {
        cout << "Select an antenna type: ";
        cin >> antenna_type;
    }
    while(antenna_type != 1 && antenna_type != 2);
    
    if(antenna_type == 2)
    {
        cout << "Specify the transmitter power of the half-wave dipole antenna: ";
        cin >> pt;
    }
    
    cout << "Enter the TX icosahedron tessellation frequency: ";
    cin >> tessellation;
    
	//Construct icosahedron with speicified tessellation frequency
	auto icosahedron = ConstructIcosahedron(tessellation);
	vector<Vect3f> vertices = icosahedron.first;
	vector<Vect3u> faces = icosahedron.second;
	
	vector<Ray> icosahedron_rays = GetRaysOnIcosahedron(tx_origin, vertices, faces);
	
    // Create scene. In this case we just have a finite plane
	Vect3d p0 = Vect3d(6,2,1);
	Vect3d p1 = Vect3d(6,6,1);
	Vect3d p2 = Vect3d(2,6,1);
	Vect3d p3 = Vect3d(2,2,1);
	vector<Vect3d> p {p0,p1,p2,p3};
	FinitePlane fp = FinitePlane(p);
	
	// Create transmitter coordinate system
    
    cout << "Start logging efield to efield.csv\n";
    ofstream efield_csv;
    efield_csv.open("efield.csv");
    efield_csv << "X,Y,Z,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag\n";
	
	
	RectCoord tx_coord {Vect3d(0,0,-1), Vect3d(1,0,0), Vect3d(0,-1,0)};
	double transmitter_power = 1.0;
	
	for (Ray ray: icosahedron_rays)
	{
		Vect3d p_int;
		if (fp.Intersects(ray, p_int))
		{
			/* Computing the electric field for direct path */
			
			// Step 1: Get ray path vector in global coordinate system
			Vect3d ray_global = p_int - ray.orig;

			// Step 2: Rotate vector from global to local coordinate system
			Vect3d ray_local = VectRectToRect(global_coord,tx_coord,ray_global);

			// Step 3: Convert vector from local rectangular to local spherical coordinate system
			Vect3dSph ray_local_sph = PointRectToSph(ray_local);

			double r_i = ray_local_sph.r;
			double theta_i = ray_local_sph.theta;
			double phi_i = ray_local_sph.phi;
			
			// Step 4: Compute analytical efield in spherical coordinate system
			Cvect3dsph efield_sph = ComputeAnalyticalEfield(antenna_type, r_i, theta_i, phi_i, global_coord, tx_coord, transmitter_power);
			
			// Step 5: Convert efield from SCS to RCS
			Cvect3d efield_rect = VectSphToRect(efield_sph, theta_i, phi_i);
			
			// Step 6: Convert efield from LCS to GCS
			Cvect3d e_i = VectRectToRect(tx_coord, global_coord, efield_rect);
			
			
			/* Computing the electric field for indirect path */
			
			// Step 1: Convert efield to facet-fixed coordinate system
			RectCoord surf_coord_sys = GetSurfCoordSys(fp, ray);
			
			//Cvect3d  = VectRectToRect(global_coord,facet_fixed,e_i);
			
			
			
			
			
			//csv_writeln(efield_csv, p_int, e_total);
		}
	}
	efield_csv.close();
}

Cvect3d ComputeRefcField(const Cvect3d& efield_i, double rel_perm, double sigma, double freq, double theta_i)
{
	complex<double> e_t;
	e_t.real(rel_perm);
	e_t.imag(sigma/(E0*2*PI*freq));
    
    double lamda = SPEED_OF_LIGHT/freq;
	
	// Permittivity of air
	complex<double> e_i(1,0);
	
	complex<double> theta_t = trans_angle(theta_i, e_i, e_t);

	complex<double> te_refc, tm_refc, te_tranc, tm_tranc;
	tm_coeff (theta_i, theta_t, e_i, e_t, lamda, width, inf_wall, tm_refc, tm_tranc);
	te_coeff (theta_i, theta_t, e_i, e_t, lamda, width, inf_wall, te_refc, te_tranc);
    
    
}

//Cvect3d ComputeTransField(double r, )
//{
	
//}

Cvect3dsph ComputeAnalyticalEfield
(
	int antenna_type,
	double r,
	double theta,
	double phi,
	const RectCoord& global_coord_sys, 
	const RectCoord& local_coord_sys,
	double pt
)
{
	Cvect3dsph e_field_sph (0,0,0);
	
	complex<double> propagation_term;
	propagation_term.real(cos(-1*k*r)/r);
	propagation_term.imag(sin(-1*k*r)/r);
	
	if(antenna_type == 1)
	{
		// Field of Hertzian Dipole in TX Spherical Coordinate System
		e_field_sph.r = 0;
		e_field_sph.theta = sin(theta)*propagation_term;
		e_field_sph.phi = 0;
	}
	else if(antenna_type == 2)
	{
		// Field of half wave dipole in TX SphereCoord
		e_field_sph.r = 0;
		e_field_sph.theta = sqrt(60*pt)*(cos(PI*cos(theta)/2.0)/sin(theta))*propagation_term;
		e_field_sph.phi = 0;
	}
	
	return e_field_sph;
}

/**
  * Generates the rays from the center to each face of the icosahedron
  * 
  * @param originPoint - the center of the icosahedron
  *        [in] vertices - vector containing all of the vertices of the icosahedron (in no particular order)
  *        [in] faces - vector containing all of the faces of the icosahedron
  *
  * @return vector of Ray structs, each Ray struct represents a ray from center of the icosahedron to one of its faces
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
		Vect3d center_Triangle = Vect3d(center_x, center_y, center_z);
		Ray ray {originPoint, normalize(center_Triangle)};
		outbound_rays.push_back(ray);
	}
	
	return outbound_rays;
}

/**
  * Group 'icosahedron_rays' into ray tubes
  *
  * @param ico_rays - a set of rays fired from the center of the icosahedron
  *        faces - vector containing all of the faces of the icosahedron
  *        num_vertices - total number of vertices
  *
  * @return 2D vector of struct Ray where each vector is a single ray tube
  *
  * @note
  *    
  *     Algorithm: Iterate from 0 to the number of vertices. Find all of the face triangles that contain the current vertex index. 
  *         For example, triangle <1,5,3> and triangle <6,12,1> both contain vertex[1], they belong in the same ray tube. We assume
  *         that the first 12 vertices is based on a regular icosahedron (tessellation=0), and any subsequent vertices appended is 
  *         either an edge or interior vertex. This assumption is made because icosahedron is generated in a iterative mannar. 
  *         Vertices 0-11 are Original Vertex, they have pentagonal wavefronts, so when we find 5 rays for that ray tube, we exit the
  *         for loop. Similarly, vertices 12 and above are either edge or interior vertex, they have hexagonal wavefronts, so we break 
  *         out when we find 6 rays.
  *  
  * @example
  *
  *     Before calling the ConstructRayTubes function, a set of rays need to be generated via GetRaysOnIcosahedron
  *
  *         auto icosahedron = ConstructIcosahedron(1);
  *         vector<Ray> ray_vects = GetRaysOnIcosahedron(makeVect3d(0,0,0), icosahedron.first, icosahedron.second);
  *         vector<vector<Ray>> ray_tubes = ConstructRayTubes(icosahedron.second, icosahedron.first.size());
  *     
  *     The above example first creates an icosahedron with tessellation = 1. It then generates a set of rays from the center to each
  *      face of the icosahedron. Lastly, the generated rays are grouped into ray tubes
*/
vector<vector<Ray>> ConstructRayTubes(const vector<Ray>& ico_rays, const vector<Vect3u>& faces, int num_vertices)
{
	vector<vector<Ray>> ray_tubes;
	
    // Iterate through the indices of the vertices
	for(int v = 0; v<num_vertices; v++)
	{
		vector<Ray> ray_group;
		for(vector<Vect3u>::const_iterator it = faces.cbegin() ; it != faces.cend(); ++it)
		{
			//Each Vect3u in 'faces' vector contains 3 indices to the vertices vector
			// Example: <0,4,1> means vertices[0], vertices[4] and vertices[1] form a triangle
			if(it->v0 == v || it->v1 == v || it->v2 == v)
			{
				ray_group.push_back(ico_rays[it - faces.begin()]);
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

RectCoord GetSurfCoordSys(FinitePlane& wall, const Ray& inc_ray)
{
	// Fail-safe procedure: ensure both ray direction and wall norm
	// 	is unit vector
	Vect3d n = normalize(wall.unit_norm);
	Vect3d ray_dir = normalize(inc_ray.dir);
	
	double theta_i = acos(dot(ray_dir,n));
	
	//double theta_i = acos((zw*inc_ray.dir)/(mag(zw)*mag(inc_ray.dir)));
	if (theta_i > PI/2)
	{
		theta_i = PI - theta_i;
		
		// From Neeraj's ray tracer code:
		// If the angle between incident ray and the surface normal is > 90 degrees,
		//	this indicates that the surface normal is inward pointing. 
		// 	Reverse the surface normal and also reverse orientation of vertices
		
		//wall.unit_norm_ = -1*wall.unit_norm_;
		//wall.d = - wall.d;
	}
	
	// Assuming both incident ray dir and zw is unit vector
	
	Vect3d zw = n;
	Vect3d yw = cross(ray_dir,n)/sin(theta_i);
	Vect3d xw = cross(yw,zw);
	
	return RectCoord{xw, yw, zw};
}

inline complex<double> operator * (const complex<double> &c1, const double &d)
{
	complex<double> c(c1.real()*d, c1.imag()*d);	
	return c;
}

inline complex<double> operator * (const double &d, const complex<double> &c1)
{
	complex<double> c(c1.real()*d, c1.imag()*d);	
	return c;
}

complex<double> trans_angle (double theta_i, complex<double> e_i, complex<double> e_t)
{	
	complex<double> t_i (theta_i, 0);
	complex<double> t_t;					

	complex<double> n_i = sqrt(e_i);
	complex<double> n_t = sqrt(e_t);
	t_t = invsin(n_i*sin(t_i)/n_t);
	return t_t;
}

inline complex<double> invsin (const complex<double> &z)
{
	//http://mathworld.wolfram.com/InverseSine.html
	// asin(z) = -i ln[iz + sqrt(1 - z^2)]
	complex<double> i(0, 1);
	complex<double> one(1, 0);
	return -1.0*i*log(i*z + sqrt(one - z*z));	
}

void tm_coeff 
(
    complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
    double lamda, double width, bool inf_wall,
	complex<double> &ref_coeff, complex<double> &tran_coeff
)
{
	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);	
	complex<double> cos_i = cos(theta_i);
	complex<double> cos_t = cos(theta_t); 
    complex<double> sin_i = sin(theta_i);

	complex<double> gamma_tm = (n_t*cos_i - n_i*cos_t)/(n_i*cos_t + n_t*cos_i);
    
    if(inf_wall)
	{
		ref_coeff = gamma_tm;
		trans_coeff = 0;
	}
	else
	{
        complex<double> q = (TWOPI*width/lamda)* sqrt(n_t - (sin_i*sin_i));
        complex<double> phi_factor = exp(-1i*2.0*q);
        complex<double> denom = 1.0 - (gamma_tm*gamma_tm*phi_factor);
        ref_coeff = (gamma_tm*(1.0 - phi_factor)) / denom;
        trans_coeff = ((1.0 - gamma_tm*gamma_tm)*phi_factor)/ denom;
	}
}

void te_coeff 
(
    complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	double lamda, double width, bool inf_wall,
	complex<double> &ref_coeff, complex<double> &tran_coeff
)
{	
	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);		
	complex<double> cos_i = cos(theta_i);
	complex<double> cos_t = cos(theta_t);
    complex<double> sin_i = sin(theta_i);

	complex<double> gamma_te = ((n_i*cos_i) - (n_t*cos_t))-((n_i*cos_i) + (n_t*cos_t));
	
	if(inf_wall)
	{
		ref_coeff = gamma_te;
		trans_coeff = 0;
	}
	else
	{
        complex<double> q = (TWOPI*width/lamda)* sqrt(n_t - (sin_i*sin_i));
        complex<double> phi_factor = exp(-1i*2.0*q);
        complex<double> denom = 1.0 - (gamma_te*gamma_te*phi_factor);
        ref_coeff = (gamma_te*(1.0 - phi_factor)) / denom;
        trans_coeff = ((1.0 - gamma_te*gamma_te)*phi_factor)/ denom;
	}
}


/* Logging Functions */
void csv_writeln(ofstream& file, Vect3d point_int, Cvect3d E)
{
    file << point_int.x << "," << point_int.y << "," << point_int.z << ",";
    file << real(E.x) << "," << imag(E.x) << ",";
    file << real(E.y) << "," << imag(E.y) << ",";
    file << real(E.z) << "," << imag(E.z) << "\n";
}

void log_rays(const vector<Ray>& rays)
{
    ofstream log_file;
    log_file.open("rays.csv");
    log_file << "x0,y0,z0,xr,yr,zr\n";
    
    for(Ray r: rays)
    {
		log_file << r.orig.x << "," << r.orig.y << "," << r.orig.z << ",";
		log_file << r.dir.x << "," << r.dir.y << "," << r.dir.z << "\n";
    }
    log_file.close();
}

void log_ray_tubes(const vector<vector<Ray>>& ray_tube)
{
    ofstream log_file;
    log_file.open("../visualization/raytube.csv");
    log_file << "x0,y0,z0,xr,yr,zr,ray_tube\n";
    int idx_ray_group = 1;
    
    for(vector<Ray> ray_group: ray_tube)
    {
        for(Ray r: ray_group)
        {
            log_file << r.orig.x << "," << r.orig.y << "," << r.orig.z << ",";
            log_file << r.dir.x << "," << r.dir.y << "," << r.dir.z << ",";
            log_file << idx_ray_group << "\n";
        }
        idx_ray_group += 1;
    }
    log_file.close();
}
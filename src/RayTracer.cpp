#include "IcosahedronMesh.hpp"
#include "Vect_Utility.hpp"
#include "Geometry.hpp"
#include "Constants.hpp"
#include "Interpolation.hpp"
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

//Generate rays and ray tubes
vector<vector<Ray>> ConstructRayTubes(const vector<Ray>& ico_rays,const vector<Vect3u>& faces, int num_vertices);
vector<Ray> GetRaysOnIcosahedron(Vect3d originPoint, const vector<Vect3f>& vertices, const vector<Vect3u>& faces);

Cvect3dsph GetAnalyticEfieldPattern
(
	int antenna_type,
	double theta,
	double phi,
	double pt = 0.0
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

RectCoord GetSurfCoordSys(RectWall& surf, const Ray& inc_ray);

Cvect3dsph ComputeRefcField
(
	const Cvect3dsph& efield_i_sph,
	double rel_perm, 
	double sigma, 
	double freq, 
	double theta_i,
	double width, 
	bool inf_wall
);

Cvect3dsph ComputeTransField
(
	const Cvect3dsph& efield_i_sph,
	double rel_perm, 
	double sigma, 
	double freq, 
	double theta_i,
	double width, 
	bool inf_wall
);

complex<double> trans_angle (double theta_i, complex<double> e_i, complex<double> e_t);

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
Point3d tx_origin = {4,4,3}; //In global coordinate system
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
	
    // Create scene. In this case we just have finite planes
	vector<Point3d> p {Point3d(6,2,1),Point3d(6,6,1),Point3d(2,6,1),Point3d(2,2,1)};
	vector<Point3d> p2 {Point3d(6,2,2),Point3d(6,6,2),Point3d(2,6,2),Point3d(2,2,2)};
	Geometry_List Scene;
	Scene.add(make_shared<RectWall>(p, 0.0, 0.1, 5));
	Scene.add(make_shared<RectWall>(p2, 0.0, 0.1, 5));
	
	RectCoord tx_coord {Vect3d(0,0,-1), Vect3d(1,0,0), Vect3d(0,-1,0)};
	double transmitter_power = 1.0;
	
	// First pass, fire rays from icosahedron
	// Return all reflected and transmitted rays + fields
	
	RectCoord current_coord_sys = global_coord;
	RectCoord next_coord_sys = tx_coord;
	
	/*
	for (Ray ray: icosahedron_rays)
	{
		Vect3d p_int;
		if (surf_0.Intersects(ray, p_int))
		{
			// Step 1: Get ray path vector in global coordinate system
			Vect3d ray_global = p_int - ray.orig;

			// Step 2: Rotate vector from global to local coordinate system
			Vect3d ray_local = VectRectToRect(current_coord_sys,next_coord_sys,ray_global);
			current_coord_sys = next_coord_sys;
			
			// Step 3: Convert vector from local rectangular to local spherical coordinate system
			Vect3dSph ray_local_sph = PointRectToSph(ray_local);

			double r = ray_local_sph.r;
			double theta = ray_local_sph.theta;
			double phi = ray_local_sph.phi;
			
			// Step 4: Apply the analytical efield pattern in spherical coordinate system
			Cvect3dsph efield_pattern = GetAnalyticEfieldPattern(antenna_type, theta, phi, transmitter_power);
			
			complex<double> propagation_term;
			propagation_term.real(cos(-1*k*r)/r);
			propagation_term.imag(sin(-1*k*r)/r);
			Cvect3dsph efield_local_sph;
			efield_local_sph.r = efield_pattern.r * propagation_term;
			efield_local_sph.theta = efield_pattern.theta * propagation_term;
			efield_local_sph.phi = efield_pattern.phi * propagation_term;
			
			// Step 5: Convert efield from SCS to RCS
			Cvect3d efield_rect = VectSphToRect(efield_local_sph, theta, phi);
			
			// Step 6: Convert coord system to next facet-fixed coord sys
			next_coord_sys = GetSurfCoordSys(surf_0, ray);
			Cvect3d efield_incdnt = VectRectToRect(current_coord_sys, next_coord_sys, efield_rect);
			ray_local = VectRectToRect(current_coord_sys,next_coord_sys,ray_local);
			ray_local_sph = PointRectToSph(ray_local);
			
			r = ray_local_sph.r;
			theta = ray_local_sph.theta;
			phi = ray_local_sph.phi;
			
			Cvect3dsph efield_incdnt_sph = VectRectToSph(efield_incdnt,theta,phi);
			Cvect3dsph efield_ref_sph = ComputeRefcField(efield_incdnt_sph, surf_rel_perm, surf_sigma, freq, theta, surf_width, true);
			Cvect3dsph efield_trans_sph = ComputeTransField(efield_incdnt_sph, surf_rel_perm, surf_sigma, freq, theta, surf_width, true);
			
			//csv_writeln(efield_csv, p_int, e_total);
		}
	}
	efield_csv.close();
	*/
}

Cvect3dsph Execute_Ray_Tracing
(
	vector<Interaction> path,
	const Geometry_List& scene, 
	const Point3d& rx_loc, 
	vector<Ray>& ref_rays,
	int depth
)
{
	
}

Cvect3dsph ComputeRefcField
(
	const Cvect3dsph& efield_i_sph,
	double rel_perm, 
	double sigma, 
	double freq, 
	double theta_i,
	double width, 
	bool inf_wall
)
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
    
	Cvect3dsph refc_field;
	refc_field.r = efield_i_sph.r;
	refc_field.phi = efield_i_sph.phi * te_refc;
	refc_field.theta = efield_i_sph.theta* tm_refc;
	
	return refc_field;
}

Cvect3dsph ComputeTransField
(
	const Cvect3dsph& efield_i_sph,
	double rel_perm, 
	double sigma, 
	double freq, 
	double theta_i,
	double width, 
	bool inf_wall
)
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
    
	Cvect3dsph trans_field;
	trans_field.r = efield_i_sph.r;
	trans_field.phi = efield_i_sph.phi * te_tranc;
	trans_field.theta = efield_i_sph.theta* tm_tranc;
	
	return trans_field;
}

Cvect3dsph GetAnalyticEfieldPattern
(
	int antenna_type,
	double theta,
	double phi,
	double pt
)
{
	Cvect3dsph e_field_sph (0,0,0);
	
	if(antenna_type == 1)
	{
		// Field of Hertzian Dipole in TX Spherical Coordinate System
		e_field_sph.r = 0;
		e_field_sph.theta = sin(theta);
		e_field_sph.phi = 0;
	}
	else if(antenna_type == 2)
	{
		// Field of half wave dipole in TX Spherical Coordinate System
		e_field_sph.r = 0;
		e_field_sph.theta = sqrt(60*pt)*(cos(PI*cos(theta)/2.0)/sin(theta));
		e_field_sph.phi = 0;
	}
	else if(antenna_type == 3)
	{
		// Field of antenna array in TX Spherical Coordinate System
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

RectCoord GetSurfCoordSys(RectWall& wall, const Ray& inc_ray)
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
	t_t = asin(n_i*sin(t_i)/n_t);
	return t_t;
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
		tran_coeff = 0;
	}
	else
	{
        complex<double> q = (TWOPI*width/lamda)* sqrt(n_t - (sin_i*sin_i));
        complex<double> phi_factor = exp(-j*2.0*q);
        complex<double> denom = 1.0 - (gamma_tm*gamma_tm*phi_factor);
        ref_coeff = (gamma_tm*(1.0 - phi_factor)) / denom;
        tran_coeff = ((1.0 - gamma_tm*gamma_tm)*phi_factor)/ denom;
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
		tran_coeff = 0;
	}
	else
	{
        complex<double> q = (TWOPI*width/lamda)* sqrt(n_t - (sin_i*sin_i));
		complex<double> j (0.0,1.0);
        complex<double> phi_factor = exp(-j*2.0*q);
        complex<double> denom = 1.0 - (gamma_te*gamma_te*phi_factor);
        ref_coeff = (gamma_te*(1.0 - phi_factor)) / denom;
        tran_coeff = ((1.0 - gamma_te*gamma_te)*phi_factor)/ denom;
	}
}
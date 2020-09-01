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

//Coordinate transformation
RectCoord RotationMatrix(const RectCoord& src, const RectCoord& target);
Vect3d GetVectInNewCoordSys(const RectCoord& rotation_matrix, const Vect3d& p);
double keep_between_zero_and_2pi(double angle);

//Output and debug functions
void csv_writeln(ofstream& file, Vect3d point_int, complex<double> Ex, complex<double> Ey, complex<double> Ez);
void log_rays(vector<Ray> rays);
void log_ray_tubes(const vector<vector<Ray>>& ray_tube);

// Functions to compute Fresnel reflection and transmission coefficients
// theta_i, theta_t is in degrees
void te_coeff (double theta_i, double theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac);

// theta_i, theta_t is in degrees
void te_coeff (complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac);

// theta_i, theta_t is in degrees
void tm_coeff (double theta_i, double theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac);

// theta_i, theta_t is in degrees
void tm_coeff (complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac);

/*
****************************
End of Function Prototypes
****************************
*/

//**********************************
//* Global Variables and Constants
double freq = 2.4e9;
double lamda = SPEED_OF_LIGHT/freq;
double k = (2*M_PI)/lamda;
RectCoord global_coord {make_Vect3d(1,0,0), make_Vect3d(0,1,0), make_Vect3d(0,0,1)};
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
	vector<vector<Ray>> ray_tubes = ConstructRayTubes(icosahedron_rays,faces, vertices.size());
    log_ray_tubes(ray_tubes);
	
    // Create scene. In this case we just have a finite plane
	Vect3d p0 = make_Vect3d(6,2,1);
	Vect3d p1 = make_Vect3d(6,6,1);
	Vect3d p2 = make_Vect3d(2,6,1);
	Vect3d p3 = make_Vect3d(2,2,1);
	vector<Vect3d> p {p0,p1,p2,p3};
	FinitePlane fp = FinitePlane(p);
	
	// Create transmitter coordinate system
	RectCoord tx_coord {make_Vect3d(0,0,-1), make_Vect3d(1,0,0), make_Vect3d(0,-1,0)};
	RectCoord Rgt = RotationMatrix(global_coord,tx_coord);
	RectCoord Rtg = RotationMatrix(tx_coord,global_coord);
	
	complex<double> propagation_term;
	double rt, theta_t, phi_t;
	double st,sp,ct,cp;
	Ray ref_ray;
    
    cout << "Start wrting to CSV\n";
    ofstream efield_csv;
    efield_csv.open("../visualization/efield.csv");
    efield_csv << "X,Y,Z,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag\n";
	
	for (Ray ray: icosahedron_rays)
	{
		if(fp.Intersects(ray, ref_ray))
		{
            Vect3d coord_global = ref_ray.orig - ray.orig;
			Vect3d coord_tx = GetVectInNewCoordSys(Rgt,coord_global);
            Vect3d coord_tx_sph = toSphereVect(coord_tx);
            
            rt = coord_tx_sph.x;
            theta_t = coord_tx_sph.y;
            phi_t = coord_tx_sph.z;
			
			
			propagation_term.real(cos(-1*k*rt)/rt);
			propagation_term.imag(sin(-1*k*rt)/rt);
            
            complex<double> E_r = 0;
            complex<double> E_theta = 0;
            complex<double> E_phi = 0;
			
            if(antenna_type == 1)
            {
                // Field of Hertzian Dipole in TX Spherical Coordinate System
                E_r = 0;
                E_theta = sin(theta_t)*propagation_term;
                E_phi = 0;
            }
            else if(antenna_type == 2)
            {
                // Field of half wave dipole in TX SphereCoord
                E_r = 0;
                E_theta = sqrt(60*pt)*(cos(M_PI*cos(theta_t)/2.0)/sin(theta_t))*propagation_term;
                E_phi = 0;
            }

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
			
            csv_writeln(efield_csv, ref_ray.orig, Ex,Ey,Ez);
		}
	}
	efield_csv.close();
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
		Vect3d center_Triangle = make_Vect3d(center_x, center_y, center_z);
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

/**
  * Computes the rotataion matrix used for converting between global and local rectangular coordinate systems. 
  *     This function is strictly used for rectanglar coordinate system rotation
  * 
  * @param src - source rectangular coordinate system
  *        target - target rectangular coordinate system
  * @return the rotation matrix in rectanglar coordinate system
  * @note
  *
  *     Source coordinate system
  *         ux = [xu yu zu]
  *         uy = [xu yu zu]
  *         uz = [xu yu zu]
  *  
  *     Target coordinate system
  *         wx = [xw yw zw]
  *         wy = [xw yw zw]
  *         wz = [xw yw zw]
  *  
  *     Source to target coord_transform_matrix
  *         [dot(ux,wx) dot(uy,wx) dot(uz,wx)]
  *         [dot(ux,wy) dot(uy,wy) dot(uz,wy)]
  *         [dot(ux,wz) dot(uy,wz) dot(uz,wz)]
  *
  * @example
  *
  *     To convert a vector from the global coordinate system to local coordinate system, we can use the function as following
  *     
  *         Rgt = RotationMatrix(global_coordinate_system, local_coordinate_system)
  *         vect_local_coord = vect_global_coord * Rgt
  *     
  *     Note that the multiplication is matrix multiplication
  *         
*/
RectCoord RotationMatrix(const RectCoord& src, const RectCoord& target)
{
    RectCoord result_matrix;

    result_matrix.ax = make_Vect3d(dot(src.ax,target.ax), dot(src.ay,target.ax), dot(src.az,target.ax));
    result_matrix.ay = make_Vect3d(dot(src.ax,target.ay), dot(src.ay,target.ay), dot(src.az,target.ay));
    result_matrix.az = make_Vect3d(dot(src.ax,target.az), dot(src.ay,target.az), dot(src.az,target.az));
    return result_matrix;
}

/**
  * Transform a 3D vector from old to new coordinate system
  * 
  * @param rotation_matrix - matrix for transforming from old to new coordinate system
  *        p - vector in old coorindate system
  * @return vector in new coordinate system

  * @example
  *
  *     To convert a vector from the global coordinate system to local coordinate system, we can use the function as following
  *     
  *         Rgt = RotationMatrix(global_coordinate_system, local_coordinate_system)
  *         vect_local_coord = GetVectInNewCoordSys(Rgt,vect_global_coord);
  *         
*/
Vect3d GetVectInNewCoordSys(const RectCoord& rotation_matrix, const Vect3d& p)
{
	Vect3d result;
	
    result.x = dot(rotation_matrix.ax, p);
    result.y = dot(rotation_matrix.ay, p);
    result.z = dot(rotation_matrix.az, p);
	return result;
}

void csv_writeln(ofstream& file, Vect3d point_int, complex<double> Ex, complex<double> Ey, complex<double> Ez)
{
    file << point_int.x << "," << point_int.y << "," << point_int.z << ",";
    file << real(Ex) << "," << imag(Ex) << ",";
    file << real(Ey) << "," << imag(Ey) << ",";
    file << real(Ez) << "," << imag(Ez) << "\n";
}

void log_rays(const vector<Ray>& rays)
{
    ofstream log_file;
    log_file.open("../visualization/rays.csv");
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


void te_coeff (double theta_i, double theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac, double thresh)
{
	// Refer to page 429 of Cheng
	//          n_i*cos(theta_i) - n_t*cos(theta_t)      eta_t*cos(theta_i) - eta_i*cos(theta_t)
	//	r,TE = ------------------------------------  =  -----------------------------------------
	//          n_i*cos(theta_i) + n_t*cos(theta_t)      eta_t*cos(theta_i) + eta_i*cos(theta_t)
	//
	//                 2 n_i*cos(theta_i) 
	//	t,TE = -----------------------------------
	//          n_i*cos(theta_i) + n_t*cos(theta_t) 

	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);	
	double cos_i = cos(theta_i);
	double cos_t = cos(theta_t); 

	complex<double> denom = n_i*cos_i + n_t*cos_t;
	ref_coeff = (n_i*cos_i - n_t*cos_t)/denom;
	tran_coeff = (2*n_i*cos_i)/ denom;

	// Surface Roughness Factor
	double rho = exp(-8*pow(((M_PI*rho_h*cos(theta_i))/lamda), 2));
	// Apply the scattering/attenuation factor only if the magnitude of the reflection coefficient 
	// is larger than thresh
	if (abs(ref_coeff) > thresh)
	{
		ref_coeff = mod_fac*((1 + rho)/2)*ref_coeff;
	}
}

void te_coeff (complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac, double thresh)
{	
	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);		
	complex<double> cos_i = cos(theta_i);
	complex<double> cos_t = cos(theta_t);

	complex<double> denom = n_i*cos_i + n_t*cos_t;
	ref_coeff = (n_i*cos_i - n_t*cos_t)/denom;
	tran_coeff = (2*n_i*cos_i)/ denom;

	// Surface Roughness Factor
	double rho = exp(-8*pow(((M_PI*rho_h*cos(theta_i.real()))/lamda), 2));
	// Apply the scattering/attenuation factor only if the magnitude of the reflection coefficient 
	// is larger than thresh
	if (abs(ref_coeff) > thresh)
	{
		ref_coeff = mod_fac*((1 + rho)/2)*ref_coeff;
	}
}

void tm_coeff (double theta_i, double theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac, double thresh)
{
	// Refer to page 431 of Cheng
	//          n_i*cos(theta_t) - n_t*cos(theta_i)
	//	r,TM = ----------------------------------- 
	//          n_i*cos(theta_t) + n_t*cos(theta_i) 
	//
	//                 2 n_i*cos(theta_i) 
	//	t,TM = -----------------------------------
	//          n_i*cos(theta_t) + n_t*cos(theta_i) 
	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);	
	double cos_i = cos(theta_i);
	double cos_t = cos(theta_t); 

	complex<double> denom = n_i*cos_t + n_t*cos_i;
	ref_coeff = (n_t*cos_i - n_i*cos_t)/denom;
	tran_coeff = (2*n_i*cos_i)/ denom;

	// Surface Roughness Factor
	double rho = exp(-8*pow(((M_PI*rho_h*cos(theta_i))/lamda), 2));
	// Apply the scattering/attenuation factor only if the magnitude of the reflection coefficient 
	// is larger than thresh
	if (abs(ref_coeff) > thresh)
	{
		ref_coeff = mod_fac*((1 + rho)/2)*ref_coeff;
	}
}

void tm_coeff (complex<double> theta_i, complex<double> theta_t, 
	complex<double> rel_perm_i, complex<double> rel_perm_t,
	complex<double> &ref_coeff, complex<double> &tran_coeff,
	double rho_h, double lamda, double mod_fac, double thresh)
{
	complex<double> n_i, n_t;
	n_i = sqrt(rel_perm_i);
	n_t = sqrt(rel_perm_t);	
	complex<double> cos_i = cos(theta_i);
	complex<double> cos_t = cos(theta_t); 

	complex<double> denom = n_i*cos_t + n_t*cos_i;
	ref_coeff = (n_t*cos_i - n_i*cos_t)/denom;
	tran_coeff = (2*n_i*cos_i)/ denom;

	// Surface Roughness Factor
	double rho = exp(-8*pow(((M_PI*rho_h*cos(theta_i.real()))/lamda), 2));
	// Apply the scattering/attenuation factor only if the magnitude of the reflection coefficient 
	// is larger than thresh
	if (abs(ref_coeff) > thresh)
	{
		ref_coeff = mod_fac*((1 + rho)/2)*ref_coeff;
	}
	/*
	cout << "\n===============REFLECTION COEFFICIENT BEGIN====================\n\n";
	cout << "theta_i: " << theta_i << "\n";
	cout << "theta_t: " << theta_t << "\n";
	cout << "n_i: " << n_i << "\n";
	cout << "n_t: " << n_t << "\n";
	cout << "cos_i: " << cos_i << "\n";
	cout << "cos_t: " << cos_t << "\n";
	cout << "denom: " << denom << "\n";
	cout << "ref_coeff: " << ref_coeff << "\n";
	cout << "\n===============REFLECTION COEFFICIENT END====================\n";
	*/
}

double keep_between_zero_and_2pi(double angle)
{
	return fmod(angle,2*M_PI) + (fmod(angle,2*M_PI) >= 0 ? 0 : 2*M_PI);	// This makes sure that the value is in [0:2*PI)
}
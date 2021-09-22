#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <complex>
#include "def.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string>

#pragma warning( disable : 4996)
#pragma warning( disable : 4267)

using namespace std;

typedef complex<double> cdouble;

#define MAX_CMD_ARG_LEN 200
#define MAX_VAR_TARGET_LEN 50
#define MAX_NUM_LEN 150
// bool enable_illum_zone = true;
// int use_pattern = 1;

/*****************************************************************************************/
// Structures _BEGIN

struct point_list
{
	point p;
	point_list* next;
};

struct xcylinder
{
	double r;
	double yc;
	double zc;
};

struct cylinder
{
	double r;
	point p0;
	vect u;
};

// generic surface structure
// contains a superset of all the parameters
// used to specify all supported surfaces
struct surface
{
	int id;
	surface_type s_type;
	double d;	// d for plane (ax + by + cz + d = 0), radius of the sphere
	vect normal;	// (a, b, c) of the plane equation, not specified for sphere
	double rel_perm;
	double sigma;				// conductivity, E'' = sigma/(W*E0)
	double rho_h;				// standard deviation of the surface roughness
	int num_vert;				// number of vertices of a finite plane/polygon
	//point v[4];					// vertices of the finite plane
	double shortest_side;		// shortest side of the surface
	point_list* v;
	point c;					// center of the sphere
	cylinder cyl;				// stores parameters for the cylindrical surfaces that is being approximated by the planar surface
	bool planar;				// true => planar representation is exact
};

/*
struct surface_list
{
surface surf;
surface* next_surf;
};
*/

// Definition of a 3D plane
struct plane
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'
	double d;
	vect unit_normal;
	double rel_perm;
	double sigma;
	double rho_h;
};

/*
struct vert_list
{
point p;
vert_list* v;
};
*/

struct finite_plane
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'

	int id;
	double d;
	vect unit_normal;
	double rel_perm;
	// number of vertices (plane can have any number of vertices greater than 2)
	int num_vert;
	//point p1, p2, p3, p4;	
	// points must be specified in order
	// point vl [4];
	point_list* v;
	double sigma;
	double rho_h;
};

// NSOOD: Illumination zone _START

struct plane_lite // without Wall properties
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'
	double d;
	vect unit_normal;
};

struct finite_nplane_lite // without Wall properties
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'

	double d;
	vect unit_normal;
	//	int num_vert;
	point_list* vert; 	// points must be specified in order
};

struct finite_plane_lite // without Wall properties
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'

	double d;
	vect unit_normal;
	int num_vert;
	point_list* v;
	//	point vl [4]; 	// points must be specified in order
};

struct fp_lite_list_node
{
	finite_plane_lite fp;
	fp_lite_list_node* next;
};

struct finite_plane_lite3 // without Wall properties
{
	//  Equation of a plane is
	//	ax + by + cz + d = 0
	//	unit_normal represents (a,b,c)
	//	and d represents the constant term 'd'

	double d;
	vect unit_normal;
	int num_vert;
	point vl [3]; 	// points must be specified in order
};

struct rect_pyramid
{
	// Specifies the illumination zone that corresponds to an image source
	// Do not need to store the relative permittivity and conductivity
	// - The apex is the location of the image
	// - The vertices of the base are obtained by intersecting the line formed
	//   by the apex and vertices of the plane through which image is obtained
	//   and the bounding cube.
	// - NSOOD: Aug, 2015 update
	//	 - Don't need the bounding cube anymore
	//   - Updating the pyramid-plane intersection method
	//   - Also don't need to explicitly store faces as finite_plane_lite3 objects
	point apex;
	finite_plane_lite base;
	fp_lite_list_node* face_list;
};

struct cube	// actually a more general rectangular prism
{
	// Is made up of 6 planes (lite)
	finite_plane_lite finite_top;
	finite_plane_lite finite_bottom;
	finite_plane_lite finite_left;
	finite_plane_lite finite_right;
	finite_plane_lite finite_near;
	finite_plane_lite finite_far;
	plane_lite top;
	plane_lite bottom;
	plane_lite left;
	plane_lite right;
	plane_lite near;
	plane_lite far;
};

// NSOOD: Illumination zone _END

// Definition of a Sphere
struct sphere
{
	double r;
	point c;
	double rel_perm;
};

struct efield
{	
	vect dir;
	complex<double> mag;
};

struct c_vect
{
	complex<double> x;
	complex<double> y;
	complex<double> z;
};

struct c_vect2
{
	complex<double> theta;
	complex<double> phi;	
};

struct image_tree_node
{
	int surf_id;
	point image_loc;
	image_tree_node* parent;
	image_tree_node* child;
	image_tree_node* sibling;
};

struct path_node
{
	point p;	// Intersection point
	point ap;	// Adjusted intersection point
	point ip;	// Image point
	int surf_id;
	path_node* next;
	path_node* prev;
};

struct path_list_node
{
	path_node* path;	
	path_list_node* next;
	double path_length;
	double adj_path_length; // Adjusted path length
	string path_id;
	int num_refs;
	double tx_theta;
	double tx_phi;
	double rx_theta;
	double rx_phi;
};

// All local co-ordinate systems are defined in terms
// of the global cartesian co-ordinate system
struct rect_cord
{
	vect ax;
	vect ay;
	vect az;
};

struct sph_cord
{
	vect ar;
	vect at;
	vect ap;
};

struct point_sph
{
	double r;
	double theta;
	double phi;
};

// Structures _END
/*****************************************************************************************/
// Some useful helper functions for complex numbers
inline c_vect init ()
{
	c_vect v;
	v.x.real(0);
	v.x.imag(0);
	v.y.real(0);
	v.y.imag(0);
	v.z.real(0);
	v.z.imag(0);	
	return v;
}

inline c_vect2 init2 ()
{
	c_vect2 v;
	v.phi.real(0);
	v.phi.imag(0);
	v.theta.real(0);
	v.theta.imag(0);	
	return v;
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

inline complex<double> invsin (const complex<double> &z)
{
	//http://mathworld.wolfram.com/InverseSine.html
	// asin(z) = -i ln[iz + sqrt(1 - z^2)]
	complex<double> i(0, 1);
	complex<double> one(1, 0);
	return -1*i*log(i*z + sqrt(one - z*z));	
}

inline efield operator + (const efield &e1, const efield &e2)
{
	efield e;
	complex<double> x = e1.mag*e1.dir.x + e2.mag*e2.dir.x;
	complex<double> y = e1.mag*e1.dir.y + e2.mag*e2.dir.y;
	complex<double> z = e1.mag*e1.dir.z + e2.mag*e2.dir.z;

	// TODO: Read up on complex vectors and how their magnitude and direction is calculated.
	e.mag = sqrt(x*x + y*y + z*z);
	complex<double> x_dir = (x/e.mag);
	complex<double> y_dir = (y/e.mag);
	complex<double> z_dir = (z/e.mag);

	e.dir.x = x_dir.real();
	e.dir.y = y_dir.real();
	e.dir.z = z_dir.real();

	return e;
}

inline c_vect operator + (const c_vect &v1, const c_vect &v2)
{
	c_vect v;
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	v.z = v1.z + v2.z;
	return v;
}

inline c_vect operator - (const c_vect &v1, const c_vect &v2)
{
	c_vect v;
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;
	return v;
}

inline c_vect operator * (const double d, const c_vect &v1)
{
	c_vect v;
	v.x = d*v1.x;
	v.y = d*v1.y;
	v.z = d*v1.z;
	return v;
}

inline c_vect operator * (const c_vect &v1, const double d)
{
	c_vect v;
	v.x = d*v1.x;
	v.y = d*v1.y;
	v.z = d*v1.z;
	return v;
}

inline c_vect operator * (const vect &v1, const complex<double> c)
{
	c_vect v;
	v.x = c*v1.x;
	v.y = c*v1.y;
	v.z = c*v1.z;
	return v;
}

inline c_vect operator * (const complex<double> c, const vect &v1)
{
	c_vect v;
	v.x = c*v1.x;
	v.y = c*v1.y;
	v.z = c*v1.z;
	return v;
}

inline c_vect operator * (const c_vect &v1, const complex<double> c)
{
	c_vect v;
	v.x = c*v1.x;
	v.y = c*v1.y;
	v.z = c*v1.z;
	return v;
}

inline c_vect operator * (const complex<double> c, const c_vect &v1)
{
	c_vect v;
	v.x = c*v1.x;
	v.y = c*v1.y;
	v.z = c*v1.z;
	return v;
}

inline complex<double> operator * (const c_vect &v1, const c_vect &v2)
{
	complex<double> dot_product = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	return dot_product;
}

inline complex<double> operator * (const c_vect &v1, const vect &v2)
{
	complex<double> dot_product = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	return dot_product;
}

inline complex<double> operator * (const vect &v2, const c_vect &v1)
{
	complex<double> dot_product = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	return dot_product;
}

inline double sq_mag (const c_vect &v)
{
	double magnitude = abs(v.x)*abs(v.x) +
		abs(v.y)*abs(v.y) +
		abs(v.z)*abs(v.z);
	return magnitude;
}

void init_point 
	(
	point & p,
	double x,
	double y,
	double z
	)
{
	p.x = x;
	p.y = y;
	p.z = z;
}

/*****************************************************************************************/
// Printing a 3D point 
inline ostream& operator << (ostream &output, const point &p)
{
	output << p.x << ", " << p.y << ", " << p.z;
	return output;
}

inline ostream& operator << (ostream &output, const point_sph &p)
{
	output << p.r << ", " << p.phi << ", " << p.theta;
	return output;
}

// Printing a vector
inline ostream& operator << (ostream &output, const vect &v)
{
	output << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return output;
}

// Printing a vector with complex projections
inline ostream& operator << (ostream &output, const c_vect &v)
{
	/*	output << v.x << " x +\n" 
	<< v.y << " y +\n" 
	<< v.z << " z\n";
	*/
	output << v.x.real() << "\t" << v.x.imag() << "\t"
		<< v.y.real() << "\t" << v.y.imag() << "\t"
		<< v.z.real() << "\t" << v.z.imag() << "\t";
	return output;
}

// Printing a vector with complex projections
inline ostream& operator << (ostream &output, const c_vect2 &v)
{
	output << v.phi << " phi +\n" 
		<< v.theta << " theta\n";
	return output;
}

inline ostream& operator << (ostream &output, image_tree_node* iter)
{
	if (iter != NULL)
	{
		// Potential path id string
		char temp[4];
		//itoa(iter->surf_id, temp, 10);
        sprintf(temp, "%d", iter->surf_id);
		string id = "R-";
		if (iter->surf_id == -1)
		{
			id = id + "T";
		}
		else
		{
			id = id + temp;
		}

		image_tree_node* piter = iter->parent;
		while (piter != NULL)
		{
			//itoa(piter->surf_id, temp, 10);
            sprintf(temp, "%d", piter->surf_id);
            if (piter->surf_id == -1)
			{
				id = id + "-T";
			}
			else
			{
				id = id + "-" + temp;
			}			
			piter = piter->parent;
		}
		id = string(id.rbegin(),id.rend());
		cout << id << endl;
		output << iter->image_loc << "\n" << iter->surf_id;	
	}
	else
	{
		output << "Printing a NULL pointer\n";
	}
	return output;
}

/*****************************************************************************************/
// HELPER FUNCTIONS _BEGIN

bool is_line_comment (string s);
int get_num_surfaces (char* input_file_name);
//void get_surface_info ();
void string_to_point (string s, point& p);
void string_to_point (string s, vect& v);
//void string_to_vect (string s, vect& v);
void string_to_complex (string s, complex<double> & c);
void string_to_vertices (string s, point_list*& v, int& num_vert);
void string_to_surface_type(string param_value, surface_type &st);
void get_param_from_line (string &param_name, string &param_value, string s);
void get_params (point & trans_loc, point &rec_loc, 
	char* input_file_name, double &resolution, double &power,
	surface surf[], int &rec_level, double &start_freq,
	double &end_freq, double &incr_freq, string &sweep_param,
	int num_surfaces, int &image_tree_int, double &rel_perm_global,
	double &sigma_global, double &rho_h_global, int &dump_tree_int, 
	int &dump_paths_int, double &refcoef_fac, int &use_pattern, 
	double &patch_red_fac, double &bw, double &ref_thresh,
	string &tx_fname, string &rx_tcut_fname, string &rx_pcut_fname,
	int &rem_red_paths, int &dump_paths_info, double &length_interval,
	double &angle_interval, int &enable_illum_zone_int, int &adj_path_geom);

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

complex<double> trans_angle (double theta_i, complex<double> e_i, complex<double> e_t);

// Functions to compute intersections
bool intersection_ray_plane (plane p, ray r, ray &ref_ray, double &t_value);
bool intersection_ray_plane_lite (plane_lite p, ray r, double &t_value);
bool intersection_ray_finite_plane (finite_plane p, ray r, ray &ref_ray, double &t_value);
bool intersection_ray_finite_plane_lite (finite_plane p, ray r, double &t_value);
bool intersection_ray_finite_plane_lite (finite_plane_lite p, ray r, double &t_value);
bool intersection_ray_finite_plane_lite (finite_plane_lite3 p, ray r, double &t_value);
bool intersection_ray_inf_finite_plane_lite (finite_plane_lite p, ray r, double &t_value);
bool intersection_line_seg_line_seg (point p1, point p2, point p3, point p4, double& t_value, double& s_value);
//bool intersection_line_seg_triangle (point apex, point a, point b, ray r, double &t_value);
bool intersection_ray_sphere (sphere s, ray r, ray &ref_ray, double &t_value);
bool intersection_ray_xcylinder (xcylinder &xc, ray &r, double &t_value, point &ap);
bool intersection_ray_cylinder (cylinder &c, ray &r, double &t1, double &t2, point &p1, point &p2);
vect sphere_normal (point p, point center);
bool intersection_finite_plane_pyramid (finite_plane p, rect_pyramid py);
bool intersection_finite_plane_ext_zone (finite_plane p, rect_pyramid py);

vect theta_unit_vect(vect v);
void my_delete (path_node* node);
void my_delete (path_list_node* node);
void my_delete (image_tree_node* node);
void my_delete (point_list* node);
void my_delete (fp_lite_list_node* node);


void image_source_loc (finite_plane p, point src, point &img);
void create_image_tree (surface surf[], int num_surfs, int num_refls, 
	image_tree_node* root, cube c, bool enable_illum_zone);
void build_root_children (image_tree_node* parent, surface surf[], int num_surfs, int& num_node);
void build_node_children (image_tree_node* parent, surface surf[], int num_surfs,
	cube c, int& num_node, int& num_node_level, bool enable_illum_zone);
bool next_cousin (image_tree_node* iter, image_tree_node*& niter);

void dump_image_tree (image_tree_node* root);
finite_plane surf_to_fplane (surface surf);
finite_plane_lite surf_to_fplane_lite (surface surf);
void dump_image_tree_rec (image_tree_node* root);
void extract_paths (image_tree_node* root, surface surf[], int num_surfs, point rec_loc, 
	point trans_loc, path_list_node*& list_of_paths, int& num_paths);
void extract_paths_rec 
	(
	image_tree_node* curr_node,
	surface surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_list_node*& path_list_iter1,
	path_list_node*& path_list_iter2,
	int& path_count
	);
bool check_path (image_tree_node* node, surface surf[], int num_surfs, point rec_loc, 
	point trans_loc, path_node*& path);

void compute_path_lengths (path_list_node* path_list, point trans_loc, point rec_loc, rect_cord coord, surface surf[], int num_surfs);
void dump_paths (path_list_node* path_list);
c_vect compute_path_gain_vector_with_rotation (path_list_node* path_list, point trans_loc, point_sph trans_angles, point rec_loc, point_sph rec_angles, 
	surface surf[], int num_surfs, double k, double rho_h, double lamda,
	double freq, double trans_power, double theta_gain[], double phi_gain[], int refs_l, int refs_u, double refcoef_fac,
	int &use_pattern, double theta_gain_rx_tc[], double theta_gain_rx_pc[],
	double phi_gain_rx_tc[], double phi_gain_rx_pc[], double bw, double refcoef_thresh,
	rect_cord trans_cord, rect_cord rec_cord, int adj_path_geom);
c_vect2 efield_pattern (double theta, double phi, double trans_power, double theta_gain[], double phi_gain[]);
c_vect sph_to_rect (c_vect2 v_sph, double theta, double phi);
c_vect2 rect_to_sph (c_vect v_rect, double theta, double phi);
c_vect rect_to_rect (c_vect v_rect_old, rect_cord old_sys, rect_cord new_sys);
vect rect_to_rect (vect v_rect_old, rect_cord old_sys, rect_cord new_sys);
point rect_to_rect (point v_rect_old, rect_cord old_sys, rect_cord new_sys);
rect_cord get_local_cord (surface surf);
void get_global_cord (rect_cord &global);
void get_sph_ref (point rect_ref, point_sph &sph_ref);
double keep_between_zero_and_2pi(double angle);

//double gain (double antenna_pattern[], double theta, double phi, int freq_index);
double gain (double antenna_pattern[], double theta, double phi);
double gain_theta (double theta, double phi);
double gain_theta_rx (double theta, double phi, double theta_cut[], double phi_cut[]);

double delta (double a1, double b1, double a2, double b2);
int get_loc (string s, point &rec);
int get_loc_and_angles (string s, point &rec, point_sph &angles);
int get_ant_info (string str, point &loc, point_sph &angle, rect_cord &ant_cord);

void  read_pattern_files (double theta[], double phi[]);
void  read_pattern_files_nff (double theta[], double phi[], string fname, const int size);
void  read_pattern_for_freq (double gain_data[], string fname, int freq_index);
string freq_index_to_string (int freq_index);
int freq_to_freq_index (double freq);

void loop (int path [], int path_size, int num_surfs, int index,
	finite_plane surf[], point trans_loc, point rec_loc, path_list_node*& list_of_paths);
void gen_paths 
	(
	int path_nodes [],
	int num_path_nodes,
	finite_plane surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_list_node*& list_of_paths
	//	int& path_count
	);

bool check_direct_path (finite_plane surf[], int num_surfs, point trans_loc, point rec_loc);

void setup_cube 
	(
	cube &c,
	double xmin, 
	double xmax, 
	double ymin, 
	double ymax, 
	double zmin, 
	double zmax
	);

void init_finite_plane_lite
	(
	finite_plane_lite & fp,
	// Points pa, pb, pc, pd are in cyclic or acyclic order 
	point pa,
	point pb,
	point pc,
	point pd
	);

void init_finite_plane_lite
	(
	finite_plane_lite & fp,	
	// Points pa, pb, pc are in cyclic or acyclic order (always satisfied for 3 points)
	point pa,
	point pb,
	point pc
	);

void init_finite_plane_lite
	(
	finite_plane_lite3 & fp,
	// Points pa, pb, pc are in cyclic or acyclic order (always satisfied for 3 points)
	point pa,
	point pb,
	point pc
	);

void plane_normal (finite_plane_lite & fp);
void plane_normal (finite_plane_lite3 & fp);
bool iszero(double f);
void reduce_surface (surface & s, double fac);
double compute_att_fac(double a, double theta_i, double k, double bw, bool isTE);
double sinc(const double& x);
double shortest_side(surface surf);

bool compare_paths(path_list_node* path1, path_list_node* path2, double length_interval, double angle_interval);
bool compare_points(const point& p1, const point& p2);
void remove_duplicate_paths(path_list_node* iter, double length_interval, double angle_interval);
void dump_paths_info_func (path_list_node* path_list, point trans_loc, point rec_loc);
void compute_departure_angles(vect ray_dir, rect_cord coord, double &theta, double &phi);
void adjust_ref_point(point tx, point rx, point ref, surface surf);

//bool intersection_ray_xcylinder (xcylinder &xc, ray &r, double &t_value, point &ap);
//bool intersection_ray_cylinder (cylinder &c, ray &r, double &t1, double &t2, point &p1, point &p2);

int check_points_on_plane(finite_plane plane);
// For command line help
void printHelp();
void printSurfaces(finite_plane* surf_list_fp, int num_surfaces);

// Global constants
const int NUM_PHI = 180;
const int NUM_THETA = 91;
const int NUM_VALS = NUM_PHI*NUM_THETA;
const int NUM_FREQ = 11; // There are 11 frequency points - 2GHz, 2.5GHz, 3GHz, ... , 7GHz
ofstream pathinfo;

// HELPER FUNCTIONS _END
/*****************************************************************************************/
// MAIN _BEGIN

int main(int argc, char* argv[])
{

	cout.precision(15);

	// CLOCK STUFF BEGIN
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	cout << "Start Time: " << asctime (timeinfo) << "\n";	
	clock_t start, stop;
	start = 0;
	double t = 0.0;
	// Start timer 
	assert((start = clock())!=-1);
	// CLOCK STUFF END

	// Embed version number 
	int version = 0;
	int subversion = 2;
	cout << "Version " << version << "." << subversion << endl << endl;
	
	char* input_file_name;	
	if (argc < 2)
	{
		// Throw error;
		cout << "Usage:\nraytracer <input_file_name>\n";
		cout << "raytracer <input_file_name> <receiver_locations>\n";
		cout << "raytracer -help\n";
		return 0;
	}
	else
	{
		input_file_name = argv[1];
	}

	if (strcmp(argv[1], "-help") == 0)
	{
		printHelp();
		return(0);
	}

	//char* rec_val_override;	
	/*	int rec_val_override = 0;
	if (argc >= 4)
	{		
	rec_val_override = atoi(argv[3]);
	}
	*/
	// Declare scene parameters
	point trans_loc, rec_loc;
	point_sph trans_angles, rec_angles;
	rect_cord trans_cord, rec_cord;
	double trans_power = 1;
	double resolution = 1;
	double start_freq, end_freq, incr_freq;
	string sweep_param;
	int num_surfaces = get_num_surfaces (input_file_name);	
	int rec_level = 0;
	int image_tree_int = 0;
	int dump_tree_int = 0;
	int dump_paths_int = 0;
	double rel_perm_global = -1;	// Make it some negative number, when it is set it will be positive ...
	double sigma_global = -1;		// ... and this way it can be tested whether these are set or not
	double rho_h_global = -1;
	bool debug = false;	// Used whether to print information on all of the surfaces for debug purposes
	double refcoef_fac = 1;		// Default value is 1 which does not modify the reflection coefficient, if -1 then use the attenuation factor derived from the scattered fields from a finite sized plate
	double patch_red_fac = LARGE_DOUBLE;
	int use_pattern = 1;	// By default receiver pattern is used
	double bw = 15;		// By default set the beamwidth to 15 degrees
	double ref_thresh = 0;	// By default set the threshold to 0
	string tx_fname = "tx_pattern.txt";
	string rx_tc_fname = "rx_pattern_tc.txt";
	string rx_pc_fname = "rx_pattern_pc.txt";
	int rem_red_paths = 0;
	int dump_paths_info = 0;
	double length_interval = 1;
	double angle_interval = 2.5*PI/180;
	int enable_illum_zone_int = 1;	// Controls whether illumination zones are used
	int adj_path_geom = 0;			// By default path is not adjusted 
		
	surface* surf_list = new surface[num_surfaces];
	for (int i=0; i < num_surfaces; i++)
	{
		surf_list[i].rel_perm = -1;	// Make it some negative number, when it is set it will be positive ...	
		surf_list[i].sigma = -1;	// ... and this way it can be tested whether these are set or not
		surf_list[i].rho_h = -1;
	}

	get_params (trans_loc, rec_loc, input_file_name, resolution, 
		trans_power, surf_list, rec_level, start_freq, end_freq, 
		incr_freq, sweep_param, num_surfaces, image_tree_int, 
		rel_perm_global, sigma_global, rho_h_global, dump_tree_int,
		dump_paths_int, refcoef_fac, use_pattern, patch_red_fac, bw,
		ref_thresh, tx_fname, rx_tc_fname, rx_pc_fname, rem_red_paths,
		dump_paths_info, length_interval, angle_interval, enable_illum_zone_int,
		adj_path_geom);

	bool enable_illum_zone = true;
	if (enable_illum_zone_int == 0)
	{
		enable_illum_zone = false;
	}
	

	cout << "tx_fname: " << tx_fname << endl;
	cout << "rx_tc_fname: " << rx_tc_fname << endl;
	cout << "rx_pc_fname: " << rx_pc_fname << endl;

	// If true, image tree is used to store nodes
	// else, a linear list is used
	bool image_tree = false;
	if (image_tree_int == 0)
	{
		image_tree = false;
	}
	else
	{
		image_tree = true;
	}	

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Override any parameters which are specified on the command line
	////////////////////////////////////////////////////////////////////////////////////////////////
	int i = 2;  // Counter to access a certain string of argv
	// Starts at 2 since the .exe string is the 0th string and the env file is the 1th string
	int argumentsValid = 0; // 0 if all of the arguments are correct, -1 otherwise

	char* rec_locs_file_name = NULL;
	char* trans_locs_file_name = NULL;

	// Parse command line arguments
	while(argv[i] != NULL)          // Check for the next string argument
	{
		char* argument = NULL;   // Temporary variable to copy the strings entered by user
		char* variableTarget = NULL;  // Holds the variable to be modified with a new value
		int argIndex = 0;       // Index of the argument that is extracted from the array of strings
		// variable 'i' is the index of the strings

		argument = new char[strlen(argv[i]) + 1];
		variableTarget = new char[strlen(argv[i]) + 1];

		strcpy(argument, argv[i]);

		// No spaces are in these arguments
		// Error check that argument starts with '-'
		if(argument[0] != '-')
		{
			fprintf(stderr, "Arguments should be of the form -variable=NUM \n");
			fprintf(stderr, "Argument '%s' has no '-' \n", argument);
			argumentsValid = -1;
			// put in a break here and exit the program with error
			break;         
		}

		//printf("%s \n", argument);

		// Copy out the variable to be changeed from the argument
		// j starts at 1 to remove the '-'. If 2 dashes are used, j would start at 2
		for(argIndex = 1; argument[argIndex] != '=' && argument[argIndex] != '\0'; argIndex++)     
		{
			variableTarget[argIndex-1] = argument[argIndex];      // [j-1] since there is one '-' and shifting

		}
		variableTarget[argIndex-1] = '\0';	// Place the NULL terminator, argIndex-1 since argument[argIndex] is '='

		// The value argument[argIndex] is now '='. So the character one after the value is what
		// we are interested in. Check that the argIndex i an = otherwise it is an invalid param

		// Error check that argument has an '=' after the variable to be modified
		if (argument[argIndex] != '=' && strcmp(variableTarget, "debug") != 0 )
		{
			fprintf(stderr, "Argument should be of the form -variable=NUM \n");
			fprintf(stderr, "Argument '%s' has no '=' \n", argument);
			argumentsValid = -1;

			delete[] argument;
			delete[] variableTarget;
			argument = NULL;
			variableTarget = NULL;
			break;
		}
		//printf("%s \n", variableTarget);

		// String compare the variableTarget with known variable names and change them accordingly

		if (strcmp(variableTarget, "someIntValue") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];

				// Error check that number has a valid number
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9'))
				{
					fprintf(stderr, "Argument '%s' has an invalid integer \n", argument);
					argumentsValid = -1;                    
					break;
				}

				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber
			int newIntVal = atoi(cmdNumber);   
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%d \n", newIntVal); 
		}
		else if (strcmp(variableTarget, "refs") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];

				// Error check that number has a valid number
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9'))
				{
					fprintf(stderr, "Argument '%s' has an invalid integer \n", argument);
					argumentsValid = -1;                    
					break;
				}

				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber
			rec_level = atoi(cmdNumber);   
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%d \n", rec_level); // Debug to check integer was correctly extracted
		}        
		else if (strcmp(variableTarget, "someFloatValue") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			unsigned int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];                
				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber            
			int dotCount = 0;

			for (j = 0 ; j < strlen(cmdNumber); j++)
			{
				// Error check that number is a valid double
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9') && cmdNumber[j] != '.')
				{
					fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
					argumentsValid = -1;                    
					break;
				}

				if (cmdNumber[j] == '.')
				{
					dotCount++;
					if (dotCount > 1)
					{
						fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
						argumentsValid = -1;                    
						break;
					}
				}
			}
			double newFloatVal = atof(cmdNumber);
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%.10lf \n", newFloatVal); // Debug to check integer was correctly extracted
		}
		else if (strcmp(variableTarget, "grel_perm") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			unsigned int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];

				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber

			int dotCount = 0;
			for (j = 0 ; j < strlen(cmdNumber); j++)
			{
				// Error check that number is a valid double
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9') && cmdNumber[j] != '.')
				{
					fprintf(stderr, "Argument '%s' has an invalid double due to character %c\n", argument, cmdNumber[j]);
					fprintf(stderr, "Extracted '%s' \n", cmdNumber);
					argumentsValid = -1;                    
					break;
				}

				if (cmdNumber[j] == '.')
				{
					dotCount++;
					if (dotCount > 1)
					{
						fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
						argumentsValid = -1;                    
						break;
					}
				}
			}
			rel_perm_global = atof(cmdNumber);  
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%.10lf \n", rel_perm_global); // Debug to check integer was correctly extracted
		}
		else if (strcmp(variableTarget, "gsigma") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			unsigned int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];   
				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber

			int dotCount = 0;
			for (j = 0 ; j < strlen(cmdNumber); j++)
			{
				// Error check that number is a valid double
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9') && cmdNumber[j] != '.')
				{
					fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
					argumentsValid = -1;                    
					break;
				}

				if (cmdNumber[j] == '.')
				{
					dotCount++;
					if (dotCount > 1)
					{
						fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
						argumentsValid = -1;                    
						break;
					}
				}
			}

			sigma_global = atof(cmdNumber);  
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%.10lf \n", sigma_global); // Debug to check integer was correctly extracted
		}
		else if (strcmp(variableTarget, "grho_h") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdNumber = new char[strlen(argv[i] + 1)];
			int k = 0;
			unsigned int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdNumber[j] = argument[k];   
				j++;    // Increment the number string idnex
			}
			cmdNumber[j] = 0;     // append the NULL terminator to the cmdNumber

			int dotCount = 0;
			for (j = 0 ; j < strlen(cmdNumber); j++)
			{
				// Error check that number is a valid double
				if (!(cmdNumber[j] >= '0' && cmdNumber[j] <= '9') && cmdNumber[j] != '.')
				{
					fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
					argumentsValid = -1;                    
					break;
				}

				if (cmdNumber[j] == '.')
				{
					dotCount++;
					if (dotCount > 1)
					{
						fprintf(stderr, "Argument '%s' has an invalid double \n", argument);
						argumentsValid = -1;                    
						break;
					}
				}
			}

			rho_h_global = atof(cmdNumber);  
			delete[] cmdNumber;
			cmdNumber= NULL;
			//printf("%.10lf \n", rho_h_global); // Debug to check integer was correctly extracted
		}
		else if( strcmp(variableTarget, "rec_loc") == 0)
		{
			// Extract the string value from the argument
			char* cmdVal = new char[strlen(argv[i] + 1)];
			int k = 0;
			int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdVal[j] = argument[k];         
				j++;    // Increment the number string idnex
			}
			cmdVal[j] = 0;            // append the NULL terminator 

			// No error checking needed for a string
			rec_locs_file_name = new char[strlen(cmdVal) + 1];
			strcpy (rec_locs_file_name, cmdVal);

			delete[] cmdVal;
			cmdVal = NULL;

			//printf("%s \n", rec_locs_file_name); // Debug to check integer was correctly extracted        
		}
		else if( strcmp(variableTarget, "trans_loc") == 0)
		{
			// Extract the string value from the argument
			char* cmdVal = new char[strlen(argv[i] + 1)];
			int k = 0;
			int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdVal[j] = argument[k];         
				j++;    // Increment the number string idnex
			}
			cmdVal[j] = 0;            // append the NULL terminator 

			// No error checking needed for a string
			trans_locs_file_name = new char[strlen(cmdVal) + 1];
			strcpy (trans_locs_file_name, cmdVal);

			delete[] cmdVal;
			cmdVal = NULL;

			//printf("%s \n", rec_locs_file_name); // Debug to check integer was correctly extracted        
		}
		else if( strcmp(variableTarget, "someStringValue") == 0)
		{
			// Extract the numerical value from the argument
			char* cmdVal = new char[strlen(argv[i] + 1)];
			int k = 0;
			int j = 0;

			// before this for loop starts, argument[argIndex] is '=' so we want argument[argIndex+1]
			for(k = argIndex+1; argument[k] != '\0'; k++)
			{
				cmdVal[j] = argument[k];             
				j++;    // Increment the number string idnex
			}
			cmdVal[j] = 0;            // append the NULL terminator 

			// No error checking needed for a string

			char* newStringValue = new char[strlen(cmdVal) + 1];
			strcpy (newStringValue, cmdVal);

			delete[] cmdVal;
			cmdVal = NULL;
			//printf("%s \n", newStringValue); // Debug to check integer was correctly extracted        
		}
		else if (strcmp(variableTarget, "debug") == 0)
		{
			// Print information on all of the surfaces later
			debug = true;
		}
		else
		{
			// variable that was entered was invalid
			fprintf(stderr, "Argument '%s' has an invalid variable name %s \n", argument, variableTarget);
		}

		// deallocate memory using new
		delete[] argument;
		delete[] variableTarget;
		argument = NULL;
		variableTarget = NULL;
		i++;    // increment i to get the next argument
	}


	if (argumentsValid == -1)
	{
		fprintf(stderr, "Use raytracer.exe -help \n");
		// Exit main due to invalid command line argument program
		return(-1);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	// If rel_perm/sigma/rho_h for a surface is -1, use the global value for it
	// If the global value is not set, throw an error
	for (int i=0; i < num_surfaces; i++)
	{
		if (surf_list[i].rel_perm == -1)
		{
			if (rel_perm_global == -1)
			{
				cout << "Relative permittivity is not specified properly! Specify the global value. \n";
				return 0;
			}
			else
			{
				surf_list[i].rel_perm = rel_perm_global;
			}
		}

		if (surf_list[i].sigma == -1)
		{
			if (sigma_global == -1)
			{
				cout << "Conductivity is not specified properly! Specify the global value.\n";
				return 0;
			}
			else
			{
				surf_list[i].sigma = sigma_global;
			}
		}

		if (surf_list[i].rho_h == -1)
		{
			if (rho_h_global == -1)
			{
				cout << "Surface roughness is not specified properly! Specify the global value.\n";
				return 0;
			}
			else
			{
				surf_list[i].rho_h = rho_h_global;
			}
		}
	}

	// Calculate the bouding cube
	cube outer_cube;
	double xmin = 1e99;
	double xmax = -1e99;
	double ymin = 1e99;
	double ymax = -1e99;
	double zmin = 1e99;
	double zmax = -1e99;

	// Create an array of all the finite planes instead of using surfaces
	// to represent them to avoid unneccessary conversations later
	finite_plane* surf_list_fp = new finite_plane[num_surfaces];
	for (int i = 0; i < num_surfaces; i++)
	{
		surf_list_fp[i] = surf_to_fplane (surf_list[i]);
		// NSOOD: Illumination zone _START
		// NSOOD: Aug 1, 2015 - convert vertices from array of points to point_list		
		point_list* iter;
		for (int j = 0; j < surf_list[i].num_vert; j++)	// cycle through all 4 vertices of the plane
		{
			iter = surf_list[i].v;
			for (int k = 0; k < j; k++)
			{
				iter = iter->next;	// NSOOD: Add error handling
			}

			if (xmin > iter->p.x)
			{
				xmin = iter->p.x;
			}
			if (xmax < iter->p.x)
			{
				xmax = iter->p.x;
			}

			if (ymin > iter->p.y)
			{
				ymin = iter->p.y;
			}
			if (ymax < iter->p.y)
			{
				ymax = iter->p.y;
			}

			if (zmin > iter->p.z)
			{
				zmin = iter->p.z;
			}
			if (zmax < iter->p.z)
			{
				zmax = iter->p.z;
			}			
		}
		// NSOOD: Illumination zone _END
	}

	// Check that the vertices of a surface specify a correct plane
	for (int i = 0; i < num_surfaces; i ++)
	{
		/*
		// ax + by + cz + d = 0
		double a = surf_list_fp[i].unit_normal.x;
		double b = surf_list_fp[i].unit_normal.y;
		double c = surf_list_fp[i].unit_normal.z;
		double d = surf_list_fp[i].d;

		for (int j = 0 ; j < surf_list_fp[i].num_vert; j ++)
		{
		// Check equation of plane on each point
		double value = a * surf_list_fp[i].vl[j].x + 
		b * surf_list_fp[i].vl[j].y + 
		c * surf_list_fp[i].vl[j].z + d;

		if (fabs(value) > 0.00001)
		{
		// A point does not match on the plane
		cout << "The points on surface " << i << " do not specify a plane." << endl;
		return(-1);
		}
		}*/

		/*		int status = check_points_on_plane(surf_list_fp[i]);

		if (status != 0)
		{
		cout << "The points on surface " << i << " do not specify a plane." << endl;
		return(-1);		
		}*/

	}

	if(debug)
	{
		printSurfaces(surf_list_fp, num_surfaces);
	}

	// Add buffer on every side
	double BUFFER = 10; // TODO: this should be ideally proportionally to the largest value
	xmin = xmin - BUFFER;
	ymin = ymin - BUFFER;
	zmin = zmin - BUFFER;
	xmax = xmax + BUFFER;
	ymax = ymax + BUFFER;
	zmax = zmax + BUFFER;

	// Create the 6 planes needed for the bounding cube	
	setup_cube (outer_cube, xmin, xmax, ymin, ymax, zmin, zmax);
	//	cout << "XMIN: " << xmin << "\tYMIN: " << ymin << "\tZMIN: " << zmin << "\n";
	//	cout << "XMAX: " << xmax << "\tYMAX: " << ymax << "\tZMAX: " << zmax << "\n\n";


	// Adjust units on the frequency values from MHz to Hz
	start_freq *= 1e6;
	end_freq *= 1e6;
	incr_freq *= 1e6;

	double freq = start_freq;

	// Parameters for calculation of magnitude of the electric field
	double eta = sqrt(U0/E0);	
	double lamda = SPEED_OF_LIGHT/freq;
	double wave_number = TWOPI/lamda;

	// Override rec value if specified via command line
	/*	if (rec_val_override != 0)
	{
	rec_level = rec_val_override;
	}
	*/	
	cout << "NUM SURFACES: " << num_surfaces << endl;
	cout << "FREQUENCY: " << freq << endl;
	cout << "RECURSION LEVEL: " << rec_level << endl;
	cout << "TRANSMITTER LOCATION: " << trans_loc << endl;
	cout << "RECEIVER LOCATION: " << rec_loc << endl;
	cout << "RADIATED POWER: " << trans_power << endl;

	// Read antenna pattern files into an array
	double theta_gain [NUM_VALS];
	double phi_gain [NUM_VALS];
	read_pattern_files_nff (theta_gain, phi_gain, tx_fname, NUM_VALS);

	const int TCUT_SIZE = 181;
	double theta_gain_rx_tc [TCUT_SIZE];
	double phi_gain_rx_tc [TCUT_SIZE];
	read_pattern_files_nff (theta_gain_rx_tc, phi_gain_rx_tc, rx_tc_fname, TCUT_SIZE);

	const int PCUT_SIZE = 361;
	double theta_gain_rx_pc [PCUT_SIZE];
	double phi_gain_rx_pc [PCUT_SIZE];	
	read_pattern_files_nff (theta_gain_rx_pc, phi_gain_rx_pc, rx_pc_fname, PCUT_SIZE);

	
	if (sweep_param.compare("dist") == 0)
	{
		// Perform a distance sweep
		// In addition to the input file need to read the file with receiver locations
		// If it is not specified, throw an error
		if ((rec_locs_file_name == NULL)||(trans_locs_file_name == NULL))
		{
			// Throw error;
			cout << "Usage: raytracer <input_file_name> -rec_loc=<receiver_locations_file_name> -trans_loc=<transmitter_locations_file_name>\n";
			return 0;
		}		

		// File with results for received power
		ofstream pathgain, fieldinfo;

		char fieldinfobuffer [65536];
		fieldinfo.rdbuf()->pubsetbuf(fieldinfobuffer,65536);

		char pathgainbuffer [65536];
		pathgain.rdbuf()->pubsetbuf(pathgainbuffer,65536);

		pathgain.open ("results.txt");
		if (dump_paths_info == 1)
		{
			pathinfo.open ("path_info.txt");
		}
		fieldinfo.open ("fields.txt");

		ifstream rec_locs,trans_locs;
		trans_locs.open(trans_locs_file_name);
		string s;
		while (!trans_locs.eof())
		{
			getline (trans_locs, s);
			//int check = get_loc_and_angles (s, trans_loc,trans_angles);
			int check = get_ant_info(s, trans_loc, trans_angles, trans_cord);
			if (check == 0)
			{
				//	return 0;
				break;
			}
			// trans_loc is set, perform rest of the analysis

			image_tree_node* root = new image_tree_node;
			root->surf_id = -1;
			root->image_loc = trans_loc;
			root->parent = NULL;
			root->child = NULL;
			root->sibling = NULL;

			//cout << "++++++++++++++ Transmitted at: " << trans_loc << "++++++++++++++ Receiver at: " << rec_loc << endl;
			cout << "Building Image Tree ... " << endl;
			create_image_tree (surf_list, num_surfaces, rec_level, root, outer_cube, enable_illum_zone);
			cout << "Image Tree built ... " << endl;
			//dump_image_tree (root);
			if (dump_tree_int == 1)
			{
				dump_image_tree_rec (root);
			}

			rec_locs.open(rec_locs_file_name);
			while (!rec_locs.eof())
			{
				getline (rec_locs, s);
				// check = get_loc_and_angles (s, rec_loc, rec_angles);
				check = get_ant_info(s, rec_loc, rec_angles, rec_cord);
				if (check == 0)
				{
					//	return 0;
					break;
				}
				// rec_loc is set, perform rest of the analysis			
				cout << endl << "++++++++++++++ Receiver at: " << rec_loc << endl;
				if (dump_paths_int == 1)
				{
				//	pathinfo << "++++++++++++++ Receiver at: " << rec_loc << endl;
				}

				path_list_node* paths;
				int num_paths = 0;
				extract_paths (root, surf_list, num_surfaces, rec_loc, trans_loc, paths, num_paths);			
				compute_path_lengths (paths, trans_loc, rec_loc, trans_cord, surf_list, num_surfaces);

				if (rem_red_paths != 0)
				{
					remove_duplicate_paths(paths, length_interval, angle_interval);	// Remove any duplicate paths that may exist
				}

				if (dump_paths_int == 1)
				{
					dump_paths (paths);
				}
				double rho_h = 0;	// Initialized over here, will be appropriately assigned as per the facet that is being used in the function

				double path_gain_db_from_vect;

				int refs_l = 0;
				//int refs_u = rec_level;

				fieldinfo << trans_loc.x << "\t" << trans_loc.y << "\t" << trans_loc.z << "\t"
						<< rec_loc.x << "\t" << rec_loc.y << "\t" << rec_loc.z << "\t";
				pathgain << trans_loc.x << "\t" << trans_loc.y << "\t" << trans_loc.z << "\t"	
						<< rec_loc.x << "\t" << rec_loc.y << "\t" << rec_loc.z << "\t";

				if (dump_paths_info == 1)
				{
					dump_paths_info_func(paths, trans_loc, rec_loc);
				}

				double ray_traced_power;
				double field_intensity;
				double field_db;
				c_vect field;
				c_vect curr_field;

				for (int refs_u=0; refs_u <= rec_level; refs_u++) 
				{
					if (refs_u == 0)
					{
						field = compute_path_gain_vector_with_rotation (paths, trans_loc, trans_angles, rec_loc, rec_angles, surf_list, 
						num_surfaces, wave_number, rho_h, lamda, freq, trans_power, theta_gain, phi_gain, refs_u, refs_u, refcoef_fac,
						use_pattern, theta_gain_rx_tc, theta_gain_rx_pc, phi_gain_rx_tc, phi_gain_rx_pc, bw, ref_thresh,
						trans_cord, rec_cord, adj_path_geom);
						if (dump_paths_int == 1)
						{
						//	pathinfo << "k = " << refs_u << ",\t" << field    
						//	<< pow(abs(field.x),2) << "\t" << pow(abs(field.y),2) << "\t" << pow(abs(field.z),2) << "\t" << sq_mag(field) << endl;
						}
					}
					else
					{
						curr_field = compute_path_gain_vector_with_rotation (paths, trans_loc, trans_angles, rec_loc, rec_angles, surf_list, 
						num_surfaces, wave_number, rho_h, lamda, freq, trans_power, theta_gain, phi_gain, refs_u, refs_u, refcoef_fac,
						use_pattern, theta_gain_rx_tc, theta_gain_rx_pc, phi_gain_rx_tc, phi_gain_rx_pc, bw, ref_thresh,
						trans_cord, rec_cord, adj_path_geom);
						field = field+curr_field;
						if (dump_paths_int == 1)
						{
						//	pathinfo << "k = " << refs_u << ",\t" << curr_field 
						//	<< pow(abs(curr_field.x),2) << "\t" << pow(abs(curr_field.y),2) << "\t" << pow(abs(curr_field.z),2) << "\t" << sq_mag(curr_field) << endl;
						}
					}

					ray_traced_power = sq_mag(field)*lamda*lamda/(2*eta*4*PI);
					if (ray_traced_power > 0)
					{
						path_gain_db_from_vect = 10*log10(ray_traced_power);
					}
					else
					{
						path_gain_db_from_vect = -500;
					}

					//field_intensity = abs(field.z)*abs(field.z)*lamda*lamda/(2*eta*4*PI);
					field_intensity = sq_mag(field)*lamda*lamda/(2*eta*4*PI);
					if (field_intensity > 0)
					{
						field_db = 10*log10(field_intensity*1e3);
					}
					else
					{
						field_db = -500;
					}

					
					fieldinfo << field.x.real() << "\t" << field.x.imag() << "\t"
						<< field.y.real() << "\t" << field.y.imag() << "\t"
						<< field.z.real() << "\t" << field.z.imag() << "\t";
					
					pathgain << path_gain_db_from_vect << "\t" << field_db << "\t";
				}
				fieldinfo << endl;
				pathgain << endl;
				if (dump_paths_int == 1)
				{
				//	pathinfo << endl;
				}

				// Delete Path List
				if (paths != NULL) my_delete (paths);
			}
			// Delete Image Tree
			rec_locs.close();
			if (root != NULL) my_delete (root);
			cout << "++++++++++++++ Transmitted at: " << trans_loc << endl;
		}
		trans_locs.close();
		pathgain.close();
		if (dump_paths_int == 1)
		{
			pathinfo.close();
		}
		fieldinfo.close();
	}

	// Stop timer
	stop = clock();
	t = (double) (stop-start)/CLOCKS_PER_SEC;
	cout << "RUN TIME: " << t << "\n";

	rawtime;	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	cout << "Stop Time: " << asctime (timeinfo) << "\n";

	return 0;
}
// MAIN _END
/*****************************************************************************************/
// Variable number of nested foor loops using recursion
void loop (int path [], int path_size, int num_surfs, int index,
	finite_plane surf[], point trans_loc, point rec_loc, path_list_node*& list_of_paths)
{
	for (int i = 0; i < num_surfs; i++)
	{
		if ((index == 0) || ((index > 0) && (i != path[index-1])))
		{
			path [index] = i;
			if (index < path_size - 1)
			{
				loop (path, path_size, num_surfs, index+1, surf, trans_loc, rec_loc, list_of_paths);
			}
			else
			{	// The is the last level loop
				gen_paths ( path, path_size, surf, num_surfs, rec_loc, trans_loc,
					list_of_paths
					//					int& path_count
					);
			}
		}
	}
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

vect theta_unit_vect(vect v)
{
	vect theta_v;
	double ax = v.x;
	double ay = v.y;
	double az = v.z;

	double theta = acos(az/sqrt(ax*ax + ay*ay + az*az));
	double phi = atan2(ay, ax);			// Note: This is atan2 NOT atan

	theta_v.x = cos(theta)*cos(phi);
	theta_v.y = cos(theta)*sin(phi);
	theta_v.z = -sin(theta);
	return theta_v;
}

void my_delete (path_node* node)
{
	if (node != NULL)
	{
		if (node->next != NULL)
		{
			my_delete (node->next);
		}
		delete node;
	}
}

void my_delete (image_tree_node* node)
{
	if (node != NULL)
	{
		if (node->child != NULL)
		{
			my_delete (node->child);
		}
		if (node->sibling != NULL)
		{
			my_delete (node->sibling);
		}
		delete node;
	}
}

void my_delete (path_list_node* node)
{
	if (node != NULL)
	{
		if (node->path != NULL)
		{
			my_delete (node->path);
		}
		if (node->next != NULL)
		{
			my_delete (node->next);
		}
		delete node;
	}
}

void my_delete (point_list* node)
{
	if (node != NULL)
	{
		if (node->next != NULL)
		{
			my_delete (node->next);
		}
		delete node;
	}
}

void my_delete (fp_lite_list_node* node)
{
	if (node != NULL)
	{
		if (node->next != NULL)
		{
			my_delete (node->next);
		}
		delete node;
	}
}

/*****************************************************************************************/
// FILE PARSING FUNCTIONS _BEGIN
bool is_line_comment (string s)
{
	int len = s.length();
	bool is_comment = false;
	bool is_not_comment = false;
	string test_char;
	for (int i = 0; (i < len) && (is_comment == false)
		&& (is_not_comment == false); i++)
	{
		test_char = s.substr(i,1);
		// check to see if this line a comment
		if(test_char == "/")
		{		
			is_comment = true;
		}
		// check for all possible types of spaces in the beginning
		else if ((test_char == "\t") ||
			(test_char == "\n") ||
			(test_char == " "))
		{
			// we dont know if it is a comment or not yet
		}
		else
		{
			// it is not a comment, no need to check more characters
			is_not_comment = true;
		}
	}
	return is_comment;
}

void string_to_point (string s, point& p)
{
	// s should be of the form "(x, y, z)"
	string x, y, z;
	int open_par_loc = s.find_first_of('(');
	int first_comma = s.find_first_of(',');
	int second_comma = s.find_last_of(',');
	int close_par_loc = s.find_first_of(')');

	x = s.substr(open_par_loc + 1, first_comma - open_par_loc - 1);
	y = s.substr(first_comma + 1, second_comma - first_comma - 1);
	z = s.substr(second_comma + 1, close_par_loc - second_comma - 1);

	const char* xc = x.c_str();
	const char* yc = y.c_str();
	const char* zc = z.c_str();

	p.x = atof(xc);
	p.y = atof(yc);
	p.z = atof(zc);
}

void string_to_vect (string s, vect& v)
{
	// s should be of the form "(x, y, z)"
	string x, y, z;
	int open_par_loc = s.find_first_of('(');
	int first_comma = s.find_first_of(',');
	int second_comma = s.find_last_of(',');
	int close_par_loc = s.find_first_of(')');

	x = s.substr(open_par_loc + 1, first_comma - open_par_loc - 1);
	y = s.substr(first_comma + 1, second_comma - first_comma - 1);
	z = s.substr(second_comma + 1, close_par_loc - second_comma - 1);

	const char* xc = x.c_str();
	const char* yc = y.c_str();
	const char* zc = z.c_str();

	v.x = atof(xc);
	v.y = atof(yc);
	v.z = atof(zc);
}

void string_to_complex (string s, complex<double>& c)
{
	// s should be of the form "(real, imag)"
	string real, imag;
	int open_par_loc = s.find_first_of('(');
	int first_comma = s.find_first_of(',');
	int close_par_loc = s.find_first_of(')');

	real = s.substr(open_par_loc + 1, first_comma - open_par_loc - 1);
	imag = s.substr(first_comma + 1, close_par_loc - first_comma - 1);

	const char* realc = real.c_str();
	const char* imagc = imag.c_str();

	c.real(atof(realc));
	c.imag(atof(imagc));
}

void string_to_vertices (string s, point_list*& v, int& num_vert)
{	
	string part, left;
	size_t open_par_loc = 0;
	size_t close_par_loc = 0;
	left = s;

	open_par_loc = left.find_first_of('(');
	close_par_loc = left.find_first_of(')');

	num_vert = 0;
	point_list* iter;
	point_list* temp_pl;
	while ((open_par_loc != string::npos) && (close_par_loc != string::npos))
	{				
		part = left.substr(open_par_loc, close_par_loc - open_par_loc);
		left = left.substr(close_par_loc + 1, left.length());
		temp_pl = new point_list;
		point temp_point;
		string_to_point(part, temp_point);
		temp_pl->p = temp_point;
		temp_pl->next = NULL;
		if (num_vert == 0)
		{
			v = temp_pl;
			iter = temp_pl;
		}
		else
		{
			iter->next = temp_pl;
			iter = iter->next;
		}
		num_vert++;
		
		open_par_loc = left.find_first_of('(');
		close_par_loc = left.find_first_of(')');
	}
}

void string_to_surface_type(string param_value, surface_type &st)
{
	if (param_value.compare("PLANE") == 0)	
	{
		st = PLANE;
	}
	else if(param_value.compare("FINITE_PLANE") == 0)
	{
		st = FINITE_PLANE;
	}
	else if(param_value.compare("SPHERE") == 0)
	{
		st = SPHERE;
	}	
}

void get_param_from_line (string &param_name, string &param_value, string s)
{
	int equal_loc = string::npos;     
	equal_loc = s.find_first_of('=');
	param_name = "";
	param_value = "";
	if (equal_loc != string::npos)
	{
		param_name = s.substr(0, equal_loc);
		int semi_co_loc = s.find_first_of(';');
		param_value = s.substr(equal_loc + 1, semi_co_loc - equal_loc - 1);
		char* ws = "\t ";
		param_name.erase(0,param_name.find_first_not_of(ws));
		int tr_sp_loc = param_name.find_first_of(ws);
		if (tr_sp_loc != string::npos)
			param_name.erase(param_name.find_first_of(ws));	

		param_value.erase(0,param_value.find_first_not_of(ws));	
		//cout << param_name << endl << param_value << endl << "\n";
	}
}

void get_params 
	(		
	point &trans_loc,
	point &rec_loc, 
	char* input_file_name,
	double &resolution,
	double &power,
	surface surf[],
	int &rec_level,
	double &start_freq,
	double &end_freq,
	double &incr_freq,
	string &sweep_param,
	int num_surfaces,
	int &image_tree_int,
	double &rel_perm_global,
	double &sigma_global,
	double &rho_h_global,
	int &dump_tree_int, 
	int &dump_paths_int,
	double &refcoef_fac,
	int &use_pattern, 
	double &patch_red_fac,
	double &bw,
	double &ref_thresh,
	string &tx_fname,
	string &rx_tcut_fname,
	string &rx_pcut_fname,
	int &rem_red_paths,
	int &dump_paths_info,
	double &length_interval,
	double &angle_interval,
	int &enable_illum_zone_int,
	int &adj_path_geom
	)
{
	int i = 0;
	bool in_section = false;
	bool crossed_first_surface = false;
	bool increment_pointer = true;

	ifstream input_data;
	input_data.open(input_file_name);
	string s;
	while (!input_data.eof())
	{
		getline(input_data, s);		// read the environment file line by line
		int len = s.length();
		bool is_comment = is_line_comment(s);		
		if (!is_comment && (len > 0))
		{
			// TODO: for now assumed file is well formed. In future add checks
			// check to see what info line has
			string param_name;
			string param_value;
			get_param_from_line (param_name, param_value, s);

			if ((s.substr(0, 1)).compare("{") == 0)
			{
				in_section = true;			
			}
			else if ((s.substr(0, 1)).compare("}") == 0)
			{
				in_section = false;
			}
			else if (((s.substr(0, 1)).compare(" ") == 0) ||
				((s.substr(0, 1)).compare("\t") == 0))
			{
				//	cout << "Ignore this line";
			}

			// Assign value to each parameter
			// TODO: ignore number of surfaces for now - dealing with just 1 plane
			// will be useful later one when there can be an arbitrary amount of surfaces
			if (param_name.compare("transmitter_location") == 0)
			{			
				string_to_point(param_value, trans_loc);
			}
			else if (param_name.compare("receiver_location") == 0)
			{				
				string_to_point(param_value, rec_loc);
			}
			else if (param_name.compare("angular_resolution") == 0)
			{
				resolution = atof(param_value.c_str());
			}
			else if (param_name.compare("radiated_power") == 0)
			{
				power = atof(param_value.c_str());
			}
			else if (param_name.compare("recursion_level") == 0)
			{
				rec_level = atoi(param_value.c_str());
			}
			else if (param_name.compare("start_frequency") == 0)
			{
				start_freq = atof(param_value.c_str());
			}
			else if (param_name.compare("end_frequency") == 0)
			{
				end_freq = atof(param_value.c_str());
			}
			else if (param_name.compare("frequency_incr") == 0)
			{
				incr_freq = atof(param_value.c_str());
			}
			else if (param_name.compare("sweep_param") == 0)
			{
				sweep_param = param_value;
			}
			else if (param_name.compare("image_tree") == 0)
			{
				image_tree_int = atoi(param_value.c_str());
			}
			else if (param_name.compare("dump_tree") == 0)
			{
				dump_tree_int = atoi(param_value.c_str());
			}
			else if (param_name.compare("dump_paths") == 0)
			{
				dump_paths_int = atoi(param_value.c_str());
			}
			else if (param_name.compare("dump_paths_info") == 0)
			{
				dump_paths_info = atoi(param_value.c_str());
			}
			else if (param_name.compare("remove_duplicate_paths") == 0)
			{
				rem_red_paths = atoi(param_value.c_str());
			}
			else if (param_name.compare("enable_illum_zone") == 0)
			{
				enable_illum_zone_int = atoi(param_value.c_str());
			}
			else if (param_name.compare("adj_path_geom") == 0)
			{
				adj_path_geom = atoi(param_value.c_str());
			}
			else if (param_name.compare("g_rel_perm") == 0)
			{
				rel_perm_global = atof(param_value.c_str());
			}
			else if (param_name.compare("g_conductivity") == 0)
			{
				sigma_global = atof(param_value.c_str());
			}
			else if (param_name.compare("g_roughness") == 0)
			{
				rho_h_global = atof(param_value.c_str());
			}
			else if (param_name.compare("refcoef_fac") == 0)
			{
				refcoef_fac = atof(param_value.c_str());
			}
			else if (param_name.compare("use_pattern") == 0)
			{
				use_pattern = atoi(param_value.c_str());
			}
			else if (param_name.compare("patch_red_fac") == 0)
			{
				patch_red_fac = atof(param_value.c_str());
			}
			else if (param_name.compare("beam_width") == 0)
			{
				bw = atof(param_value.c_str());
			}
			else if (param_name.compare("refcoef_thresh") == 0)
			{
				ref_thresh = atof(param_value.c_str());
			}
			else if (param_name.compare("tx_fname") == 0)
			{
				tx_fname = param_value;
			}
			else if (param_name.compare("rx_tcut_fname") == 0)
			{
				rx_tcut_fname = param_value;
			}
			else if (param_name.compare("rx_pcut_fname") == 0)
			{
				rx_pcut_fname = param_value;
			}
			else if (param_name.compare("length_interval") == 0)
			{
				length_interval = atof(param_value.c_str());
			}
			else if (param_name.compare("angle_interval") == 0)
			{
				angle_interval = atof(param_value.c_str());
			}

			if (in_section)
			{
				if (param_name.compare("vertices") == 0)
				{				
					string_to_vertices(param_value, surf[i].v, surf[i].num_vert);
                    // By default a surface's planar field is set to true
                    surf[i].planar = true;
				}				
				else if (param_name.compare("rel_perm") == 0)
				{
					// read relative permittivity					
					surf[i].rel_perm = atof(param_value.c_str());
				}
				else if (param_name.compare("conductivity") == 0)
				{
					// read conductivity					
					surf[i].sigma = atof(param_value.c_str());					
				}
				else if (param_name.compare("roughness") == 0)
				{
					// read roughness					
					surf[i].rho_h = atof(param_value.c_str());					
				}
				else if (param_name.compare("num_vertices") == 0)
				{					
					surf[i].num_vert = atoi(param_value.c_str());
				}
				else if(param_name.compare("outward_normal") == 0)
				{
					// read normal
					string_to_vect(param_value, surf[i].normal);					 
				}
				else if (param_name.compare("constant_d") == 0)
				{
					// read point				
					surf[i].d =  atof(param_value.c_str());
				}				
				else if (param_name.compare("surface_type") == 0)
				{
					string_to_surface_type(param_value, surf[i].s_type);
				}
                else if (param_name.compare("cyl.point") == 0)
                {
                    string_to_point(param_value, surf[i].cyl.p0);
                }
                else if (param_name.compare("cyl.vect") == 0)
                {
                    string_to_vect(param_value, surf[i].cyl.u);
                }
                else if (param_name.compare("cyl.r") == 0)
                {
                    surf[i].cyl.r = atof(param_value.c_str());
                    surf[i].planar = false;
                    // By default a surface's planar field is set to true
                    // It is set to false when the curved surface associated with it is defined
                }

				crossed_first_surface = true;
				increment_pointer = true;
			}
			else
			{
				if (crossed_first_surface == true)
				{
					if(increment_pointer == true)
					{
						surf[i].id = i;
						i = i+1;
						increment_pointer = false;
					}
				}
			}
		}	
	}
	input_data.close();

	// Automatically compute the normal and constant d from the veritices 
	// of the surface. 
	// Also for the Image based raytracer all surfaces are of type FINITE_PLANE
	for (int j = 0; j < num_surfaces; j++)
	{
		// Assign every surface's surface_type to FINITE_PLANE
		surf[j].s_type = FINITE_PLANE;

		// Compute outward normal and constant d for each plane
        point temp;
        temp = surf[j].v->next->p - surf[j].v->p;
		vect n1 = point2vect(temp);	// 2nd pt - 1st pt
        temp = surf[j].v->next->next->p - surf[j].v->next->p;
		vect n2 = point2vect(temp);
		surf[j].normal = normalize(cross(n1, n2));
		surf[j].d = -1*surf[j].normal * point2vect(surf[j].v->p);
		//surf[j].num_vert = 4;	// NSOOD: Hard-coded for now, should make this data-driven

//		cout << "TEST" << endl;
//		cout << surf[j].normal << "\t" << surf[j].d << endl << endl;
//		reduce_surface (surf[j], patch_red_fac);	// when fac is 1, you get the mid-point
		surf[j].shortest_side = shortest_side(surf[j]);

	}
}

int get_num_surfaces (char* input_file_name)
{
	int num_surfaces = 0;
	ifstream input_data;
	input_data.open(input_file_name);
	string s;
	bool found = false;
	while (!input_data.eof() && !found)
	{
		getline(input_data, s);		// read the environment file line by line
		int len = s.length();
		bool is_comment = is_line_comment(s);		
		if (!is_comment && (len > 0))
		{
			string param_name;
			string param_value;
			get_param_from_line (param_name, param_value, s);

			if((s.substr(0, 7)).compare("surface") == 0)
			{
				num_surfaces++;
			}
		}
	}
	input_data.close();
	return num_surfaces;
}

int get_loc (string s, point &rec)
{
	// String s is of the form
	// x_val, y_val, z_val
	int comma_loc = string::npos;
	string rem_string;
	bool fail = false;

	comma_loc = s.find_first_of(',');
	if (comma_loc != string::npos)
	{
		rec.x = atof((s.substr(0,comma_loc)).c_str());
		rem_string = s.substr(comma_loc+1, string::npos);

		comma_loc = string::npos;
		comma_loc = rem_string.find_first_of(',');
		if (comma_loc != string::npos)
		{
			rec.y = atof((rem_string.substr(0,comma_loc)).c_str());
			rec.z = atof((rem_string.substr(comma_loc+1, string::npos)).c_str());
		}
		else
		{
			fail = true;			
		}
	}
	else
	{
		fail = true;
	}

	if (fail)
	{
		cout << "ERROR: File specifying receiver locations is badly formed\n";
		return 0;
	}

	return 1;
}

int get_ant_info (string str, point &loc, point_sph &angle, rect_cord &ant_cord)
{
  // ********************************
  // Convert std string to char*
  // Method 1:
  //char *cstr = &*str.begin(); 

  // Method 2:
  //char *cstr = &str[0u];

  // Method 3:
  // Post by ildjarn @ http://stackoverflow.com/questions/7352099/stdstring-to-char 
  // More info http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char/4152881#4152881  
  vector<char> chars(str.c_str(), str.c_str() + str.size() + 1u);	
  char *cstr = &chars[0];

  int count = 1;
  char *pch = strtok (cstr," ,");
  while (pch != NULL)
  {
  //  printf ("%s\n",pch);

	switch (count)
	{
        case 1: loc.x = atof(pch);
				break;
        case 2: loc.y = atof(pch);
				break;
        case 3: loc.z = atof(pch);
				break;
        case 4: angle.theta = atof(pch);
				break;
        case 5: angle.phi = atof(pch);
                break;              
        case 6: ant_cord.ax.x = atof(pch);
				break;
		case 7: ant_cord.ax.y = atof(pch);
				break;
		case 8: ant_cord.ax.z = atof(pch);
				break;
		case 9: ant_cord.ay.x = atof(pch);
				break;
		case 10: ant_cord.ay.y = atof(pch);
				break;
		case 11: ant_cord.ay.z = atof(pch);
				break;
		case 12: ant_cord.az.x = atof(pch);
				break;
		case 13: ant_cord.az.y = atof(pch);
				break;
		case 14: ant_cord.az.z = atof(pch);
				break;
    }

    pch = strtok (NULL, " ,");
	count++;
  }

  if (count < 6)
  {
    cout << "ERROR: File specifying antenna locations is badly formed\n";
	return 0;
  }

  return 1;
}

int get_loc_and_angles (string s, point &rec, point_sph &angles)
{
	// String s is of the form
	// old format: x_val,y_val,z_val,theta_ant,phi_ant
	// new format: x_val,y_val,z_val,theta_ant,phi_ant,ux,uy,uz,vx,vy,vz,wx,wy,wz
	
	int comma_loc;
	string rem_string;
	bool fail = false;

	comma_loc = s.find_first_of(','); 
	if (comma_loc != string::npos)
	{
		rec.x = atof((s.substr(0,comma_loc)).c_str());
		rem_string = s.substr(comma_loc+1, string::npos);

		comma_loc = rem_string.find_first_of(',');
		if (comma_loc != string::npos)
		{
			rec.y = atof((rem_string.substr(0,comma_loc)).c_str());
			rem_string = rem_string.substr(comma_loc+1, string::npos);
			comma_loc = rem_string.find_first_of(',');
			if (comma_loc != string::npos)
			{
				rec.z = atof((rem_string.substr(0,comma_loc)).c_str());
				rem_string = rem_string.substr(comma_loc+1, string::npos);
				comma_loc = rem_string.find_first_of(',');
				if (comma_loc != string::npos)
				{
					angles.theta = atof((rem_string.substr(0,comma_loc)).c_str());
					rem_string = rem_string.substr(comma_loc+1, string::npos);
					comma_loc = rem_string.find_first_of(',');
					angles.phi = atof(rem_string.c_str());
				}
				else
				{
					fail = true;			
				}
			}
			else
			{
				fail = true;			
			}
		}
		else
		{
			fail = true;			
		}
	}
	else
	{
		fail = true;
	}

	if (fail)
	{
		cout << "ERROR: File specifying receiver locations is badly formed\n";
		return 0;
	}

	return 1;
}


//FILE PARSING FUNCTIONS _END
/*****************************************************************************************/

/*****************************************************************************************/
// Functions to calculate the Fresnel coefficients
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
	double rho = exp(-8*pow(((PI*rho_h*cos(theta_i))/lamda), 2));
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
	double rho = exp(-8*pow(((PI*rho_h*cos(theta_i.real()))/lamda), 2));
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
	double rho = exp(-8*pow(((PI*rho_h*cos(theta_i))/lamda), 2));
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
	double rho = exp(-8*pow(((PI*rho_h*cos(theta_i.real()))/lamda), 2));
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

/*****************************************************************************************/
// Functions to compute INTERSECTIONS
bool intersection_ray_plane (plane p, ray r, ray &ref_ray, double &t_value)
{
	// hit point has the following elements that need to be set in this function
	// point pt - the intersection point
	// ray ref_ray - the reflected ray

	// For a plane: ax + by + cz + d = 0
	// and a ray: r = r_o + t*r_d
	// point of intersection is t = - (r_o*n + d) / r_d*n
	// where n = (a,b,c)
	// The reflected ray is given by
	// r_r = r_i - 2*n (r_i*n)
	double t = 0;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{	
		// point of intersection is also the origin of the reflected ray
		ref_ray.orig = r.orig + t*r.dir;
		ref_ray.dir = r.dir + (-2*r.dir*p.unit_normal)*p.unit_normal;
		t_value = t;
		return true;
	}
	return false;
}

bool intersection_ray_plane_lite (plane_lite p, ray r, double &t_value)
{
	// hit point has the following elements that need to be set in this function
	// point pt - the intersection point
	// ray ref_ray - the reflected ray

	// For a plane: ax + by + cz + d = 0
	// and a ray: r = r_o + t*r_d
	// point of intersection is t = - (r_o*n + d) / r_d*n
	// where n = (a,b,c)
	// The reflected ray is given by
	// r_r = r_i - 2*n (r_i*n)
	double t = 0;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		t_value = t;
		return true;
	}
	return false;
}

bool intersection_ray_inf_finite_plane_lite (finite_plane_lite p, ray r, double &t_value)
{
	// hit point has the following elements that need to be set in this function
	// point pt - the intersection point
	// ray ref_ray - the reflected ray

	// For a plane: ax + by + cz + d = 0
	// and a ray: r = r_o + t*r_d
	// point of intersection is t = - (r_o*n + d) / r_d*n
	// where n = (a,b,c)
	// The reflected ray is given by
	// r_r = r_i - 2*n (r_i*n)
	double t = 0;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		t_value = t;
		return true;
	}
	return false;
}

bool intersection_ray_finite_plane (finite_plane p, ray r, ray &ref_ray, double &t_value)
{
	// ALGORITHM: To find whether the ray interesects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray interesects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	//const double TH = 0.0001;
	const double TH = 1e-7;
	double t = 0;
	point p_int;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		p_int = r.orig + t*r.dir;
		t_value = t;	
		double sum_angles = 0;
		double mag1, mag2;
		point p1, p2;
		bool p_is_vertex = false;

		point* vl = new point[p.num_vert];
		point_list* iter = p.v;
		for (int i=0; i < p.num_vert; i++)
		{
			vl[i] = iter->p;
			iter = iter->next;
		}

		for (int i=0; (i < p.num_vert) && !p_is_vertex; i++)
		{
			p1 = vl[i] - p_int;
			p2 = vl[(i+1)%p.num_vert] - p_int;

			mag1 = sqrt(p1*p1);
			mag2 = sqrt(p2*p2);
			if (mag1*mag2 <= SMALL_DOUBLE)
			{
				// This is true when the intersection point is one of the vertices
				// Consider it inside
				p_is_vertex = true;
			}
			else
			{
				sum_angles += acos((p1*p2)/(mag1*mag2));
			}	
		}

		// Delete vl, don't need it anymore!
		delete[] vl;

		if (p_is_vertex || ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH)))
		{
			ref_ray.orig = p_int;
			ref_ray.dir = r.dir + (-2*r.dir*p.unit_normal)*p.unit_normal;			
			return true;
		}
	}
	
	return false;
}

bool intersection_ray_finite_plane_lite (finite_plane_lite p, ray r, double &t_value)
{
	// ALGORITHM: To find whether the ray interesects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray interesects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	//const double TH = 0.0001;
	const double TH = 1e-7;
	double t = 0;
	point p_int;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		p_int = r.orig + t*r.dir;
		t_value = t;	
		double sum_angles = 0;
		double mag1, mag2;
		point p1, p2;
		bool p_is_vertex = false;

		point* vl = new point[p.num_vert];
		point_list* iter = p.v;
		for (int i=0; i < p.num_vert; i++)
		{
			vl[i] = iter->p;
			iter = iter->next;
		}

		for (int i = 0; (i < p.num_vert) && !p_is_vertex; i ++)
		{
			p1 = vl[i] - p_int;
			p2 = vl[(i+1)%p.num_vert] - p_int;

			mag1 = sqrt(p1*p1);
			mag2 = sqrt(p2*p2);
			//if (mag1*mag2 <= SMALL_DOUBLE)
			if ((mag1 < TH) || (mag2 < TH))
			{
				// This is true when the intersection point is one of the vertices
				// Consider it outside
				//p_is_vertex = true;
				return false;
			}
			else
			{
				sum_angles += acos((p1*p2)/(mag1*mag2));
			}	
		}

		// Delete vl, don't need it anymore!
		delete[] vl;

		//if (p_is_vertex || ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH)))
		if ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH))
		{			
			return true;
		}
	}
	return false;
}

bool intersection_ray_finite_plane_lite (finite_plane p, ray r, double &t_value)
{
	// ALGORITHM: To find whether the ray interesects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray interesects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	//const double TH = 0.0001;
	const double TH = 1e-7;
	double t = 0;
	point p_int;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		p_int = r.orig + t*r.dir;
		t_value = t;	
		double sum_angles = 0;
		double mag1, mag2;
		point p1, p2;
		bool p_is_vertex = false;

		point* vl = new point[p.num_vert];
		point_list* iter = p.v;
		for (int i=0; i < p.num_vert; i++)
		{
			vl[i] = iter->p;
			iter = iter->next;
		}

		for (int i = 0; (i < p.num_vert) && !p_is_vertex; i ++)
		{
			p1 = vl[i] - p_int;
			p2 = vl[(i+1)%p.num_vert] - p_int;

			mag1 = sqrt(p1*p1);
			mag2 = sqrt(p2*p2);
			if (mag1*mag2 <= SMALL_DOUBLE)
			{
				// This is true when the intersection point is one of the vertices
				// Consider it inside
				p_is_vertex = true;
			}
			else
			{
				sum_angles += acos((p1*p2)/(mag1*mag2));
			}	
		}

		// Delete vl, don't need it anymore!
		delete[] vl;

		if (p_is_vertex || ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH)))
		{			
			return true;
		}
	}
	return false;
}

bool intersection_ray_finite_plane_lite (finite_plane_lite3 p, ray r, double &t_value)
{
	// ALGORITHM: To find whether the ray interesects with the finite plane or not
	// compute the sum of the angles between the test point and every pair of edge
	// points. This sum will be 2PI only if the ray interesects with the finite
	// plane. The angle sum will tend to 0 the further away the test point becomes.
	//const double TH = 0.0001;
	const double TH = 1e-7;
	double t = 0;
	point p_int;
	double denom = r.dir*p.unit_normal;
	if (denom != 0)
	{
		t = (-1*(p.unit_normal*r.orig + p.d))/denom;
	}
	if (t > 0)
	{			
		p_int = r.orig + t*r.dir;
		t_value = t;	
		double sum_angles = 0;
		double mag1, mag2;
		point p1, p2;
		bool p_is_vertex = false;

		for (int i = 0; (i < p.num_vert) && !p_is_vertex; i ++)
		{
			p1 = p.vl[i] - p_int;
			p2 = p.vl[(i+1)%p.num_vert] - p_int;

			mag1 = sqrt(p1*p1);
			mag2 = sqrt(p2*p2);
			if (mag1*mag2 <= SMALL_DOUBLE)
			{
				// This is true when the intersection point is one of the vertices
				// Consider it inside
				p_is_vertex = true;
			}
			else
			{
				sum_angles += acos((p1*p2)/(mag1*mag2));
			}	
		}
		if (p_is_vertex || ((sum_angles < TWOPI + TH) && (sum_angles > TWOPI - TH)))
		{			
			return true;
		}
	}
	return false;
}

bool intersection_ray_sphere (sphere s, ray r, ray &ref_ray, double &t_value)
{	
	// solve the quad eqation at^2 + bt + c = 0
	// a = r_d * r_d
	// b = 2 * r_d * (r_o  - s_c)
	// c = (r_o - s_c)*(r_o - s_c) - s_r^2
	// Also need the normal at the point of intersection
	// n = (2(x - xc), 2(y - yc), 2(z - zc))
	vect n;
	double a, b, c, t0, t1;
	t0 = 0;
	t1 = 1;
	a = r.dir*r.dir;
	b = 2*r.dir*(r.orig - s.c);
	c = (r.orig - s.c)*(r.orig - s.c) - s.r*s.r;
	double disc = b*b - 4*a*c;
	if (disc >= 0)
	{
		// Then there is an intersection
		t0 = -(b + sqrt(disc))/(2*a);
		if (t0 >= 0)
		{
			ref_ray.orig = r.orig + t0*r.dir;
			n = sphere_normal(ref_ray.orig, s.c);
			ref_ray.dir = r.dir + (-2*r.dir*n)*n;
			t_value = t0;
			return true;
		}
		else
		{
			t1 = -(b - sqrt(disc))/(2*a);
			if (t1 >= 0)
			{
				ref_ray.orig = r.orig + t1*r.dir;
				n = sphere_normal(ref_ray.orig, s.c);
				ref_ray.dir = r.dir + (-2*r.dir*n)*n;				
				t_value = t1;
				return true;
			}		
		}		
	}
	return false;
}

vect sphere_normal (point p, point center)
{
	vect n;
	point temp = 2*(p - center);
	n.x = temp.x;
	n.y = temp.y;
	n.z = temp.z;
	return n;
}

// Compute the reflection of a given point in a plane
void image_source_loc (finite_plane p, point src, point &img)
{		
	double a = p.unit_normal.x;
	double b = p.unit_normal.y;
	double c = p.unit_normal.z;
	double d = p.d;

	// To find the reflection of the point src in the plane p do:
	// 1) Find parametrized equation of the straight line normal to
	//    the plane and passing through src

	double param_val =  -1*(a*src.x + b*src.y + c*src.z + d)/(a*a + b*b + c*c);

	img.x = src.x + 2*param_val*a;
	img.y = src.y + 2*param_val*b;
	img.z = src.z + 2*param_val*c;
}

void create_image_tree 
	(
	surface surf[], 
	int num_surfs, 
	int num_refls, 
	image_tree_node* root,
	cube c,
	bool enable_illum_zone
	)
{
	// Keep track of how many nodes have been built
	int num_node = 0;
	int num_node_level = 0;

	// Build the first level of nodes under root
	int level = 1;
	cout << " Building Level " << level << " ... " << endl;

	build_root_children (root, surf, num_surfs, num_node);
	cout << "\tBuilt " << num_surfs << endl;

	// Build the rest of the levels corresponding to each higher-order reflection	
	image_tree_node* first_node_iter = root->child; 
	image_tree_node* iter;
	image_tree_node* niter;
	int left = num_refls - 1;	// One level already built	
	bool done = false;
	while (left > 0)
	{
		level++;
		cout << " Building Level " << level << " ... " << endl;
		iter = first_node_iter;
		done = false;
		while (!done && iter != NULL)
		{
			build_node_children (iter, surf, num_surfs, c, num_node, num_node_level, enable_illum_zone);			

			bool found_cousin = false;
			int count = 0;

			if (iter->sibling == NULL)	// All siblings done, check if parent has more siblings
			{
				found_cousin = next_cousin (iter, niter);
				if (found_cousin) { iter = niter; }
				else { done = true; }				
			}
			else
			{
				iter = iter->sibling;
			}
		}		

		left = left - 1;
		if ((left > 0) && first_node_iter != NULL)
		{
			image_tree_node* ncousin;
//			while ((first_node_iter->child == NULL) && (first_node_iter->sibling != NULL))
			while (first_node_iter->child == NULL)
			{
				// if sibling is null ...
				if (first_node_iter->sibling == NULL)
				{
					// ... then find the next cousin
					if(next_cousin(first_node_iter, ncousin))
					{
						first_node_iter = ncousin;
					}
				}
				else // else go to the next sibling
				{
					first_node_iter = first_node_iter->sibling;
				}
			}

			first_node_iter = first_node_iter->child;
		}

		cout << "\tBuilt " << num_node_level << endl; //<< " nodes in the image tree ... \n\n"
		num_node_level = 0;
	}
	cout << "Built " << num_node << " nodes in the image tree ... \n\n";
}

bool next_cousin (image_tree_node* iter, image_tree_node*& niter)
{
	bool found = false;
	bool done = false;
	bool finish_siblings = false;

	int levels = 0;
	int levels_rem = 0;

	image_tree_node* temp_iter;

	while (!found && !done)
	{
		if (iter->parent != NULL)
		{
			iter = iter->parent;
			levels++;

			finish_siblings = false;

			// Check all siblings at this level
			// to find a non-NULL child at the desired level
			while (!finish_siblings && !found)
			{
				if (iter->sibling != NULL)
				{
					iter = iter->sibling;
					levels_rem = levels;

					temp_iter = iter;

					while (levels_rem > 0)
					{
						if (temp_iter->child != NULL)
						{
							temp_iter = temp_iter->child;
							levels_rem--;
						}
						else
						{
							break;
						}
					}

					if (levels_rem == 0)
					{
						// Found right node
						iter = temp_iter;
						found = true;
					}				
				}
				else
				{
					finish_siblings = true;
				}
			}
		}
		else
		{
			done = true;
		}
	}

	if (found) 	{ niter = iter; }

	return found;
}

void build_root_children (image_tree_node* parent, surface surf[], int num_surfs, int& num_node)
{
	image_tree_node* first_child = new image_tree_node;
	first_child->surf_id = surf[0].id;
	point image_loc;
	// need to convert surface into finite plane
	finite_plane p = surf_to_fplane(surf[0]);
	image_source_loc(p , parent->image_loc, image_loc);

	first_child->image_loc = image_loc;
	first_child->parent = parent;
	first_child->sibling = NULL;
	first_child->child = NULL;
	parent->child = first_child;

	num_node++;
	// Print out first_child
	//	cout << first_child << "\n\n";

	image_tree_node* prev_child = first_child;
	image_tree_node* child;

	for (int i = 1; i < num_surfs ; i++)
	{
		p = surf_to_fplane(surf[i]);
		image_source_loc (p, parent->image_loc, image_loc);

		child = new image_tree_node;
		child->surf_id = surf[i].id;
		child->image_loc = image_loc;
		child->parent = parent;
		child->child = NULL;
		child->sibling = NULL;

		prev_child->sibling = child;
		prev_child = child;

		num_node++;
		// Print out other child
		//		cout << child << "\n\n";
	}
}


void build_node_children
	(
	image_tree_node* parent,
	surface surf[],
	int num_surfs,
	cube c,
	int& num_node,
	int& num_node_level,
	bool enable_illum_zone
	)
{
	if (parent != NULL)
	{
		//int nodes_in_level = 0;		// keeps track of how many nodes are in this level

		bool is_first = true;
		bool add_node = true;	// get this from the pyramid plane intersection
//		bool found_base = false;
		image_tree_node* prev_child;
		image_tree_node* child;

		//cout << "Node: " << num_node << endl;

		rect_pyramid zone;
		if (enable_illum_zone)
		{
			// Define Illumination zone
			finite_plane parent_plane = surf_to_fplane(surf[parent->surf_id]);			
			zone.apex = parent->image_loc;
			zone.base = surf_to_fplane_lite(surf[parent->surf_id]);
			// NSOOD: Aug 2015 - Don't need the base obtained using the bounding cube
			// Instead just need the parent_plane to be the base (with modified intersection tests)
			// Since bounding box is not used, zone and ext_zone are the same	
						
			// Create triangular planes representing the faces of the pyramid
			point_list* iter = zone.base.v;
			zone.face_list = new fp_lite_list_node;
			fp_lite_list_node* iterf = zone.face_list;
			for (int i = 0; i < zone.base.num_vert-1; i++)
			{
				init_finite_plane_lite(iterf->fp, zone.apex, iter->p, iter->next->p);
				iterf->next = new fp_lite_list_node;
				iterf = iterf->next;
				iter = iter->next;
			}
			// create last face, iterf is now pointing to the last node of face_list
			init_finite_plane_lite(iterf->fp, zone.apex, iter->p, zone.base.v->p);
			iterf->next = NULL;
		}

		finite_plane p;
		for (int i = 0; i < num_surfs ; i++)
		{
			if (parent->surf_id != surf[i].id)
			{	
				p = surf_to_fplane(surf[i]);

				if(enable_illum_zone) // else it is initialized to true always
				{
					// Print out node info					
					//N cout << "\t\t" << parent->image_loc << endl;
					//N cout << "\t\tSurface ID: " << p.id << endl;
					
					add_node = intersection_finite_plane_pyramid (p, zone);
					if(add_node)
					{
						bool temp = intersection_finite_plane_ext_zone (p, zone);
						if (temp == true)
						{
							add_node = false;
						}
					}					
					
				}					
				//				cout << "\n\n" << add_node << "\n";
				//N cout << "\t\tadd_node: " << add_node << endl << endl;
								
				if (add_node)
				{					
					point image_loc;					
					image_source_loc(p , parent->image_loc, image_loc);


					if (is_first)
					{

						image_tree_node* first_child = new image_tree_node;
						first_child->surf_id = surf[i].id;
						first_child->image_loc = image_loc;
						first_child->parent = parent;
						first_child->sibling = NULL;
						first_child->child = NULL;
						parent->child = first_child;

						prev_child = first_child;											
						is_first = false;

						num_node++;
						//nodes_in_level++;
						num_node_level++;
						// Print out first_child
						//						cout << "First child ... " << parent->surf_id << "\n";
						//						cout << first_child << "\n\n";
					}
					else
					{
						p = surf_to_fplane(surf[i]);
						image_source_loc (p, parent->image_loc, image_loc);

						try
						{
							child = new image_tree_node;
						}
						catch (bad_alloc&)
						{
							cout << "Error: Dynamic Memory Allocation Failed. Program Terminating." << endl;
							cin.get();
						}

						child->surf_id = surf[i].id;
						child->image_loc = image_loc;
						child->parent = parent;
						child->child = NULL;
						child->sibling = NULL;

						prev_child->sibling = child;
						prev_child = child;

						num_node++;
						//nodes_in_level++;
						num_node_level++;
						// Print out other child
						//						cout << child << "\n\n";
					}					
				}
			}
		}

		// Delete illumination zone pyramid - zone
		if (enable_illum_zone)
		{
			my_delete(zone.face_list);
		}

		//		cout << "\n\n";
		//		cout << "\tBuilt " << nodes_in_level << endl;
	}	
}

/*		
//void build_node_children (image_tree_node* parent, surface surf[], int num_surfs, cube c)
void build_node_children (image_tree_node* parent, surface surf[], int num_surfs)
{
if (parent != NULL)
{
if (parent->surf_id != surf[0].id)
{
image_tree_node* first_child = new image_tree_node;
first_child->surf_id = surf[0].id;
point image_loc;
// need to convert surface into finite plane
finite_plane p = surf_to_fplane(surf[0]);
image_source_loc(p , parent->image_loc, image_loc);

first_child->image_loc = image_loc;
first_child->parent = parent;
first_child->sibling = NULL;
first_child->child = NULL;
parent->child = first_child;

image_tree_node* prev_child = first_child;
image_tree_node* child;

for (int i = 1; i < num_surfs ; i++)
{
if (parent->surf_id != surf[i].id)
{
p = surf_to_fplane(surf[i]);
image_source_loc (p, parent->image_loc, image_loc);

child = new image_tree_node;
child->surf_id = surf[i].id;
child->image_loc = image_loc;
child->parent = parent;
child->child = NULL;
child->sibling = NULL;

prev_child->sibling = child;
prev_child = child;
}
}
}
else if (num_surfs >1)
{
image_tree_node* first_child = new image_tree_node;
first_child->surf_id = surf[1].id;
point image_loc;
// need to convert surface into finite plane
finite_plane p = surf_to_fplane(surf[1]);
image_source_loc(p , parent->image_loc, image_loc);

first_child->image_loc = image_loc;
first_child->parent = parent;
first_child->sibling = NULL;
first_child->child = NULL;
parent->child = first_child;

image_tree_node* prev_child = first_child;
image_tree_node* child;

for (int i = 2; i < num_surfs ; i++)
{
p = surf_to_fplane(surf[i]);
image_source_loc (p, parent->image_loc, image_loc);

child = new image_tree_node;
child->surf_id = surf[i].id;
child->image_loc = image_loc;
child->parent = parent;
child->child = NULL;
child->sibling = NULL;

prev_child->sibling = child;
prev_child = child;	
}
}
}
}
*/	

finite_plane surf_to_fplane (surface surf)
{
	finite_plane p;
	p.id = surf.id;
	p.d = surf.d;
	p.rel_perm = surf.rel_perm;
	p.sigma = surf.sigma;
	p.rho_h = surf.rho_h;
	p.unit_normal = surf.normal;
	p.num_vert = surf.num_vert;
	/*for (int j = 0; j < p.num_vert; j++)
	{
		p.vl[j] = surf.v[j];
	}*/
	p.v = surf.v;	// both the plane and surface pointer point to the same point_list structure
	return p;
}

finite_plane_lite surf_to_fplane_lite (surface surf)
{
	finite_plane_lite p;	
	p.d = surf.d;	
	p.unit_normal = surf.normal;
	p.num_vert = surf.num_vert;	
	p.v = surf.v;	// both the plane and surface pointer point to the same point_list structure
	return p;
}

void dump_image_tree (image_tree_node* root)
{
	image_tree_node* first_node_iter = root;
	image_tree_node* iter;
	bool done = false;
	int count = 0;
	int level = 0;

	while (first_node_iter != NULL)
	{
		cout << "+++++++++ Level " << level << " ... " << endl;
		level++;
		iter = first_node_iter;		

		done = false;
		while (!done)
		{
			count++;
			cout << "IMAGE # " << count << endl;
			cout << iter << "\n\n";

			bool found_sibling = false;
			int count = 0;

			if (iter->sibling == NULL)	// All siblings done, check if parent has more siblings
			{
				while ((iter->parent != NULL) && !found_sibling)
				{

					if (iter->parent->sibling != NULL)
					{
						iter = iter->parent->sibling->child;
						found_sibling = true;
					}
					else
					{
						iter = iter->parent;
						count++;
					}
				}

				if (iter->parent == NULL)
				{
					done = true;
				}
				else
				{
					while (count >0)
					{
						iter = iter->child;
						count--;
					}
				}
			}
			else
			{
				iter = iter->sibling;
			}
		}		

		first_node_iter = first_node_iter->child;
	}
	cout << "NUM NODES: " << count << endl;
}

void dump_image_tree_rec (image_tree_node* root)
{
	if (root->child != NULL)
		dump_image_tree_rec(root->child);
	if (root->sibling != NULL)
		dump_image_tree_rec(root->sibling);

	cout << root << "\n\n";
}

void extract_paths_rec 
	(
	image_tree_node* curr_node,
	surface surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_list_node*& path_list_iter1,
	path_list_node*& path_list_iter2,
	int& path_count
	)
{
	if (curr_node->child != NULL)
		extract_paths_rec (curr_node->child, surf, num_surfs, rec_loc, trans_loc, path_list_iter1, path_list_iter2, path_count);
	if (curr_node->sibling != NULL)
		extract_paths_rec (curr_node->sibling, surf, num_surfs, rec_loc, trans_loc, path_list_iter1, path_list_iter2, path_count);


	path_node* path;
	bool path_exists = check_path (curr_node, surf, num_surfs, rec_loc, trans_loc, path);
	if (path_exists)
	{
		path_count++;
		//N		cout << "PATH # " << path_count << endl;

		// Add to linked-list				
		if (path_count > 1)
		{				
			path_list_iter2 = new path_list_node;
			path_list_iter2->path = path;
			path_list_iter2->next = NULL;
			path_list_iter1->next = path_list_iter2;
			path_list_iter1 = path_list_iter1->next;
		}
		else
		{
			path_list_iter2->path = path;
			path_list_iter2->next = NULL;
		}

		// DEBUGGING
		/*		path_node* iter = path_list_iter2->path;
		while (iter != NULL)
		{
		cout << iter->p << endl;
		cout << iter->surf_id << endl << endl;
		iter = iter->next;
		}
		*/		// DEBUGGING
	}
	else
	{
		my_delete (path);
	}
}

void extract_paths 
	(
	image_tree_node* root, 
	surface surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_list_node*& list_of_paths,
	int& path_count
	)
{
	image_tree_node* first_node_iter = root;
	//	image_tree_node* iter;
	bool done = false;
	bool path_exists = false; 
	path_count = 0;
	int level = 0;

	// Create The list of paths
	path_list_node* path_list = new path_list_node; 
	path_list_node* path_list_iter1 = path_list;
	path_list_node* path_list_iter2 = path_list;

	cout << "Extracting Paths ... " << endl;

	extract_paths_rec (first_node_iter, surf, num_surfs, rec_loc, trans_loc, path_list_iter1, path_list_iter2, path_count);

	/*
	while (first_node_iter != NULL)
	{
	cout << "\n\n";
	cout << "+++++++++ Level " << level << " ... " << endl;
	level++;
	iter = first_node_iter;		



	done = false;
	while (!done)
	{			
	path_node* path;
	path_exists = check_path (iter, surf, num_surfs, rec_loc, trans_loc, path);
	if (path_exists)
	{
	path_count++;
	cout << "PATH # " << path_count << endl;

	// Add to linked-list				
	if (path_count > 1)
	{				
	path_list_iter2 = new path_list_node;
	path_list_iter2->path = path;
	path_list_iter2->next = NULL;
	path_list_iter1->next = path_list_iter2;
	path_list_iter1 = path_list_iter1->next;
	}
	else
	{
	path_list_iter2->path = path;
	path_list_iter2->next = NULL;
	}

	// DEBUGGING
	path_node* iter = path_list_iter2->path;
	while (iter != NULL)
	{
	cout << iter->p << endl;
	cout << iter->surf_id << endl << endl;
	iter = iter->next;
	}
	// DEBUGGING
	}
	else
	{
	my_delete (path);
	}

	bool found_sibling = false;
	int count = 0;

	if (iter->sibling == NULL)	// All siblings done, check if parent has more siblings
	{
	while ((iter->parent != NULL) && !found_sibling)
	{

	if (iter->parent->sibling != NULL)
	{
	iter = iter->parent->sibling->child;
	found_sibling = true;
	}
	else
	{
	iter = iter->parent;
	count++;
	}
	}

	if (iter->parent == NULL)
	{
	done = true;
	}
	else
	{
	while (count >0)
	{
	iter = iter->child;
	count--;
	}
	}
	}
	else
	{
	iter = iter->sibling;
	}
	}		

	first_node_iter = first_node_iter->child;
	}
	//*/
	cout << "NUM PATHS: " << path_count << endl << endl;
	if (path_count == 0)
	{
		//	my_delete (path_list);
		path_list = NULL;
	}		
	list_of_paths = path_list;	
}

bool check_path 
	(
	image_tree_node* node, 
	surface surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_node*& path
	)
{
	image_tree_node* iter = node;
	point end_point = rec_loc;

	int node_count = 1;
	path_node* path_iter2 = new path_node;
	path_iter2->next = NULL;
	path_iter2->surf_id = -1;
	path_node* path_iter1 = path_iter2;

	bool pass = true;
	while (iter != NULL)
	{
		ray line;
		line.orig = iter->image_loc;
		point dir = end_point - iter->image_loc;
		line.dir = point2vect (dir);

		double t_val;
		// convert the surface that the ray is supposed to reflect from into a finite_plane
		if (iter->surf_id != -1)
		{
			finite_plane p = surf_to_fplane (surf[iter->surf_id]);
			pass = intersection_ray_finite_plane_lite (p, line, t_val);
			if (t_val > 1)
			{
				pass = false;
			}
			end_point = line.orig + t_val*line.dir;	// point of intersection, this is where the ray reflects

			if (node_count > 1)
			{				
				//node_count++;
				path_iter2 = new path_node;									
				path_iter2->p = end_point;
				path_iter2->ip = iter->image_loc;
				path_iter2->surf_id = iter->surf_id;
				path_iter2->next = NULL;
				path_iter2->next = path_iter1;
				path_iter1 = path_iter2;
			}
			else
			{
				path_iter2->p = end_point;
				path_iter2->ip = iter->image_loc;
				path_iter2->surf_id = iter->surf_id;
				path_iter2->next = NULL;
			}
		}
		if (!pass)
		{
			path = path_iter2;
			return false;
		}
		iter = iter->parent;
		node_count++;
	}

	path = path_iter2;
	// TODO: DO intersections with other planes here NEERAJ

	node_count = 1;
	point start = trans_loc;
	point end;
	bool intersect = false;
	int prev_surf_id;
	ray segment;
	finite_plane p;
	double t_val;

	if (path_iter2->surf_id != -1)
	{
		while (path_iter2 != NULL)
		{			
			end = path_iter2->p;

			segment.orig = start;
			point dir = end - start;
			segment.dir = point2vect (dir);

			if (mag(segment.dir) < 1e-14)
			{
				return false;
			}

			for (int i = 0; i < num_surfs ; i++)
			{
				if (node_count == 1)
				{
					if (path_iter2->surf_id != surf[i].id)
					{
						p = surf_to_fplane(surf[i]);				
						intersect = intersection_ray_finite_plane_lite (p, segment, t_val);
						if (intersect && (t_val < 1))
						{
							return false;
						}					
					}
				}
				else
				{
					if ((path_iter2->surf_id != surf[i].id) &&
						(prev_surf_id != surf[i].id))
					{
						p = surf_to_fplane(surf[i]);					
						intersect = intersection_ray_finite_plane_lite (p, segment, t_val);
						if (intersect && (t_val < 1))
						{
							return false;
						}					
					}
				}				
			}

			start = end;
			node_count++;
			prev_surf_id = path_iter2->surf_id;
			path_iter2 = path_iter2->next;
		}

		end = rec_loc;		
		segment.orig = start;
		point dir = end - start;
		segment.dir = point2vect (dir);

		if (mag(segment.dir) < 1e-14)
		{
			return false;
		}

		for (int i = 0; i < num_surfs ; i++)
		{
			if (prev_surf_id != surf[i].id)
			{
				p = surf_to_fplane(surf[i]);					
				intersect = intersection_ray_finite_plane_lite (p, segment, t_val);
				if (intersect && (t_val < 1))
				{
					return false;
				}					
			}
		}
	}

	else // Check direct path
	{
		end = rec_loc;		
		segment.orig = start;
		point dir = end - start;
		segment.dir = point2vect (dir);
		for (int i = 0; i < num_surfs ; i++)
		{	
			p = surf_to_fplane(surf[i]);			
			intersect = intersection_ray_finite_plane_lite (p, segment, t_val);
			if (intersect && (t_val < 1))
			{
				return false;
			}			
		}
	}

	return true;
}

void gen_paths 
	(
	int path_nodes [],
	int num_path_nodes,
	finite_plane surf[], 
	int num_surfs, 
	point rec_loc,
	point trans_loc,
	path_list_node*& list_of_paths
	//	int& path_count
	)
{
	//	path_list_node* temp_list_node = new path_list_node; // Do this after checkin the path
	path_node* path = new path_node;
	path_node* path_iter = path;
	path_node* path_iter_prev = path;
	path_node* temp_path = NULL;

	point source = trans_loc;
	//	image_source_loc (finite_plane p, point src, point &img);
	for (int i = 0; i < num_path_nodes; i++)
	{
		int surf_id = path_nodes[i];
		image_source_loc (surf[surf_id], source, path_iter->ip);
		path_iter->surf_id = surf_id;

		source = path_iter->ip;

		// Update the prev field beofre path_iter is incremented
		if (i == 0)
		{
			path_iter->prev = NULL;
		}
		else
		{
			path_iter->prev = path_iter_prev;
			path_iter_prev = path_iter_prev->next;
		}

		// Update the next field and allocate memory for next node if needed
		if (i == num_path_nodes-1) // Last node
		{
			path_iter->next = NULL;
		}
		else
		{
			temp_path = new path_node;
			path_iter->next = temp_path;
			path_iter = path_iter->next;
		}		
	}

	// Check to see if the path actually exists	
	path_node* iter = path_iter;	// path_iter points at the end
	point end_point = rec_loc; 
	bool pass = true;

	while (iter != NULL)
	{
		ray line;
		line.orig = iter->ip;
		point dir = end_point - iter->ip;
		line.dir = point2vect (dir);

		double t_val;
		pass = intersection_ray_finite_plane_lite (surf[iter->surf_id], line, t_val);
		if (!pass) { break;	}
		if (t_val > 1)
		{
			pass = false;
			break;
		}
		end_point = line.orig + t_val*line.dir;	// point of intersection, this is where the ray reflects
		iter->p = end_point;

		iter = iter->prev;
	}

	if (pass)
	{
		// Check if any path segments are blocked
		path_iter = path;
		int node_count = 1;
		point start = trans_loc;
		point end;
		bool intersect = false;
		int prev_surf_id;
		ray segment;
		double t_val;

		while (path_iter != NULL)
		{			
			end = path_iter->p;

			segment.orig = start;
			point dir = end - start;
			segment.dir = point2vect (dir);

			if (mag(segment.dir) < 1e-14)
			{
				pass = false;
				break;
			}

			for (int i = 0; i < num_surfs ; i++)
			{
				if (node_count == 1)
				{
					if (path_iter->surf_id != surf[i].id)
					{			
						intersect = intersection_ray_finite_plane_lite (surf[i], segment, t_val);
						if (intersect && (t_val < 1))
						{
							pass = false;
							break;
						}					
					}
				}
				else
				{
					if ((path_iter->surf_id != surf[i].id) &&
						(prev_surf_id != surf[i].id))
					{
						intersect = intersection_ray_finite_plane_lite (surf[i], segment, t_val);
						if (intersect && (t_val < 1))
						{
							pass = false;
							break;
						}					
					}
				}				
			}

			start = end;
			node_count++;
			prev_surf_id = path_iter->surf_id;
			path_iter = path_iter->next;
		} // END WHILE

		if (pass)
		{
			end = rec_loc;		
			segment.orig = start;
			point dir = end - start;
			segment.dir = point2vect (dir);

			if (mag(segment.dir) < 1e-14)
			{
				pass = false;
			}

			if (pass)
			{
				for (int i = 0; i < num_surfs ; i++)
				{
					if (prev_surf_id != surf[i].id)
					{
						intersect = intersection_ray_finite_plane_lite (surf [i], segment, t_val);
						if (intersect && (t_val < 1))
						{
							pass = false;
						}					
					}
				}
			}
		}
	}

	// Add to the list
	if (pass)
	{
		path_list_node* temp_list_node = new path_list_node;
		temp_list_node->path = path;
		list_of_paths->next = temp_list_node;
		// update list_of_paths which is the pointer to the last element in the list
		list_of_paths = list_of_paths->next;
		// TODO: After recursive loop function is finished,
		// set last list_of_paths node's next to NULL

		//iter = path;
		//cout << "+++++++++ Path # " << endl;
		//while (iter != NULL)
		//{
		//	cout << iter->p << endl;
		//	cout << iter->ip << endl;
		//	cout << iter->surf_id << endl << endl;
		//	iter = iter->next;
		//}
		//cout << "******** Path # " << endl;
		//while (path_iter_prev != NULL)
		//{
		//	cout << path_iter_prev->p << endl;
		//	cout << path_iter_prev->ip << endl;
		//	cout << path_iter_prev->surf_id << endl << endl;
		//	path_iter_prev = path_iter_prev->prev;
		//}
	}
	else
	{	// Delete memory
		my_delete (path);
	}
}

bool check_direct_path 
	(
	finite_plane surf[], 
	int num_surfs, 
	point trans_loc, 
	point rec_loc
	)
{
	// This function returns true if the direct path exists and false otherwise
	bool intersect = false;
	double t_val;
	ray segment;

	segment.orig = trans_loc;
	point dir = rec_loc - trans_loc;
	segment.dir = point2vect (dir);
	for (int i = 0; i < num_surfs ; i++)
	{
		intersect = intersection_ray_finite_plane_lite (surf[i], segment, t_val);
		if (intersect && (t_val < 1))
		{
			return false;
		}			
	}

	return true;
}

void dump_paths (path_list_node* path_list)
{
	int path_id = 1;	

	cout << "Dumping all paths ... " << endl;
	while (path_list != NULL)
	{
		path_node* iter = path_list->path;
		cout << "+++++++++ Path # " << path_id << " - " << path_list->path_id << endl;
		while (iter != NULL)
		{
			//cout << iter->p << endl;
			cout << "Reflections: " << path_list->num_refs << endl;
			cout << (iter->p).x << "\t" << (iter->p).y << "\t" << (iter->p).z << endl;
			//		cout << iter->surf_id << endl << endl;
			iter = iter->next;
		}
		cout << "Path Length: " << path_list->path_length << endl << endl;
		path_list = path_list->next;
		path_id++;
	}
}

void dump_paths_info_func (path_list_node* path_list, point trans_loc, point rec_loc)
{

	cout << "Dumping all paths ... " << endl;
	while (path_list != NULL)
	{
		pathinfo << trans_loc.x << "\t" << trans_loc.y << "\t" << trans_loc.z << "\t"
						<< rec_loc.x << "\t" << rec_loc.y << "\t" << rec_loc.z << "\t";
		
		pathinfo << path_list->path_id << "\t";
		pathinfo << path_list->num_refs << "\t";
		pathinfo << path_list->path_length << "\t";
		pathinfo << path_list->adj_path_length << "\t";
		pathinfo << (path_list->tx_theta)*180/PI << "\t";
		pathinfo << (path_list->tx_phi)*180/PI << "\t";		
		
		pathinfo << "\n";

		cout << "Path Length: " << path_list->path_length << endl << endl;
		path_list = path_list->next;

	}
}

void compute_path_lengths (path_list_node* path_list, point trans_loc, point rec_loc, rect_cord tx_cord, surface surf[], int num_surfs)
{	
	double w = 3.9;
	double t_value;
	bool result;

	int gplane_id = num_surfs-1;
	if (num_surfs == 4)
	{
		gplane_id = 0;
	}

	// NSOOD:12
	//bool intersection_ray_xcylinder (xcylinder xc, ray r, double &t_value, point &ap)
	xcylinder xc;	
	xc.yc = 0;
	xc.zc = 1.8;
	xc.r = sqrt(xc.zc*xc.zc + w*w);
	
	ray r;	// The direction and origin of the ray should be assigned based on each reflection point in the loop	
	
	int count = 0;
	char temp [4];
	point last_ref_point, temp_pt;
	vect tx_ray_dir;
	int refs;
	while (path_list != NULL)
	{
		count++;
		string path_id = "T-";
		path_node* iter = path_list->path;
		if (iter->surf_id != -1)
		{
			r.orig = iter->p;
			r.dir = surf[iter->surf_id].normal;

			if (iter->surf_id != gplane_id)
				result = intersection_ray_xcylinder (xc, r, t_value, iter->ap);
			else
				iter->ap = iter->p;
			

			double adj_path_length = dist(trans_loc, iter->ap);

            temp_pt = iter->p - trans_loc;
			tx_ray_dir = point2vect(temp_pt);

			double path_length = dist(trans_loc, iter->p);
			refs = 1;
			while (iter != NULL)
			{
				if (iter->next == NULL)
				{
					path_length += dist(iter->p, rec_loc);
					adj_path_length += dist(iter->ap, rec_loc);
					//P					cout << dist(iter->ip, rec_loc) << endl;
					last_ref_point = iter->p;
				}
				else
				{
					r.orig = iter->next->p;
					r.dir = surf[iter->next->surf_id].normal;
					
					if (iter->next->surf_id != gplane_id)
						result = intersection_ray_xcylinder (xc, r, t_value, iter->next->ap);
					else
						iter->next->ap = iter->next->p;

					adj_path_length += dist(iter->ap, iter->next->ap);

					path_length += dist(iter->p, iter->next->p);
					refs++; // Updating the reflection counter
				}

				//itoa (iter->surf_id, temp, 10);
                sprintf(temp, "%d", iter->surf_id);
                path_id = path_id + temp + "-";

				iter = iter->next;
			}
			path_list->path_length = path_length;
			path_list->adj_path_length = adj_path_length;
			path_list->num_refs = refs;
		}
		else	// direct path
		{
            temp_pt = rec_loc - trans_loc;
			tx_ray_dir = point2vect(temp_pt);
			path_list->path_length = dist(trans_loc, rec_loc);
			path_list->adj_path_length = path_list->path_length;
			last_ref_point = trans_loc;
			path_list->num_refs = 0;	// Direct path undergoes 0 reflections
		}

		point rx_antenna_ref = rec_loc - last_ref_point;
		point_sph rx_antenna_ref_sph;
		get_sph_ref (rx_antenna_ref, rx_antenna_ref_sph);

		path_id += "R";
		//		cout << path_id << "\t" << path_list->path_length << "\n";
		//		cout << path_id << "\t" << path_list->path_length << "\t" << rx_antenna_ref_sph.theta * 180/PI
		//			<< "\t" << rx_antenna_ref_sph.phi * 180/PI << endl;
		path_list->path_id = path_id;

		// Calculate departure angles
		compute_departure_angles(tx_ray_dir, tx_cord, path_list->tx_theta, path_list->tx_phi);
	
		// Move to next path
		path_list = path_list->next;
	}
	//	cout << endl << endl;
	//cout << "NUM PATHS: " << count << endl;
}

c_vect compute_path_gain_vector_with_rotation
(
	path_list_node* path_list, 
	point trans_loc,
	point_sph trans_angles,
	point rec_loc,
	point_sph rec_angles,
	surface surf[], 
	int num_surfs,
	double k,
	double rho_h,
	double lamda,
	double freq,
	double trans_power,
	double theta_gain[],
	double phi_gain[],
	int refs_l,
	int refs_u,
	double refcoef_fac,
	int &use_pattern, 
	double theta_gain_rx_tc[], 
	double theta_gain_rx_pc[],
	double phi_gain_rx_tc[],
	double phi_gain_rx_pc[],
	double bw,
	double refcoef_thresh,
	rect_cord trans_cord,
	rect_cord rec_cord,
	int adj_path_geom
)
{	
	cdouble propagation_term (0, 0);
	c_vect field = init ();
	c_vect field_local = init ();

	int num_paths = 0;
	double refcoef_fac_te, refcoef_fac_tm;
	double theta_tx, phi_tx, theta_rx, phi_rx;
	double rd_ax_tx, rd_ay_tx, rd_az_tx; // x, y and z component of the intial ray direction in tx co-ordinates
	double rd_ax_rx, rd_ay_rx, rd_az_rx; // x, y and z component of the intial ray direction in rx co-ordinates
    point temp_pt;
    
	while (path_list != NULL)
	{
		if (path_list->num_refs >= refs_l && path_list->num_refs <= refs_u)
		{
			path_node* iter = path_list->path;
			path_node* last = NULL;

			double r;
			if (adj_path_geom == 1)
			{
				r = path_list->adj_path_length;				
			}
			else
			{
				r = path_list->path_length;
			}
			
			propagation_term.real(cos(-1*k*r)/r);
			propagation_term.imag(sin(-1*k*r)/r);

			point last_ref_point;
			init_point(last_ref_point, LARGE_DOUBLE, LARGE_DOUBLE, LARGE_DOUBLE);
			if (iter->surf_id != -1) // All paths except the direct path
			{
				// Multiply this with all the reflection coefficients
				point first, second;
				first = trans_loc;
				second = iter->p;
				cdouble te_refc, tm_refc, te_tranc, tm_tranc;

				// Get direction of the initial ray
                temp_pt = second - first;
				vect ray_dir = point2vect (temp_pt);

				// Temp reference points
				point temp_rect_ref;
				point_sph temp_sph_ref;

				// Define the reference point at which the spherical co-ordinate
				// system is established, rect_ref, sph_ref
				point rect_ref = second - first;
				point_sph sph_ref;
				get_sph_ref (rect_ref, sph_ref);

				// The vertical angle theta is between the ray_dir vector and the az component
				// of the transmitter co-ordinate system
				// acos produces a result that is between 0 and pi
				rd_az_tx = ray_dir*trans_cord.az;
				theta_tx = acos(rd_az_tx/(mag(trans_cord.az)*mag(ray_dir)));
				rd_ax_tx = ray_dir*trans_cord.ax;
				rd_ay_tx = ray_dir*trans_cord.ay;
				// atan2 produces a result in the range [-pi,pi]
				// when result is negative add 2pi, hence result will be in the range [0,2pi]
				phi_tx = atan2(rd_ay_tx, rd_ax_tx);
				if (phi_tx < 0)
				{
					phi_tx = phi_tx + TWOPI;
				}

				// Get the initial Electric Field from the antenna pattern function (Tx Pattern Applied)
				//c_vect2 init_field_sph = efield_pattern (keep_between_zero_and_2pi(sph_ref.theta+trans_angles.theta), keep_between_zero_and_2pi(sph_ref.phi+trans_angles.phi), trans_power, theta_gain, phi_gain);
				c_vect2 init_field_sph = efield_pattern (theta_tx, phi_tx, trans_power, theta_gain, phi_gain);

				// Apply Rx Antenna Pattern
				// Need to find last reflection point
				
				path_node* temp_iter = iter;
				while (temp_iter->next != NULL)
				{
					temp_iter = temp_iter->next;
				}
				last_ref_point = temp_iter->p;

                temp_pt = last_ref_point - rec_loc;
				vect rx_ray_dir = point2vect(temp_pt);
				
				if (use_pattern == 1)
				{
					// The vertical angle theta is between the ray_dir vector and the az component
					// of the transmitter co-ordinate system
					// acos produces a result that is between 0 and pi
					rd_az_rx = rx_ray_dir*rec_cord.az;
					theta_rx = acos(rd_az_rx/(mag(rec_cord.az)*mag(rx_ray_dir)));
					rd_ax_rx = rx_ray_dir*rec_cord.ax;
					rd_ay_rx = rx_ray_dir*rec_cord.ay;
					// atan2 produces a result in the range [-pi,pi]
					// when result is negative add 2pi, hence result will be in the range [0,2pi]
					phi_rx = atan2(rd_ay_rx, rd_ax_rx);
					if (phi_rx < 0)
					{
						phi_rx = phi_rx + TWOPI;
					}
					
					pathinfo << trans_loc.x << "\t" << trans_loc.y << "\t" << trans_loc.z << "\t"
						<< rec_loc.x << "\t" << rec_loc.y << "\t" << rec_loc.z << "\t";
		
					pathinfo << path_list->path_id << "\t";
					pathinfo << path_list->num_refs << "\t";
					pathinfo << path_list->path_length << "\t";
					pathinfo << path_list->adj_path_length << "\t";
					pathinfo << theta_rx*180/PI << "\t";
					pathinfo << phi_rx*180/PI << "\t";					
					//pathinfo << (path_list->rx_phi)*180/PI << "\t";
					//pathinfo << (path_list->rx_theta)*180/PI << "\t";
					
					pathinfo << "\n";

					init_field_sph.theta *= gain_theta_rx(theta_rx, phi_rx, theta_gain_rx_tc, theta_gain_rx_pc);
					init_field_sph.phi *= gain_theta_rx(theta_rx, phi_rx, phi_gain_rx_tc, phi_gain_rx_pc);

					/*
					init_field_sph.theta *= gain_theta_rx (keep_between_zero_and_2pi(rx_antenna_ref_sph.theta+rec_angles.theta), 
														   keep_between_zero_and_2pi(rx_antenna_ref_sph.phi+rec_angles.phi), 
														   theta_gain_rx_tc, theta_gain_rx_pc);
					init_field_sph.phi *= gain_theta_rx (keep_between_zero_and_2pi(rx_antenna_ref_sph.theta+rec_angles.theta), 
													   keep_between_zero_and_2pi(rx_antenna_ref_sph.phi+rec_angles.phi), 
													   phi_gain_rx_tc, phi_gain_rx_pc);
					*/
				}
				else // isotropic
				{
					init_field_sph.theta *= 1;
				}			

				c_vect2 field_sph_cur, field_sph_next;
				c_vect field_rect_cur = init ();
				c_vect field_rect_next = init ();
				field_sph_cur = init_field_sph;

				rect_cord global; 
				get_global_cord (global);

				rect_cord old_sys = global;
				rect_cord new_sys;

			//	if (dump_path_info == 1) cout << "----------------------" << endl;

				while (iter != NULL)
				{
					/**************************************************/
					second = iter->p;
					
		/*			if (dump_path_info == 1)
					{
						cout << first << " -> " << second << endl;
						cout << iter->surf_id << endl;
					}
		*/
					// Find angle of incidence
					// arccos (A.B/|A||B|)
                    temp_pt = first - second;
					ray_dir = point2vect (temp_pt);

					vect n = surf[iter->surf_id].normal;
					
//					if (dump_path_info == 1) cout << "Plane normal: " << n << endl;

					double theta_i = acos((n*ray_dir)/(mag(n)*mag(ray_dir)));

					if (theta_i > PI/2)
					{
						theta_i = PI - theta_i;
						n = -1*n;
						surf[iter->surf_id].normal = n;
						surf[iter->surf_id].d = -surf[iter->surf_id].d;
					}

					//cout << angle << endl;
//					if (dump_path_info == 1) cout << F180OVERPI*theta_i << endl;

					cdouble e_t;
					e_t.real(surf[iter->surf_id].rel_perm);
					e_t.imag(surf[iter->surf_id].sigma/(E0*2*PI*freq));
					cdouble e_i(1,0);	// Air
					cdouble theta_t = trans_angle (theta_i, e_i, e_t);
					rho_h = surf[iter->surf_id].rho_h;

					// NSOOD: 2016-01-25: Compute attenuation factor applied due to finiteness of the facet
					// START
					// double a = shortest_side(surf[iter->surf_id]);
					//double bw = 15;					
					if (refcoef_fac == -1) // When refcoef_fac is set to -1, use attenuation factor derived from scattered field
					{
						refcoef_fac_te = compute_att_fac(surf[iter->surf_id].shortest_side, theta_i, k, bw, true);
						refcoef_fac_tm = compute_att_fac(surf[iter->surf_id].shortest_side, theta_i, k, bw, false);
					}
					else
					{
						refcoef_fac_te = refcoef_fac;
						refcoef_fac_tm = refcoef_fac;
					}
					// END

					tm_coeff (theta_i, theta_t, e_i, e_t, tm_refc, tm_tranc, rho_h, lamda, refcoef_fac_tm, refcoef_thresh);
					te_coeff (theta_i, theta_t, e_i, e_t, te_refc, te_tranc, rho_h, lamda, refcoef_fac_te, refcoef_thresh);

					//PG_term_local *= refc;
					/**************************************************/			

					// Find the reference point in global co-ordinates
					temp_rect_ref = second - first;
					get_sph_ref (temp_rect_ref, temp_sph_ref);
					// Then convert to facet-fixed
					rect_ref = rect_to_rect (temp_rect_ref, global, old_sys);
					get_sph_ref (rect_ref, sph_ref);

					// E-field represented by init_field_rect is in the co-ordinate frame
					// local to the previous facet or antenna axis for the first pass
					// To convert this to the current facet-fixed co-ordinates				
					// Step: 1 Convert to current cartesian co-ordinates
					field_rect_cur = sph_to_rect (field_sph_cur, sph_ref.theta, sph_ref.phi);

					// Step: 2a Convert to next facet-fixed cartesian co-ordinates				
					new_sys = get_local_cord (surf[iter->surf_id]);
					field_rect_next = rect_to_rect (field_rect_cur, old_sys, new_sys);

					// Step: 2b Convert reference point from current to next co-ordinate frame			
					rect_ref = rect_to_rect (rect_ref, old_sys, new_sys);
					get_sph_ref (rect_ref, sph_ref);

					// Step: 3 Convert to next spherical co-ordinates at point sph_ref
					field_sph_next = rect_to_sph (field_rect_next, sph_ref.theta, sph_ref.phi);

					// At this point the theta component of the field gets the TM coeff
					// and the phi component gets the TE coeff
					field_sph_next.phi *= te_refc;
					field_sph_next.theta *= tm_refc;

					// Update co-ordinate systems, field_sph_cur
					old_sys = new_sys;
					field_sph_cur = field_sph_next;

					// Update first, iter
					first = second;
					if (iter->next == NULL)
					{
						last = iter;
					}
					iter = iter->next;
				}

//				if (dump_path_info == 1) cout << "----------------------" << endl << endl << endl;

				second = rec_loc;
				// Find the reference point in global co-ordinates
				temp_rect_ref = second - first;
				get_sph_ref (temp_rect_ref, temp_sph_ref);
				// Then convert to facet-fixed
				rect_ref = rect_to_rect (temp_rect_ref, global, old_sys);
				get_sph_ref (rect_ref, sph_ref);

				field_rect_cur = sph_to_rect (field_sph_cur, sph_ref.theta, sph_ref.phi);
				field_rect_next = rect_to_rect (field_rect_cur, old_sys, global);
				field_local = propagation_term * field_rect_next;
			}
			else  // Direct path
			{
				point ray_dir = rec_loc - trans_loc;
				point_sph sph_ref;
				get_sph_ref (ray_dir, sph_ref);

				c_vect2 init_field_sph = efield_pattern (keep_between_zero_and_2pi(sph_ref.theta+trans_angles.theta), keep_between_zero_and_2pi(sph_ref.phi+trans_angles.phi), trans_power, theta_gain, phi_gain);
				//// Apply Rx Antenna Pattern
				point rx_antenna_ref = trans_loc - rec_loc;
				point_sph rx_antenna_ref_sph;
				get_sph_ref (rx_antenna_ref, rx_antenna_ref_sph);
				//init_field_sph.phi *= gain (phi_gain, rx_antenna_ref_sph.theta, rx_antenna_ref_sph.phi, freq_index);
				//init_field_sph.theta *= gain (theta_gain, rx_antenna_ref_sph.theta, rx_antenna_ref_sph.phi, freq_index);
				//			init_field_sph.phi *= gain (phi_gain, rx_antenna_ref_sph.theta, rx_antenna_ref_sph.phi);
				//			init_field_sph.theta *= gain (theta_gain, rx_antenna_ref_sph.theta, rx_antenna_ref_sph.phi);
				if (use_pattern == 1)
				{
					init_field_sph.theta *= gain_theta_rx (keep_between_zero_and_2pi(rx_antenna_ref_sph.theta+rec_angles.theta),
														   keep_between_zero_and_2pi(rx_antenna_ref_sph.phi+rec_angles.phi), 
														   theta_gain_rx_tc, theta_gain_rx_pc);
					init_field_sph.phi *= gain_theta_rx (keep_between_zero_and_2pi(rx_antenna_ref_sph.theta+rec_angles.theta),
													   keep_between_zero_and_2pi(rx_antenna_ref_sph.phi+rec_angles.phi), 
													   phi_gain_rx_tc, phi_gain_rx_pc);
				}
				else // isotropic
				{
					init_field_sph.theta *= 1;
				}

				field_local = sph_to_rect (init_field_sph, sph_ref.theta, sph_ref.phi);
				field_local = propagation_term * field_local;
			}

			// pathinfo << "\t" << path_list->path_id << "\t" << path_list->num_refs << "\t" 
			//	<< path_list->path_length << "\t" << field_local << endl;
			// for Xingqi
		//	pathinfo << last_ref_point.x << "\t" << last_ref_point.y << "\t" << last_ref_point.z << "\t" 
		//		     << rec_loc.x << "\t" << rec_loc.y << "\t" << rec_loc.z << "\t" << field_local << endl;
			num_paths++;

			field = field + field_local;
		}
		path_list = path_list->next;		
	}

//	pathinfo << "NUM PATHS: " << num_paths << endl;

	return field;
}

c_vect2 efield_pattern 
	(
	double theta, 
	double phi, 
	double trans_power,
	double theta_gain[], 
	double phi_gain[]
)
{
	c_vect2 efield = init2();
	//	cout << "/***************************************/\n";
	double eta = sqrt (U0/E0);
	double Pt = trans_power;				// Transmitted Power
	double field_coeff = sqrt(eta * Pt/TWOPI);
	//efield.theta = field_coeff*sin(theta);		// TODO: Read this from a file
	//	efield.phi = field_coeff*0;			// TODO: Read this from a file

	// HACK!!!!!
	// modulate rays with elevation cut pattern based on their theta value.
	//	efield.theta = field_coeff*gain_theta(theta, phi);
	//	efield.phi = field_coeff*0;



	//	cout << "efield: " << efield << endl;

	efield.theta = field_coeff*gain(theta_gain, theta, phi);
	efield.phi = field_coeff*gain(phi_gain, theta, phi);	
	//	cout << "/***************************************/\n";
	return efield;

}

rect_cord get_local_cord (surface surf)
{
	rect_cord local;
    point temp_pt;
    temp_pt = surf.v->next->p - surf.v->p;
	local.az = normalize (surf.normal);
	local.ax = normalize (point2vect (temp_pt));
	local.ay = cross (local.az, local.ax);
	return local;
}

c_vect sph_to_rect (c_vect2 v_sph, double theta, double phi)
{
	cdouble	Er (0, 0);
	cdouble Et = v_sph.theta;
	cdouble Ep = v_sph.phi;

	c_vect v_rect;
	v_rect.x = sin(theta)*cos(phi)*Er + cos(theta)*cos(phi)*Et - sin(phi)*Ep;
	v_rect.y = sin(theta)*sin(phi)*Er + cos(theta)*sin(phi)*Et + cos(phi)*Ep;
	v_rect.z = cos(theta)*Er - sin(theta)*Et;
	return v_rect;
}

c_vect2 rect_to_sph (c_vect v_rect, double theta, double phi)
{		
	c_vect2 v_sph;
	cdouble Ex = v_rect.x;
	cdouble Ey = v_rect.y;
	cdouble Ez = v_rect.z;

	v_sph.phi = cos(phi)*Ey - sin(phi)*Ex;
	v_sph.theta = cos(theta)*cos(phi)*Ex + cos(theta)*sin(phi)*Ey - sin(theta)*Ez;

	return v_sph;
}

c_vect rect_to_rect (c_vect v_rect_old, rect_cord old_sys, rect_cord new_sys)
{	
	cdouble Ex = v_rect_old.x;
	cdouble Ey = v_rect_old.y;
	cdouble Ez = v_rect_old.z;
	vect ax = old_sys.ax;
	vect ay = old_sys.ay;
	vect az = old_sys.az;
	vect axp = new_sys.ax;
	vect ayp = new_sys.ay;
	vect azp = new_sys.az;

	c_vect v_rect_new;
	v_rect_new.x = (axp*ax)*Ex + (axp*ay)*Ey + (axp*az)*Ez;
	v_rect_new.y = (ayp*ax)*Ex + (ayp*ay)*Ey + (ayp*az)*Ez;
	v_rect_new.z = (azp*ax)*Ex + (azp*ay)*Ey + (azp*az)*Ez;

	return v_rect_new;
}

vect rect_to_rect (vect v_rect_old, rect_cord old_sys, rect_cord new_sys)
{
	double Ax = v_rect_old.x;
	double Ay = v_rect_old.y;
	double Az = v_rect_old.z;
	vect ax = old_sys.ax;
	vect ay = old_sys.ay;
	vect az = old_sys.az;
	vect axp = new_sys.ax;
	vect ayp = new_sys.ay;
	vect azp = new_sys.az;

	vect v_rect_new;
	v_rect_new.x = (axp*ax)*Ax + (axp*ay)*Ay + (axp*az)*Az;
	v_rect_new.y = (ayp*ax)*Ax + (ayp*ay)*Ay + (ayp*az)*Az;
	v_rect_new.z = (azp*ax)*Ax + (azp*ay)*Ay + (azp*az)*Az;

	return v_rect_new;
}

point rect_to_rect (point v_rect_old, rect_cord old_sys, rect_cord new_sys)
{
	double Ax = v_rect_old.x;
	double Ay = v_rect_old.y;
	double Az = v_rect_old.z;
	vect ax = old_sys.ax;
	vect ay = old_sys.ay;
	vect az = old_sys.az;
	vect axp = new_sys.ax;
	vect ayp = new_sys.ay;
	vect azp = new_sys.az;

	point v_rect_new;
	v_rect_new.x = (axp*ax)*Ax + (axp*ay)*Ay + (axp*az)*Az;
	v_rect_new.y = (ayp*ax)*Ax + (ayp*ay)*Ay + (ayp*az)*Az;
	v_rect_new.z = (azp*ax)*Ax + (azp*ay)*Ay + (azp*az)*Az;

	return v_rect_new;
}

void get_global_cord (rect_cord &global)
{
	global.ax.x = 1;
	global.ax.y = 0;
	global.ax.z = 0;
	global.ay.x = 0;
	global.ay.y = 1;
	global.ay.z = 0;
	global.az.x = 0;
	global.az.y = 0;
	global.az.z = 1;
}

void get_sph_ref (point rect_ref, point_sph &sph_ref)
{
	double Ax = rect_ref.x;
	double Ay = rect_ref.y;
	double Az = rect_ref.z;

	sph_ref.theta = acos(Az/sqrt(Ax*Ax + Ay*Ay + Az*Az));
	sph_ref.phi = atan2(Ay, Ax);			// Note: This is atan2 NOT atan
	sph_ref.phi = keep_between_zero_and_2pi(sph_ref.phi);			// This makes sure that the value is in [0:2*PI)
	sph_ref.r = sqrt(rect_ref*rect_ref);	
}

double keep_between_zero_and_2pi(double angle)
{
	return fmod(angle,2*PI) + (fmod(angle,2*PI) >= 0 ? 0 : 2*PI);			// This makes sure that the value is in [0:2*PI)
}

double gain 
	(
	double antenna_pattern[], 
	double theta, 
	double phi
	//	int freq_index
	)
{
	double theta_d = theta * 180/PI;
	double phi_d = phi * 180/PI;

	double incr = 2;
	int theta_p = (int)(theta_d/incr);
	int phi_p = (int)(phi_d/incr);

	int index1 = theta_p * NUM_PHI + phi_p;
	int index2 = index1 +  1;
	int index3 = (theta_p + 1) * NUM_PHI + phi_p;
	int index4 = index3 + 1;

	// The distances of the point from the 4 vertices are the weights
	double w1 = delta (theta_d*incr, phi_d*incr, theta_p, phi_p);
	double w2 = delta (theta_d*incr, phi_d*incr, theta_p, phi_p + incr);
	double w3 = delta (theta_d*incr, phi_d*incr, theta_p + incr, phi_p);
	double w4 = delta (theta_d*incr, phi_d*incr, theta_p + incr, phi_p + incr);

/*	double g1 = pow(10, (antenna_pattern[index1])/20);
	double g2 = pow(10, (antenna_pattern[index2])/20);
	double g3 = pow(10, (antenna_pattern[index3])/20);
	double g4 = pow(10, (antenna_pattern[index4])/20);
	return (w1*g1 + w2*g2 + w3*g3 + w4*g4)/(w1 + w2 + w3 + w4);
*/
	
	double g1 = antenna_pattern[index1];
	double g2 = antenna_pattern[index2];
	double g3 = antenna_pattern[index3];
	double g4 = antenna_pattern[index4];
	
	return pow(10,(w1*g1 + w2*g2 + w3*g3 + w4*g4)/(w1 + w2 + w3 + w4)*(1.0/20.0));

	
}

double gain_theta
	(	
	double theta, 
	double phi	
	)
{
	double theta_d = theta * 180/PI;
	double phi_d = phi * 180/PI;

	double theta_cut [181] = { 
		-13.93860334	,
		-13.930978	,
		-13.61302835	,
		-13.38261092	,
		-13.13774462	,
		-12.91396461	,
		-12.65623829	,
		-12.44851413	,
		-12.0403135	,
		-11.63751742	,
		-11.44309703	,
		-11.22813491	,
		-10.90963247	,
		-10.68522962	,
		-10.58305684	,
		-10.39968914	,
		-10.18939869	,
		-10.059426	,
		-9.807356177	,
		-9.590216497	,
		-9.484780934	,
		-9.261140861	,
		-9.061201224	,
		-8.947455279	,
		-8.744426905	,
		-8.621018745	,
		-8.467026661	,
		-8.286252727	,
		-8.197367961	,
		-7.972002002	,
		-7.900493958	,
		-7.837990825	,
		-7.788273221	,
		-7.742254277	,
		-7.762619964	,
		-7.73348038	,
		-7.681238998	,
		-7.843990546	,
		-7.929511353	,
		-8.017964377	,
		-8.264732053	,
		-8.546703363	,
		-8.76210215	,
		-9.146652339	,
		-9.502983477	,
		-9.967046843	,
		-10.38031681	,
		-10.78983224	,
		-11.12533948	,
		-11.27882703	,
		-11.1045274	,
		-10.62400828	,
		-9.861018729	,
		-8.691543294	,
		-7.474412535	,
		-6.261170392	,
		-5.030578599	,
		-3.720697438	,
		-2.561968402	,
		-1.474000879	,
		-0.380559178	,
		0.603758844	,
		1.561558416	,
		2.49488288	,
		3.36346343	,
		4.202505585	,
		4.992308919	,
		5.746078823	,
		6.464438382	,
		7.161709806	,
		7.823329146	,
		8.449950557	,
		9.065381396	,
		9.645531415	,
		10.16900143	,
		10.67739282	,
		11.14420941	,
		11.58800746	,
		12.01864179	,
		12.41118186	,
		12.77361762	,
		13.11430657	,
		13.41433608	,
		13.69588858	,
		13.93904836	,
		14.16811697	,
		14.3632429	,
		14.51782441	,
		14.65412356	,
		14.76043383	,
		14.82918621	,
		14.86826383	,
		14.88587236	,
		14.86306297	,
		14.81834842	,
		14.74840942	,
		14.6509433	,
		14.52126513	,
		14.36824444	,
		14.18295479	,
		13.97069071	,
		13.68281313	,
		13.40335409	,
		13.10353146	,
		12.78298071	,
		12.41383029	,
		12.03362585	,
		11.61972846	,
		11.17317662	,
		10.68716003	,
		10.171562	,
		9.602890661	,
		9.030352029	,
		8.377054218	,
		7.710301574	,
		6.954274739	,
		6.156688304	,
		5.344519407	,
		4.452408253	,
		3.481320173	,
		2.368232653	,
		1.135754856	,
		-0.180329552	,
		-1.757941411	,
		-3.481015103	,
		-5.548070597	,
		-8.088342596	,
		-11.29022873	,
		-15.75613157	,
		-21.93392041	,
		-20.71391084	,
		-15.42019526	,
		-12.29120372	,
		-10.03584387	,
		-8.410781385	,
		-7.189434067	,
		-6.242551899	,
		-5.471015888	,
		-4.913654096	,
		-4.447639815	,
		-4.136636736	,
		-3.857638908	,
		-3.644172247	,
		-3.534183646	,
		-3.52020524	,
		-3.485185072	,
		-3.475468434	,
		-3.593456628	,
		-3.701545864	,
		-3.785124334	,
		-3.996223882	,
		-4.152621522	,
		-4.340666153	,
		-4.528399722	,
		-4.739755147	,
		-5.008493788	,
		-5.210899615	,
		-5.443121124	,
		-5.666307113	,
		-5.940500937	,
		-6.103407422	,
		-6.273201641	,
		-6.387609689	,
		-6.590655911	,
		-6.70115388	,
		-6.751618849	,
		-6.860985337	,
		-6.999121878	,
		-7.108487043	,
		-7.235598728	,
		-7.33418974	,
		-7.461973879	,
		-7.647059413	,
		-7.855189027	,
		-8.149350689	,
		-8.408623198	,
		-8.688086745	,
		-9.09871628	,
		-9.457396797	,
		-9.811411746	,
		-10.22100102	,
	};

	double phi_cut [361] = { 
		14.73718406	,
		14.69843592	,
		14.62957178	,
		14.54345487	,
		14.43179798	,
		14.28928172	,
		14.12674011	,
		13.95172908	,
		13.73782171	,
		13.51165673	,
		13.24963634	,
		12.97423666	,
		12.68391096	,
		12.36928852	,
		12.02796769	,
		11.67329489	,
		11.28734084	,
		10.89716257	,
		10.4527037	,
		9.958814147	,
		9.437241024	,
		8.906252595	,
		8.322718747	,
		7.6930477	,
		7.039393652	,
		6.322868387	,
		5.561212649	,
		4.752275128	,
		3.881427654	,
		2.954071509	,
		1.962203964	,
		0.890655686	,
		-0.236301658	,
		-1.461363717	,
		-2.69250004	,
		-4.086029704	,
		-5.444558329	,
		-6.67523491	,
		-7.909830969	,
		-8.721697099	,
		-9.102865593	,
		-9.132091069	,
		-8.905260633	,
		-8.503470959	,
		-8.00614114	,
		-7.594897253	,
		-7.088955057	,
		-6.788791154	,
		-6.42977486	,
		-6.222894914	,
		-6.025389739	,
		-5.834830796	,
		-5.77763714	,
		-5.74796663	,
		-5.629687885	,
		-5.663351046	,
		-5.722087956	,
		-5.805888506	,
		-5.772653361	,
		-5.871088447	,
		-5.936882895	,
		-6.056121134	,
		-6.118579247	,
		-6.186104002	,
		-6.302816957	,
		-6.383301717	,
		-6.560726599	,
		-6.623718661	,
		-6.795716378	,
		-6.922045555	,
		-7.024810047	,
		-7.258548809	,
		-7.350886504	,
		-7.54109745	,
		-7.760838212	,
		-7.947046036	,
		-8.103032005	,
		-8.351424631	,
		-8.610908972	,
		-8.859622948	,
		-9.151186211	,
		-9.388702568	,
		-9.722899845	,
		-9.946148132	,
		-10.29255121	,
		-10.54964297	,
		-10.7591457	,
		-11.01462398	,
		-11.16132439	,
		-11.4686834	,
		-11.59883401	,
		-11.72271425	,
		-11.84463549	,
		-11.94490179	,
		-11.97876049	,
		-12.1182655	,
		-12.16568388	,
		-12.41003321	,
		-12.43711813	,
		-12.50066547	,
		-12.62212231	,
		-12.60099213	,
		-12.70690602	,
		-12.77686522	,
		-12.71223979	,
		-12.63445517	,
		-12.63559666	,
		-12.71184064	,
		-12.69600107	,
		-12.67091253	,
		-12.7074972	,
		-12.92340102	,
		-13.03006314	,
		-13.38261806	,
		-13.59852021	,
		-14.10286922	,
		-14.54147991	,
		-15.04246184	,
		-15.59004844	,
		-16.30819011	,
		-16.81569057	,
		-17.49762076	,
		-17.6437648	,
		-18.17929508	,
		-18.25159751	,
		-18.16908471	,
		-18.32160994	,
		-18.29823335	,
		-18.22458588	,
		-18.45399169	,
		-18.62274033	,
		-18.83702891	,
		-19.4817684	,
		-19.76218491	,
		-20.34904768	,
		-21.07038827	,
		-21.19759954	,
		-21.12023562	,
		-20.72848107	,
		-19.96240102	,
		-19.30742428	,
		-18.50952518	,
		-17.77572864	,
		-17.23046838	,
		-16.59919665	,
		-16.34818835	,
		-16.3588671	,
		-16.30712106	,
		-16.6133163	,
		-16.77642172	,
		-17.46662356	,
		-18.24627483	,
		-18.94934046	,
		-20.12119884	,
		-21.43118904	,
		-22.31230727	,
		-22.61270223	,
		-22.79110526	,
		-22.12840987	,
		-20.52503447	,
		-19.53221233	,
		-18.49228506	,
		-17.68556129	,
		-16.82786269	,
		-16.42892177	,
		-15.92269034	,
		-15.64478827	,
		-15.69140996	,
		-15.59153969	,
		-15.59001921	,
		-15.7385851	,
		-16.03296822	,
		-16.29187887	,
		-16.7867741	,
		-17.00414408	,
		-17.38734666	,
		-17.52889353	,
		-17.64766232	,
		-17.59623611	,
		-17.44377943	,
		-16.91188084	,
		-17.14868542	,
		-16.72011126	,
		-16.41429065	,
		-15.71071658	,
		-15.28804086	,
		-14.8409156	,
		-14.18916502	,
		-13.99462945	,
		-13.77228295	,
		-13.55044712	,
		-13.38631327	,
		-13.29717267	,
		-13.251568	,
		-13.48524819	,
		-13.59395167	,
		-13.73261976	,
		-14.08240382	,
		-14.46490699	,
		-14.95575998	,
		-15.52307053	,
		-16.26580028	,
		-17.19286929	,
		-17.95550189	,
		-19.19272439	,
		-20.46327161	,
		-21.58601999	,
		-22.93913353	,
		-23.97829742	,
		-24.04179945	,
		-23.39359953	,
		-22.35436424	,
		-21.51245286	,
		-20.40670699	,
		-19.44420838	,
		-18.66600647	,
		-17.82313395	,
		-17.17396711	,
		-16.68730712	,
		-16.13728874	,
		-15.77977111	,
		-15.5768849	,
		-15.27179531	,
		-15.20120031	,
		-15.11725349	,
		-15.10575813	,
		-14.8574103	,
		-14.79981244	,
		-14.94936481	,
		-15.07982341	,
		-15.12038831	,
		-15.24175486	,
		-15.44679397	,
		-15.56029301	,
		-15.60091476	,
		-15.86541355	,
		-15.88321274	,
		-15.79358186	,
		-15.93469151	,
		-15.86791404	,
		-15.84222745	,
		-15.84741225	,
		-15.5771492	,
		-15.53838878	,
		-15.31029497	,
		-15.17213018	,
		-14.98088871	,
		-14.76807446	,
		-14.44931945	,
		-14.25240611	,
		-14.12456776	,
		-13.85317317	,
		-13.69798981	,
		-13.41491579	,
		-13.31515586	,
		-13.01618698	,
		-12.83144022	,
		-12.60563196	,
		-12.52755382	,
		-12.20351998	,
		-12.11385315	,
		-11.765205	,
		-11.61858614	,
		-11.39807496	,
		-11.08205316	,
		-10.86633146	,
		-10.52731418	,
		-10.18314187	,
		-9.855333061	,
		-9.492964218	,
		-9.238799711	,
		-8.881450486	,
		-8.648861237	,
		-8.266037277	,
		-8.032711072	,
		-7.763121738	,
		-7.604381614	,
		-7.316906376	,
		-7.106330908	,
		-6.881834905	,
		-6.796506769	,
		-6.607052829	,
		-6.448059232	,
		-6.272701147	,
		-6.17420561	,
		-5.971979744	,
		-5.815572572	,
		-5.598310554	,
		-5.386061213	,
		-5.106682347	,
		-4.88361893	,
		-4.625775563	,
		-4.244839533	,
		-3.945211943	,
		-3.679941308	,
		-3.361624337	,
		-3.025034765	,
		-2.75406537	,
		-2.445441721	,
		-2.179190105	,
		-1.965173649	,
		-1.729327502	,
		-1.503129627	,
		-1.307802001	,
		-1.145639812	,
		-1.010866134	,
		-0.919509941	,
		-0.841841821	,
		-0.824779822	,
		-0.883599711	,
		-0.902436561	,
		-1.029536398	,
		-1.216587736	,
		-1.444616473	,
		-1.722613741	,
		-2.088416578	,
		-2.511481459	,
		-3.025236247	,
		-3.678514051	,
		-4.277667045	,
		-5.065559543	,
		-5.809837011	,
		-6.51944613	,
		-6.926722879	,
		-6.910759994	,
		-6.36364292	,
		-5.347434521	,
		-4.077605718	,
		-2.766231351	,
		-1.461884007	,
		-0.198152668	,
		1.003686519	,
		2.124581361	,
		3.163289402	,
		4.10403568	,
		4.883104762	,
		5.665509001	,
		6.483984329	,
		7.217746645	,
		7.924525974	,
		8.557831883	,
		9.178653669	,
		9.766693979	,
		10.30495637	,
		10.82863313	,
		11.30739708	,
		11.7495846	,
		12.17374829	,
		12.54852953	,
		12.89037277	,
		13.2059179	,
		13.5574646	,
		13.82080194	,
		14.03957781	,
		14.2287765	,
		14.39104082	,
		14.51278943	,
		14.61719585	,
		14.69115831	,
		14.7155619	,
		14.74283455	,
	};

	double incr = 1;
	int theta_p = (int)(theta_d/incr);
	int phi_p = (int)(phi_d/incr);

	int index1 = theta_p;
	int index2 = index1 +  1;

	int index3 = phi_p;
	int index4 = index3 + 1;

	// The distances of the point from the vertices are the weights
	double w1 = sqrt(pow(theta_d - theta_p,2) + pow(phi_d - 0,2));
	//	double w2 = index2*incr - theta_p;	
	double w3 = sqrt(pow(theta_d - 90,2) + pow(phi_d - phi_p,2));

	double g1 = pow(10, (theta_cut[index1])/20);
	double g2 = pow(10, (theta_cut[index2])/20);

	double g3 = pow(10, (phi_cut[index3])/20);
	double g4 = pow(10, (phi_cut[index4])/20);

	//	return (w1*g1 + w2*g2)/(w1 + w2);
	return (w1*g3 + w3*g1)/(w1 + w3);

}

double gain_theta_rx
(	
	double theta, 
	double phi,
	double theta_cut[],
	double phi_cut[]
)
{
	double theta_d = theta * 180/PI;
	double phi_d = phi * 180/PI;


	double incr = 1;
	int theta_p = (int)(theta_d/incr);
	int phi_p = (int)(phi_d/incr);

	int index1 = theta_p;
	int index2 = index1 +  1;

	int index3 = phi_p;
	int index4 = index3 + 1;

	// The distances of the point from the vertices are the weights
	// double w1 = sqrt(pow(theta_d - theta_p,2) + pow(phi_d - 0,2));
	//	double w2 = index2*incr - theta_p;	
	// double w3 = sqrt(pow(theta_d - 90,2) + pow(phi_d - phi_p,2));

//	double g1 = pow(10, (theta_cut[index1])/20);
	// double g2 = pow(10, (theta_cut[index2])/20);

//	double g3 = pow(10, (phi_cut[index3]-8.415267184999999)/20); // normalize max beam value to 1 (8.415267184999999 is max of phi_cut).
	// double g4 = pow(10, (phi_cut[index4])/20);

	//	return (w1*g1 + w2*g2)/(w1 + w2);
	
	double w1 = delta (theta_d, phi_d, theta_p, phi_p);
	double w2 = delta (theta_d, phi_d, theta_p, phi_p + incr);
	double w3 = delta (theta_d, phi_d, theta_p + incr, phi_p);
	double w4 = delta (theta_d, phi_d, theta_p + incr, phi_p + incr);

/*	double g1 = pow(10, (antenna_pattern[index1])/20);
	double g2 = pow(10, (antenna_pattern[index2])/20);
	double g3 = pow(10, (antenna_pattern[index3])/20);
	double g4 = pow(10, (antenna_pattern[index4])/20);
*/

	double g1 = theta_cut[index1];
	double g2 = theta_cut[index2];
	double g3 = phi_cut[index3];
	double g4 = phi_cut[index4];
	
	return pow(10,(w1*g1 + w2*g2 + w3*g3 + w4*g4)/(w1 + w2 + w3 + w4)*(1.0/20.0));
//	return g1*g3;

}


//double gain (double antenna_pattern[], double theta, double phi)
//{
//	const int num_phi = 64;
//	const int num_theta = 33;
//	const int num_vals = num_phi*num_theta;
//
//	double theta_d = theta * 180/PI;
//	double phi_d = phi * 180/PI;
//
//	double incr = 5.625;
//	int theta_p = (int)(theta_d/incr);
//	int phi_p = (int)(phi_d/incr);
//
//	int index1 = theta_p * num_phi + phi_p;
//	int index2 = index1 +  1;
//	int index3 = (theta_p + 1) * num_phi + phi_p;
//	int index4 = index3 + 1;
//
//	// The distances of the point from the 4 vertices are the weights
//	double w1 = delta (theta_d, phi_d, theta_p, phi_p);
//	double w2 = delta (theta_d, phi_d, theta_p, phi_p + incr);
//	double w3 = delta (theta_d, phi_d, theta_p + incr, phi_p);
//	double w4 = delta (theta_d, phi_d, theta_p + incr, phi_p + incr);
//
//	double g1 = pow(10, (antenna_pattern[index1])/20);
//	double g2 = pow(10, (antenna_pattern[index2])/20);
//	double g3 = pow(10, (antenna_pattern[index3])/20);
//	double g4 = pow(10, (antenna_pattern[index4])/20);
//
//	return (w1*g1 + w2*g2 + w3*g3 + w4*g4)/(w1 + w2 + w3 + w4);
//}

double delta (double a1, double b1, double a2, double b2)
{
	return sqrt(pow((a1 - a2),2) + pow((b1 - b2),2));
}

void  read_pattern_files 
	(
	double theta[], 
	double phi[]
)
{
	string fname_phi, fname_theta;
	string phi_file_ext = " GHz_gain_phi.txt";
	string theta_file_ext = " GHz_gain_theta.txt";

	for (int i = 0; i < NUM_FREQ; i++)
	{
		string freq_str = freq_index_to_string (i);
		fname_phi = freq_str + phi_file_ext;		
		fname_theta = freq_str + theta_file_ext;
		read_pattern_for_freq (phi, fname_phi, i);
		read_pattern_for_freq (theta, fname_theta, i);
	}
}

// Adds data for the missing half space
// where phi varies from 90.5 - 269.5
void  read_pattern_files_nff
(
	double theta[], 
	double phi[],
	string fname,
	const int size
)
{
	ifstream antenna_pattern;
	antenna_pattern.open(fname.c_str());
	assert (!antenna_pattern.fail());
	string s;
	int count = -1;
	int phi_count = 1;
	int first_space, second_space, third_space;
	string theta_val, s_remain, phi_val, s_rem_remain, gain, gain_phi;
	while (!antenna_pattern.eof())
	{
		getline(antenna_pattern, s);
		if (count >= 0)
		{
			// Get gain values
			// <format>:theta phi Etheta_mag Ephi_mag
			// <example>:0.0 0.0 0.0 0.0
			first_space = string::npos;
			first_space = s.find_first_of(' ');
			theta_val = s.substr(0, first_space);

			s_remain = s.substr(first_space + 1);

			second_space = string::npos;
			second_space = s_remain.find_first_of(' ');
			phi_val = s_remain.substr(0, second_space);

			s_rem_remain = s_remain.substr(second_space + 1);

			third_space = string::npos;
			third_space = s_rem_remain.find_first_of(' ');
			gain = s_rem_remain.substr(0, third_space);

			//string phase = s_rem_remain.substr(third_space + 1);
			gain_phi = s_rem_remain.substr(third_space + 1);


			if (count < size)
			{
				// The gain reported is 20 log10 (gain)
				// Therefore the magnitude of the gain is 10^(gain_from_file/20)
				//double temp_gain = pow(10,(atof(gain.c_str())/20));
				double temp_gain = atof(gain.c_str());
				theta [count] = temp_gain;
				temp_gain = atof(gain_phi.c_str());
				phi[count] = temp_gain;
				//			cout << theta_val << "\t\t" << phi_val << "\t\t" << gain << "\t\t" << phase << "\t\t" << count << "\n";
			}
		}
		count++;
	}
	antenna_pattern.close();

}

void  read_pattern_for_freq
	(
	double gain_data[],
	string fname,
	int freq_index
	)
{
	ifstream antenna_pattern;
	antenna_pattern.open(fname.c_str());
	assert (!antenna_pattern.fail());
	string s;
	int count = -1;
	while (!antenna_pattern.eof())
	{
		getline(antenna_pattern, s);
		if (count >= 0)
		{
			// Get gain values
			// <format>:theta phi gain
			// <example>:0.0 0.0 24
			int first_space = string::npos;
			first_space = s.find_first_of(' ');
			string theta = s.substr(0, first_space);

			string s_remain = s.substr(first_space + 1);

			int second_space = string::npos;
			second_space = s_remain.find_first_of(' ');
			string phi = s_remain.substr(0, second_space);

			string gain = s_remain.substr(second_space + 1);

			if (count < NUM_VALS)
			{
				// The gain reported is 20 log10 (gain)
				// Therefore the magnitude of the gain is 10^(gain_from_file/20)
				//double temp_gain = pow(10,(atof(gain.c_str())/20));
				double temp_gain = atof(gain.c_str());
				gain_data [freq_index*NUM_VALS + count] = temp_gain;
				//	cout << theta << "\t\t" << phi << "\t\t" << gain << "\t\t" << count << "\t\t" << temp_gain << "\n";
			}
		}
		count++;
	}
	antenna_pattern.close();

}

string freq_index_to_string (int freq_index)
{
	/*Map array index to frequency
	Index	Frequency
	0		2
	1		2.5
	2		3
	3		3.5
	4		4
	5		4.5
	6		5
	7		5.5
	8		6
	9		6.5
	10		7		*/	

	string freq_string = "";
	switch (freq_index)
	{
	case 0:
		freq_string = "2.0";
		break;
	case 1:
		freq_string = "2.5";
		break;
	case 2:
		freq_string = "3.0";
		break;
	case 3:
		freq_string = "3.5";
		break;
	case 4:
		freq_string = "4.0";
		break;
	case 5:
		freq_string = "4.5";
		break;
	case 6:
		freq_string = "5.0";
		break;
	case 7:
		freq_string = "5.5";
		break;
	case 8:
		freq_string = "6.0";
		break;
	case 9:
		freq_string = "6.5";
		break;
	case 10:
		freq_string = "7.0";
		break;
	default:
		freq_string = "";
	}
	return freq_string;
}

int freq_to_freq_index (double freq)
{
	int freq_index = -1;

	if (freq < 2.25e9)	// 2 GHz file
	{
		freq_index = 0;	
	}
	else if ( freq > 2.25e9 && freq < 2.75e9) // 2.5 GHz file
	{
		freq_index = 1;
	}
	else if (freq > 2.75e9 && freq < 3.25e9) // 3 GHz file
	{
		freq_index = 2;
	}
	else if (freq > 3.25e9 && freq < 3.75e9) // 3.5 GHz file
	{
		freq_index = 3;
	}
	else if (freq > 3.75e9 && freq < 4.25e9) // 4 GHz file
	{
		freq_index = 4;
	}
	else if (freq > 4.25e9 && freq < 4.75e9) // 4.5 GHz file
	{
		freq_index = 5;
	}
	else if (freq > 4.75e9 && freq < 5.25e9) // 5 GHz file
	{
		freq_index = 6;
	}
	else if (freq > 5.25e9 && freq < 5.75e9) // 5.5 GHz file
	{
		freq_index = 7;
	}
	else if (freq > 5.75e9 && freq < 6.25e9) // 6 GHz file
	{
		freq_index = 8;
	}
	else if (freq > 6.25e9 && freq < 6.75e9) // 6.5 GHz file
	{
		freq_index = 9;
	}
	else if (freq > 6.75e9) // 7 GHz file
	{
		freq_index = 10;
	}
	else	// Should never be reached
	{
		freq_index = -1;
	}

	return freq_index;
}

void setup_cube 
	(
	cube &c,
	double xmin, 
	double xmax, 
	double ymin, 
	double ymax, 
	double zmin, 
	double zmax
	)
{
	// First compute the 8 points that form the vertices
	// of the cube
	//			   p7_________p8
	//             /|		  /
	//            /	|		 / |	
	//			p5__|_______p6 |	
	//			|	|		|  |	
	//			|  p3_______|__p4
	//			| /			| /
	//			|/__________|/
	//			p1			p2
	//
	// p1 (xmin, ymin, zmin)
	// p2 (xmax, ymin, zmin)
	// p3 (xmin, ymax, zmin)
	// p4 (xmax, ymax, zmin)
	// p5 (xmin, ymin, zmax)
	// p6 (xmax, ymin, zmax)
	// p7 (xmin, ymax, zmax)
	// p8 (xmax, ymax, zmax)

	point p1, p2, p3, p4, p5, p6, p7, p8;
	init_point (p1, xmin, ymin, zmin);
	init_point (p2, xmax, ymin, zmin);
	init_point (p3, xmin, ymax, zmin);
	init_point (p4, xmax, ymax, zmin);
	init_point (p5, xmin, ymin, zmax);
	init_point (p6, xmax, ymin, zmax);
	init_point (p7, xmin, ymax, zmax);
	init_point (p8, xmax, ymax, zmax);

	// Bottom
	init_finite_plane_lite (c.finite_bottom , p1, p2, p4, p3);
	c.bottom.d = c.finite_bottom.d;
	c.bottom.unit_normal = c.finite_bottom.unit_normal;

	// Top
	init_finite_plane_lite (c.finite_top, p5, p6, p8, p7);
	c.top.d = c.finite_top.d;
	c.top.unit_normal = c.finite_top.unit_normal;

	// Near
	init_finite_plane_lite (c.finite_near, p1, p2, p6, p5);
	c.near.d = c.finite_near.d;
	c.near.unit_normal = c.finite_near.unit_normal;

	// Far
	init_finite_plane_lite (c.finite_far, p3, p4, p8, p7);
	c.far.d = c.finite_far.d;
	c.far.unit_normal = c.finite_far.unit_normal;

	// Left
	init_finite_plane_lite (c.finite_left, p1, p3, p7, p5);
	c.left.d = c.finite_left.d;
	c.left.unit_normal = c.finite_left.unit_normal;

	// Right
	init_finite_plane_lite (c.finite_right, p2, p4, p8, p6);
	c.right.d = c.finite_right.d;
	c.right.unit_normal = c.finite_right.unit_normal;	
}

/*
void init_finite_plane_lite (finite_plane_lite& fp, point_list* v)
{
	point_list* iter = v;
	while (iter != NULL)
	{
	}
}
*/

void init_finite_plane_lite
	(
	finite_plane_lite & fp,
	// Points pa, pb, pc, pd are in cyclic or acyclic order 
	/*
	pb------pc
	|		|						
	|		|					
	pa------pd
	*/
	point pa,
	point pb,
	point pc,
	point pd
	)
{
	fp.num_vert = 4;	// Fixed for now to make things simple
	fp.v = new point_list;
	fp.v->p = pa;
	fp.v->next = new point_list;
	fp.v->next->p = pb;
	fp.v->next->next = new point_list;
	fp.v->next->next->p = pc;
	fp.v->next->next->next = new point_list;
	fp.v->next->next->next->p = pd;
	fp.v->next->next->next->next = NULL;
	
	// Determine the normal to the plane and the value of the constant d in the equation
	// ax + by + cz + d = 0
	plane_normal (fp);
}

void init_finite_plane_lite
	(
	finite_plane_lite & fp,
	// Points pa, pb, pc are in cyclic or acyclic order (always satisfies for 3 points)
	/*
	pb---pc
	|   /								
	|  /
	| /
	|/
	pa
	*/
	point pa,
	point pb,
	point pc
	)
{
	fp.num_vert = 3;	// Fixed for now to make things simple
	fp.v = new point_list;
	fp.v->p = pa;
	fp.v->next = new point_list;
	fp.v->next->p = pb;
	fp.v->next->next = new point_list;
	fp.v->next->next->p = pc;
	fp.v->next->next->next = NULL;	

	// Determine the normal to the plane and the value of the constant d in the equation
	// ax + by + cz + d = 0
	plane_normal (fp);
}


void init_finite_plane_lite
	(
	finite_plane_lite3 & fp,
	// Points pa, pb, pc are in cyclic or acyclic order (always satisfies for 3 points)
	/*
	pb---pc
	|   /								
	|  /
	| /
	|/
	pa
	*/
	point pa,
	point pb,
	point pc
	)
{
	fp.num_vert = 3;	// Fixed for now to make things simple
	fp.vl[0] = pa;
	fp.vl[1] = pb;
	fp.vl[2] = pc;
	
	// Determine the normal to the plane and the value of the constant d in the equation
	// ax + by + cz + d = 0
	plane_normal (fp);
}

void plane_normal (finite_plane_lite & fp)
{
	point pa = fp.v->p;
	point pb = fp.v->next->p;
	point pc = fp.v->next->next->p;

	point temp = pb - pa;
	vect n1 = point2vect (temp);
	temp = pc - pb;
	vect n2 = point2vect (temp);

	vect n = cross (n1, n2);
	n = normalize (n);
	fp.unit_normal = n;
	fp.d = -1*n*pa;
}

void plane_normal (finite_plane_lite3 & fp)
{
	point pa = fp.vl[0];
	point pb = fp.vl[1];
	point pc = fp.vl[2];

	point temp = pb - pa;
	vect n1 = point2vect (temp);
	temp = pc - pb;
	vect n2 = point2vect (temp);

	vect n = cross (n1, n2);
	n = normalize (n);
	fp.unit_normal = n;
	fp.d = -1*n*pa;
}

bool intersection_finite_plane_pyramid (finite_plane p, rect_pyramid py)
{
	double t_val;
	double s_val;
	bool intersect;
	ray* ray_list = new ray[p.num_vert]; // TODO: Convert to static
	//ray ray_list[4];
	point_list* iter;
	for (int i = 0; i < p.num_vert; i++)
	{
		iter = p.v;
		for (int k = 0; k < i; k++)
		{
			iter = iter->next;	// NSOOD: Add error handling
		}

		ray_list[i].orig = py.apex;
		point dir = iter->p - py.apex;
		ray_list[i].dir = point2vect(dir);

		intersect = intersection_ray_finite_plane_lite (py.base , ray_list[i], t_val);

		if (intersect)
		{
			// At least one vertex of p lies within the pyramid
			// Delete ray_list before returning
			delete[] ray_list;
			return true;
		}
	}
	
	// No vertex of p lies within the pyramid
	// Line joining the apex and middle of the diagonal of the base (pa + pc)/2
	/*
	pb------pc
	|		|						
	|		|					
	pa------pd
	*/

	point* bvl = new point[py.base.num_vert];
	iter = py.base.v;
	for (int i=0; i < py.base.num_vert; i++)
	{
		bvl[i] = iter->p;
		iter = iter->next;
	}

	// NSOOD: Modified, since the old approach would not work for triangles
	// A polygon would have at least 3 vertices
	// Hence use the first 3 vertices
	ray mid_ray;
	mid_ray.orig = py.apex;
	point temp = (bvl[0] + bvl[1])*0.5;		// Mid-point of line-segment joining first and second vertex
	point mid_pt = (temp + bvl[2])*0.5;
	point dir = mid_pt - py.apex;
	mid_ray.dir = point2vect(dir);

	intersect = intersection_ray_finite_plane_lite (p, mid_ray, t_val);
	if (intersect)
	{
		// Plane crosses through the pyramid
		// Delete ray_list and bvl before returning
		delete[] ray_list;
		delete[] bvl;
		return true;
	}

	// Plane either is completely outside of the pyramid or crosses through
	// without covering the entire cross section
	// Find intersection of each line with each of the 4 faces of the pyramid
	// Don't need to do this for the base yet.

	point* vl = new point[p.num_vert];
	iter = p.v;
	for (int i=0; i < p.num_vert; i++)
	{
		vl[i] = iter->p;
		iter = iter->next;
	}

	for (int i = 0; i < p.num_vert; i++)
	{
		ray_list[i].orig = vl[i];
		point dir = vl[(i+1)%p.num_vert] - vl[i];
		ray_list[i].dir = point2vect(dir);

		fp_lite_list_node* iterf = py.face_list;
		for (int j = 0; j < py.base.num_vert; j++)
		{
			intersect = intersection_ray_inf_finite_plane_lite (iterf->fp, ray_list[i], t_val);			
			if (intersect && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO))
		//	if (intersect && (t_val < BIT_MORE_THAN_ONE && t_val > BIT_LESS_THAN_ZERO))
			{
				// Find intersection between the line joining the apex and the potential point of intersection
				// and the line joining the 2 points which when combined with the apex form the triangular face of the pyramid
				// I really didn't explain that well, did I ... ?
				intersect = intersection_line_seg_line_seg(py.apex, ray_list[i].orig+t_val*ray_list[i].dir, iterf->fp.v->next->p, iterf->fp.v->next->next->p, t_val, s_val);

				if (intersect && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO) && (s_val < BIT_LESS_THAN_ONE && s_val > BIT_MORE_THAN_ZERO))
				//if (intersect && (t_val < BIT_MORE_THAN_ONE && t_val > BIT_LESS_THAN_ZERO) && (s_val < BIT_MORE_THAN_ONE && s_val > BIT_LESS_THAN_ZERO))
				{
					// Delete ray_list before returning
					delete[] ray_list;
					delete[] bvl;
					delete[] vl;
					return true;	// If any one line intersects with any plane
				}
			}
			iterf = iterf->next;
		}
	}

	// Delete ray_list before returning
	delete[] ray_list;
	delete[] bvl;
	delete[] vl;
	return false; // No intersection was found there return false

}

bool intersection_finite_plane_ext_zone (finite_plane p, rect_pyramid py)
{
	/*	double t_val;
	bool intersect;
	//ray* ray_list = new ray[p.num_vert]; // TODO: Convert to static
	ray ray_list[4];
	point dir;

	point* vl = new point[p.num_vert];
	point_list* iter = p.v;
	for (int i=0; i < p.num_vert; i++)
	{
		vl[i] = iter->p;
		iter = iter->next;
	}

	// If the plane intersects with the base of the pyramid
	// Then false should be returned so that the node is included
	for (int i = 0; i < p.num_vert; i++)
	{
		ray_list[i].orig = vl[i];
		dir = vl[(i+1)%p.num_vert] - vl[i];
		ray_list[i].dir = point2vect(dir);

		intersect = intersection_ray_finite_plane_lite (py.base , ray_list[i], t_val);
		if (intersect && (t_val < 1))
		{
			return false;
		}
		
		ray_list[i].orig = py.base.vl[i];
		dir = py.base.vl[(i+1)%py.base.num_vert] - py.base.vl[i];
		ray_list[i].dir = point2vect(dir);

		intersect = intersection_ray_finite_plane_lite (p, ray_list[i], t_val);
		if (intersect && (t_val < 1))
		{
			return false;
		}

	}

	for (int i = 0; i < p.num_vert; i++)
	{
		ray_list[i].orig = py.apex;
		dir = vl[i] - py.apex;
		ray_list[i].dir = point2vect(dir);

		intersect = intersection_ray_finite_plane_lite (py.base , ray_list[i], t_val);
		if (t_val <= 1)	// Intersection happens outside of the pyramid
		{
			intersect = false;
		}

		if (intersect)
		{
			// At least one vertex of p lies within the pyramid			
//			cout << "\t\t\tIntersection happens within the pyramid" << endl;
//			cout << "\t\t\tIntersecting ray formed with vertex: " << i << endl;
//			cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
			return true;
		}
	}
	// No vertex of p lies within the pyramid
	// Line joining the apex and middle of the diagonal of the base (pa + pc)/2
	//
	//  pb------pc
	//  |		|						
	//  |		|					
	//  pa------pd
	
	ray mid_ray;
	mid_ray.orig = py.apex;
	point mid_pt = (py.base.vl[0] + py.base.vl[2])*0.5;
	dir = mid_pt - py.apex;
	mid_ray.dir = point2vect(dir);

	intersect = intersection_ray_finite_plane_lite (p , mid_ray, t_val);
	if (t_val >= 1)	// Intersection happens outside of the pyramid
	{
		intersect = false;
	}
	if (intersect)
	{
		// Plane crosses through the pyramid		
//		cout << "\t\t\tPlane crosses through the pyramid" << endl;		
//		cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
		return true;
	}

	// Plane either is completely outside of the pyramid or crosses through
	// without covering the entire cross section
	// Find intersection of each line with each of the 4 faces of the pyramid
	
	for (int i = 0; i < p.num_vert; i++)
	{
		ray_list[i].orig = vl[i];
		dir = vl[(i+1)%p.num_vert] - vl[i];
		ray_list[i].dir = point2vect(dir);

		intersect = intersection_ray_finite_plane_lite (py.face1 , ray_list[i], t_val);
		if ((intersect) && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO))
		{
//			cout << "\t\t\tIntersection happens with Face 1" << endl;
//			cout << "\t\t\tIntersecting ray formed with vertices: " << i << " " << (i+1)%p.num_vert << endl;
//			cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
			return true;	// If any one line intersects with any plane
		}
		intersect = intersection_ray_finite_plane_lite (py.face2 , ray_list[i], t_val);
		if ((intersect) && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO))
		{
//			cout << "\t\t\tIntersection happens with Face 2" << endl;
//			cout << "\t\t\tIntersecting ray formed with vertices: " << i << " " << (i+1)%p.num_vert << endl;
//			cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
			return true;	// If any one line intersects with any plane
		}
		intersect = intersection_ray_finite_plane_lite (py.face3 , ray_list[i], t_val);
		if ((intersect) && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO)) 
		{
//			cout << "\t\t\tIntersection happens with Face 3" << endl;
//			cout << "\t\t\tIntersecting ray formed with vertices: " << i << " " << (i+1)%p.num_vert << endl;
//			cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
			return true;	// If any one line intersects with any plane
		}
		intersect = intersection_ray_finite_plane_lite (py.face4 , ray_list[i], t_val);
		if ((intersect) && (t_val < BIT_LESS_THAN_ONE && t_val > BIT_MORE_THAN_ZERO))
		{
//			cout << "\t\t\tIntersection happens with Face 4" << endl;
//			cout << "\t\t\tIntersecting ray formed with vertices: " << i << " " << (i+1)%p.num_vert << endl;
//			cout << "\t\t\tIntersecting ray parameter value: " << t_val << endl;
			return true;	// If any one line intersects with any plane
		}		
	}
*/

	return false; // No intersection was found there return false
}


void printHelp()
{
	fprintf(stderr, "raytracer:	raytracer.exe [ENVIRONMENT_FILE] [optional_arguments=value ...] \n");
	fprintf(stderr, "ENVIRONMENT_FILE must be speciied and is the path of a valid environment file. \n");
	fprintf(stderr, "-refs=INTEGER	 |	Overwrites the number of reflections specified in the environment file. \n");
	fprintf(stderr, "-grel_perm=FLOAT|	If the relative permativity of a surface is not given, it will use the global value. \n");
	fprintf(stderr, "-gsigma=FLOAT	 |	If the conductivity of a surface is not given, it will use the global value. \n");
	fprintf(stderr, "-grho_h=FLOAT	 |	If the roughness of a surface is not given, it will use the global value. \n");
	fprintf(stderr, "-rec_loc=STRING |	Used for distance sweeps. Must be a valid file path with correct receiver co-ordinates. \n");
	fprintf(stderr, "-debug			 |  Prints information on the surface sigma, epsilon. \n");
	return;

}

int check_points_on_plane(finite_plane plane)
{
	// http://mathworld.wolfram.com/Coplanar.html
	// (p3 - p1) DOT [(p2-p1) x (p4-p3)] = 0
	//	a			b			c
	point p1 = plane.v->p;
	point p2 = plane.v->next->p;
	point p3 = plane.v->next->next->p;
	point p4 = plane.v->next->next->next->p;

	// should technically be vectors
	point a = p3 - p1;
	point b = p2 - p1;
	point c = p4 - p3;

	// b x c
	point d;
	d.x = b.y * c.z - c.y * b.z;
	d.y = - (b.x * c.z - c.x * b.z);
	d.z = b.x * c.y - b.y * c.x;

	// do dot product
	double ans = a.x * d.x + a.y * d.y + a.z * d.z;

	if (fabs(ans) < 0.0001)
	{
		return(0);
	}
	else
	{
		return(-1);
	}
}

void printSurfaces(finite_plane* surf_list_fp, int num_surfaces)
{
	for(int i = 0 ; i < num_surfaces;i++)
	{
		cout << "Surface " << i << ":" << endl;
		cout << "conductivity = " << surf_list_fp[i].sigma << endl;
		cout << "rel_perm = " << surf_list_fp[i].rel_perm << endl << endl;
		cout << "roughness = " << surf_list_fp[i].rho_h << endl << endl;
	}

	return;
}

bool iszero(double f)
{
	if (f < SMALL_DOUBLE && f > -SMALL_DOUBLE)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool intersection_line_seg_line_seg (point p1, point p2, point p3, point p4, double& t_value, double& s_value)
{
/*	t(x2 - x1) + s(x3 - x4) = x3 - x1
	t(y2 - y1) + s(y3 - y4) = y3 - y1
	t(z2 - z1) + s(z3 - z4) = z3 - z1

	a00*t +	a01*s = a02
	a10*t + a11*s = a12
	a20*t + a21*s = a22
*/

	const int DIM = 3;
	double A[DIM][DIM] = {{p2.x - p1.x, p3.x - p4.x, p3.x - p1.x},
					      {p2.y - p1.y, p3.y - p4.y, p3.y - p1.y},
			              {p2.z - p1.z, p3.z - p4.z, p3.z - p1.z}}; 

	// Pivoting to ensure diagonal entries a00, a11 are non-zero
	double temp;
	if (iszero(A[0][0]))
	{
		if (iszero(A[1][0]))
		{			
			if (iszero(A[2][0]))
			{
				// Coefficient of t_values is zero in all 3 equations
				// No intersection exists
				return false;
			}
			else
			{
				// a00 is 0, a10 is 0 but a20 is not
				// -> Swap R0 and R2
				for (int i=0; i < DIM; i++)
				{
					temp = A[0][i];
					A[0][i] = A[2][i];
					A[2][i] = temp;
				}
			}
		}
		else
		{
			// a00 is 0 but a10 is not
			// -> Swap R0 and R1
			for (int i=0; i < DIM; i++)
			{
				temp = A[0][i];
				A[0][i] = A[1][i];
				A[1][i] = temp;
			}
		}
	}

	if (iszero(A[1][1]))
	{
		if (iszero(A[2][1]))
		{
			if (iszero(A[0][1]))
			{
				// Coefficient of s_values is zero in all equations
				// No intersection exists
				return false;
			}
			else
			{
				// a01 != 0
				if (!iszero(A[1][0]))
				{
					// a10 != 0
					// Swap R0 and R1
					for (int i=0; i < DIM; i++)
					{
						temp = A[0][i];
						A[0][i] = A[1][i];
						A[1][i] = temp;
					}
				}
				else 
				{	// a10 and a11 are 0, therefore for consistency a12 should be zero
					if (!iszero(A[1][2]))
					{
						// if a12 is not zero, some ain't right
						return false;
					}

					// if it is consistent then check that last equation
					if (!iszero(A[2][0]))
					{
						// a20 != 0
						// Swap R0 and R2
						for (int i=0; i < DIM; i++)
						{
							temp = A[0][i];
							A[0][i] = A[1][i];
							A[1][i] = temp;
						}
					}
					else
					{	// a20 and a21 are 0, therefor for consistency a22 should be zer0
						if (!iszero(A[2][2]))
						{
							// if a22 is not zero, some ain't right
							return false;
						}
					}
				}
			}
		}
		else
		{
			// a11 is 0 but a21 is not
			// -> Swap R1 and R2
			for (int i=0; i < DIM; i++)
			{
				temp = A[1][i];
				A[1][i] = A[2][i];
				A[2][i] = temp;
			}
		}
	}

	// If control reaches here then, by now a00 and a11 are non-zero
	// Perform Gaussian Elimination to make a10, a20 and a21 zero
	// a10	
	if (!iszero(A[1][0]))
	{
		// R1 <- R1 - (a10/a00)*R0
		temp = A[1][0];
		for (int i=0; i < DIM; i++)
		{				
			A[1][i] = A[1][i] - (temp/A[0][0])*A[0][i];
		}
	}

	// a20	
	if (!iszero(A[2][0]))
	{
		// R2 <- R2 - (a20/a00)*R0
		temp = A[2][0];
		for (int i=0; i < DIM; i++)
		{				
			A[2][i] = A[2][i] - (temp/A[0][0])*A[0][i];
		}
	}

	// a21	
	if (!iszero(A[2][1]))
	{
		// R2 <- R2 - (a21/a11)*R1
		temp = A[2][1];
		for (int i=0; i < DIM; i++)
		{				
			A[2][i] = A[2][i] - (temp/A[1][1])*A[1][i];
		}
	}
	
	// Check for consistency of the third row (R2)
	if (!iszero(A[2][2]))
	{	
		// Not consistent => Intersection does not exist
		return false;
	}

	// Compute s_value using the second row (R1)
	s_value = A[1][2]/A[1][1];

	// Compute t_value using the first row (R0) and s_value
	t_value = (A[0][2] - A[0][1]*s_value)/A[0][0];

	return true;	// Intersection has been computed

/*
	if (iszero(a11))
	{
		if (iszero(a12))
		{
			if (iszero(a13))
			{
				// Consistent
			}
			else
			{
				// Not consistent
				return false;
			}
		}
		else
		{
			// The coefficient of t is 0, but s is not
			s_value = a13/a12;
			if (a21 == 0)
			{
				// Can't solve for t here, but should check for consistentcy
				if(iszero(a23 - a22*s_value))
				{
					// Consistent
				}
				else
				{
					// Not consistent
					return false;
				}
			}
			else
			{
				t_value = (a23 - a22*s_value)/a21;
				// Now check consistency with the last equation
				if (iszero(a33 - a31*t_value - a32*s_value))
				{
					// t_value and s_value have been computed
					// values are consistent with all 3 equations
					// intersection exists
					return true
				}
				else
				{
					// Not consistent
					return false;
				}
			}			
		}
	}
	else // a11 is not zero
	{
		if (iszero(a12))
		{
			t_value = a13/a11;
			if (a21 == 0)
			{	

	}
*/		

}

void reduce_surface (surface & s, double fac)
{
	// Triangular surface
	if (s.num_vert == 3)
	{
		point p1 = s.v->p;
		point p2 = s.v->next->p;
		point p3 = s.v->next->next->p;

		point c = (1.0/3)*(p1+p2+p3);

		point p1s = (fac*p1+c)*(1/(1+fac));
		point p2s = (fac*p2+c)*(1/(1+fac));
		point p3s = (fac*p3+c)*(1/(1+fac));

		s.v->p = p1s;
		s.v->next->p = p2s;
		s.v->next->next->p = p3s;	
	}
}

double compute_att_fac(double a, double theta_i, double k, double bw, bool isTE)
{
	double dphi = 0.01;
	double phi = -PI/2;
	double integral_full = 0;
	double integral_spec_beam = 0;
	bw = bw*PIOVER180;
	double curr_term;
	
	double phiu = -theta_i + bw/2;
	double phil = -theta_i - bw/2;
	if (phil < -PI/2) {	phil = -PI/2; }
		
	while (phi <= PI/2)
	{
		if (isTE)	// TE or perpendicular polarization
		{
			curr_term = pow(sinc(0.5*k*a*(sin(phi)+sin(theta_i)))*cos(theta_i),2);
		}
		else	// TM or parallel polarization
		{
			curr_term = pow(sinc(0.5*k*a*(sin(phi)+sin(theta_i)))*cos(phi),2);
		}
		integral_full += curr_term;
		if (phi > phil && phi < phiu)
		{
			integral_spec_beam += curr_term;
		}
		phi = phi+dphi;
	}

	return sqrt(integral_spec_beam/integral_full);
}

double sinc(const double& x)
{
	if (x == 0)
		return 1;
	else
		return sin(x)/x;
}

double shortest_side(surface surf)
{
	double side_length, temp_length; 
	side_length = LARGE_DOUBLE;
	point p1, p2;
	point_list* iter = surf.v;
	point first = surf.v->p;
	
	while(iter != NULL)
	{		
		if (iter->next != NULL)
		{
			p1 = iter->p;
			p2 = iter->next->p;
		}
		else
		{
			p1 = iter->p;
			p2 = first;
		}
		
		temp_length = dist(p1, p2);
		if (temp_length < side_length) { side_length = temp_length; }

		iter = iter->next;
	}
	return side_length;
}

bool compare_paths(path_list_node* path1, path_list_node* path2, double length_interval, double angle_interval)
{
	//double length_interval = 1;	// in meters
	//double angle_interval = 0.0436; // in radians	 //2.5;	// in degrees
	// Compares path1 and path2 and determines if they are "identical"
	// Firstly compares the size of the path (i.e. the number of reflections)
	// If they are the same, compares path lengths
	// If they are within length_interval, compare direction of departure

	angle_interval = angle_interval*PI/180;		// Convert it into radians
	if(path1->num_refs != path2->num_refs)
	{
		// Paths have different number of reflections
		return false;
	}
	else
	{
		if(abs(path1->path_length - path2->path_length) < length_interval)
		{
			// Paths have the same number of reflections and the path lengths are within length_interval
			// Compare direction of departure
			if(abs(path1->tx_phi - path2->tx_phi) < angle_interval && 
			   abs(path1->tx_theta - path2->tx_theta) < angle_interval)
			{			
				// Paths are modeling the same GO path
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}
}

/*
bool compare_paths(path_list_node* path1, path_list_node* path2)
{
	// Compares path1 and path2 and determines if they are "identical"
	// Firstly compares the size of the path (i.e. the number of reflections)
	// If they are the same, compares path lengths
	// If they are the same, compare reflection points
	if(path1->num_refs != path2->num_refs)
	{
		// Paths have different number of reflections
		return false;
	}
	else
	{
		if(iszero(path1->path_length - path2->path_length))
		{
			// Paths have the same number of reflections and the same path lengths
			// Compare reflection points
			path_node* iter1 = path1->path;
			path_node* iter2 = path2->path;
			while(iter1 != NULL)
			{
				if(compare_points(iter1->p, iter2->p))
				{
					iter1 = iter1->next;
					iter2 = iter2->next;
				}
				else // reflection point is not the same
				{
					return false;
				}
			}
			// All reflection points are the same
			return true;
		}
		else
		{
			return false;
		}
	}
}
*/

bool compare_points(const point& p1, const point& p2)
{
	if(iszero(p1.x-p2.x) && iszero(p1.y-p2.y) && iszero(p1.z-p2.z))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void remove_duplicate_paths(path_list_node* start_path, double length_interval, double angle_interval)
{
	path_list_node* iter;
	path_list_node* prev_iter;
	while(start_path != NULL)
	{
		iter = start_path->next;
		prev_iter = start_path;
		while(iter != NULL)
		{
			if(compare_paths(start_path, iter, length_interval, angle_interval))
			{
				// Means that iter is pointing to a duplicate path
				// Remove it
				cout << "Removing path ... " << iter->path_id << endl;

				if(iter->next == NULL)
				{
					prev_iter->next = NULL;
					delete iter;
					iter = NULL;
				}
				else
				{
					prev_iter->next = iter->next;
					delete iter;
					iter = prev_iter->next;
				}
			}
			else
			{
				iter = iter->next;
				prev_iter = prev_iter->next;
			}
		}
		start_path = start_path->next;
	}
}

void compute_departure_angles(vect ray_dir, rect_cord coord, double &theta, double &phi)
{
	// The vertical angle theta is between the ray_dir vector and the az component
	// of the transmitter co-ordinate system
	// acos produces a result that is between 0 and pi
	double rd_ax, rd_ay, rd_az;
	rd_az = ray_dir*coord.az;
	theta = acos(rd_az/(mag(coord.az)*mag(ray_dir)));
	rd_ax = ray_dir*coord.ax;
	rd_ay = ray_dir*coord.ay;
	// atan2 produces a result in the range [-pi,pi]
	// when result is negative add 2pi, hence result will be in the range [0,2pi]
	phi = atan2(rd_ay, rd_ax);
	if (phi < 0)
	{
		phi = phi + TWOPI;
	}
}

bool intersection_ray_xcylinder (xcylinder &xc, ray &r, double &t_value, point &ap)
{
	// The function adjusts the intersections point by projecting it onto the circular-arc	
	double a, b, c, yo, zo, yd, zd, yc, zc;		

	yd = r.dir.y;
	zd = r.dir.z;
	yo = r.orig.y;
	zo = r.orig.z;
	yc = xc.yc;
	zc = xc.zc;

	// Ensure normal is outward (points from facetized to curve)
	// compute dot product between the possibe normal (yd, zd) and 
	// (yo-yc, zo-zc) the vector pointing from the center to the origin of ray r
	// The origin of ray r should be the original intersection point (that is to be adjusted)
	double dotp = yd*(yo-yc) + zd*(zo-zc);
	if (dotp < 0) // The normal (yd, zd) is pointing inward, flip it
	{
		yd = -yd;
		zd = -zd;
		r.dir = -1*r.dir;
	}

	a = yd*yd + zd*zd;
	b = 2*(yd*(yo-yc) + zd*(zo-zc));
	c = (yo-yc)*(yo-yc) + (zo-zc)*(zo-zc) - xc.r*xc.r;
	double disc = b*b - 4*a*c;

	// In this case the discriminant must be greater than 0
	// Since there must be an intersection in this case
	if (disc >= 0)	// Should asssert this
	{
		// Then there is an intersection
		t_value = -(b - sqrt(disc))/(2*a);	// Assert that t_value should be positive
		if (t_value < 0)
		{
		//	cout << "Warning: t_value was negative, but not to worry ... perhaps!";
		}
		
		ap = r.orig + t_value*r.dir;
		return true;		
	}
	return false;
}

bool intersection_ray_cylinder (cylinder &cyl, ray &r, double &t1, double &t2, point &p1, point &p2)
{
	// ADD COMMENTS
	vect r0 = point2vect(r.orig);
	vect alpha = r.dir - (r.dir*cyl.u)*cyl.u;
	vect beta = r0 - point2vect(cyl.p0) - ((r0-cyl.p0)*cyl.u)*cyl.u;

	double a = alpha*alpha;
	double b = 2*alpha*beta;
	double c = beta*beta - cyl.r*cyl.r;
	double disc;

	if (iszero(a))
	{
		if (iszero(b))
		{
			// There is no intersection
			return false;
		}
		else
		{
			// There is only one root
			t1 = -c/b;
		}
	}
	else
	{
		disc = b*b - 4*a*c;

		// In this case the discriminant must be greater than 0
		// Since there must be an intersection in this case
		if (disc >= 0)	// Should asssert this
		{
			// Then there is an intersection
			t1 = -(b - sqrt(disc))/(2*a);
			t2 = -(b + sqrt(disc))/(2*a);
						
			p1 = r.orig + t1*r.dir;
			p2 = r.orig + t2*r.dir;
			return true;		
		}
	}

	return true;
}

void adjust_ref_point(point tx, point rx, point ref, surface surf)
{


}

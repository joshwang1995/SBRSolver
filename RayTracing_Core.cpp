#include "RayTracing_Core.hpp"

int main(int argc, char** argv)
{
	float x,y,z;
	float theta,phi;
	int subdivisions=0;
	bool useRayTubes = false;
	
	cout << "Enter x-coordinate of origin: ";
	cin >> x;
	cout << "Enter y-coordinate of origin: ";
	cin >> y;
	cout << "Enter z-coordinate of origin: ";
	cin >> z;
	
	cout << "Enter number of subdivisions for icosahedron: ";
	cin >> subdivisions;
	
	cout << "Generate ray tubes [1=true,0=false]: ";
	cin >> useRayTubes;
	
	Vect3f origin = {x,y,z};
	
	//Construct regular icosahedron, hence subdivision = 0
	auto icosahedron = ConstructIcosahedron(0);	
	vector<Vect3f> ray_vects = GetRayVectOnIcosahedron(origin, icosahedron.first, icosahedron.second);
	vector<pair<double,double>> ray_launch_angles = GetLaunchAngle(ray_vects);
	vector<vector<int>> ray_tubes = ConstructRayTubes(icosahedron.second, icosahedron.first.size());
	
	//Find all rays from regular icosahedron that intersects the line/plane
	
	if (useRayTubes)
	{
		print_rays(ray_launch_angles, ray_tubes);
	}
	else
	{
		print_rays(ray_launch_angles);
	}
}

//*************************************************************************************************************
// Exposed Functions
vector<Vect3f> GetRayVectOnIcosahedron(Vect3f originPoint, vector<Vect3f>& vertices, vector<Vect3u>& faces)
{
	vector<Vect3f> outbound_ray_vect;
	
	for (int i = 0; i < faces.size(); i++)
	{
		//Compute the cente of the triangle
        double center_x = (vertices[faces[i].v0].x + vertices[faces[i].v1].x + vertices[faces[i].v2].x) / 3.0;
        double center_y = (vertices[faces[i].v0].y + vertices[faces[i].v1].y + vertices[faces[i].v2].y) / 3.0;
        double center_z = (vertices[faces[i].v0].z + vertices[faces[i].v1].z + vertices[faces[i].v2].z) / 3.0;
		Vect3f center_Triangle = make_Vect3f(center_x, center_y, center_z);
		outbound_ray_vect.push_back(normalize(center_Triangle - originPoint));
	}
	
	return outbound_ray_vect;
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

vector<pair<double,double>> GetLaunchAngle(vector<Vect3f>& ray_vect)
{
	vector<pair<double,double>> result;
	for(Vect3f ray: ray_vect)
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
SpherePoint RectPoint2Sphere(Vect3f n)
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
void print_rays(vector<pair<double,double>>& ray_dir)
{
	for(auto i: ray_dir)
	{
		cout << "(Theta,Phi): (" << i.first << "," << i.second << ")\n";		
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
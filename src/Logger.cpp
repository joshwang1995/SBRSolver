// Header file used to store constants
#ifndef LOGGER
#define LOGGER

#include "Vect_Utility.hpp"
#include "Geometry.hpp"
#include <vector>
#include <fstream>

/* Logging Functions */

void csv_writeln(std::ofstream& file, Vect3d point_int, Cvect3d E)
{
    file << point_int.x << "," << point_int.y << "," << point_int.z << ",";
    file << real(E.x) << "," << imag(E.x) << ",";
    file << real(E.y) << "," << imag(E.y) << ",";
    file << real(E.z) << "," << imag(E.z) << "\n";
}

void log_rays(const std::vector<Ray>& rays)
{
    std::ofstream log_file;
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
#endif
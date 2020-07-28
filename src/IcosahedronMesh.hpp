#pragma once
#include <vector>
#include <map>
#include "Vect_Utility.hpp"

std::pair<std::vector<Vect3f>, std::vector<Vect3u>> ConstructIcosahedron(int subdivisions);
std::vector<Vect3u> Subdivide(std::vector<Vect3f>& vertices, std::vector<Vect3u> triangles);

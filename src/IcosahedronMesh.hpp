#ifndef ICOSAHEDRON
#define ICOSAHEDRON

#include <vector>
#include <map>
#include "Vect_Utility.hpp"

std::pair<std::vector<Vect3f>, std::vector<Vect3u>> ConstructIcosahedron(int subdivisions);

#endif
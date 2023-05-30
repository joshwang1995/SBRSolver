#pragma once

#include <vector>
#include <map>
#include "common/VecMatDef.h"

std::pair<std::vector<Vec3>, std::vector<Idx3>> ConstructIcosahedron(int subdivisions);
std::vector<Vec3>* GenerateRaysOnIcosahedron(int tessellation);
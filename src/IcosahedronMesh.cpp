#include "IcosahedronMesh.hpp"

int vertex_for_edge(std::map<std::pair<int, int>, int>& lookup, std::vector<Vect3f>& vertices, int first, int second)
{
	std::pair<int, int> key(first, second);
	
	if (key.first > key.second)
		std::swap(key.first, key.second);
	
	auto inserted = lookup.insert({ key, static_cast<int>(vertices.size()) });
	
	if (inserted.second)
	{
		auto& edge0 = vertices[first];
		auto& edge1 = vertices[second];
		auto point = normalize(edge0 + edge1);
		vertices.push_back(point);
	}

	return inserted.first->second;
}

std::vector<Vect3u> Subdivide(std::vector<Vect3f>& vertices, std::vector<Vect3u> triangles)
{
	std::map<std::pair<int32_t, int32_t>, int32_t> lookup;
	std::vector<Vect3u> result;

	for (auto&& each : triangles)
	{
		Vect3u mid;	
		mid.v0 = vertex_for_edge(lookup, vertices, each.v0, each.v1);
		mid.v1 = vertex_for_edge(lookup, vertices, each.v1, each.v2);
		mid.v2 = vertex_for_edge(lookup, vertices, each.v2, each.v0);
		result.push_back({ each.v0, mid.v0, mid.v2 });
		result.push_back({ each.v1, mid.v1, mid.v0 });
		result.push_back({ each.v2, mid.v2, mid.v1 });
		result.push_back({ mid.v0, mid.v1, mid.v2 });
	}

	return result;
}

std::pair<std::vector<Vect3f>, std::vector<Vect3u>> ConstructIcosahedron(int subdivisions)
{
	// X,Z,N are coordinates of a icosahedron with edge-length = 2
	// To change radius of the icosahedron
	// multiply all coordinates by r/(2*sin(2*pi/5))
	// Source: https://math.stackexchange.com/questions/441327/coordinates-of-icosahedron-vertices-with-variable-radius
	
	const double X = .525731112119133606f;
	const double Z = .850650808352039932f;
	const double N = 0.f;
	
	std::vector<Vect3f> initVertices =
	{
		{-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
		{N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
		{Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
	};

	std::vector<Vect3u> iniTriangles =
	{
		{0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
		{8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
		{7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
		{6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
	};

	std::vector<Vect3f> vertices = initVertices;
	std::vector<Vect3u> triangles = iniTriangles;

	for (int i = 0; i < subdivisions; ++i)
	{
		triangles = Subdivide(vertices, triangles);
	}

	return{ vertices, triangles };
}
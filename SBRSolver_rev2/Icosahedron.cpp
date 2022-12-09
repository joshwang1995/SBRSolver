#include "icosahedron.h"

int vertex_for_edge(std::map<std::pair<int, int>, int>& lookup, std::vector<Vec3>& vertices, int first, int second)
{
	std::pair<int, int> key(first, second);
	
	if (key.first > key.second)
		std::swap(key.first, key.second);
	
	auto inserted = lookup.insert({ key, static_cast<int>(vertices.size()) });
	
	if (inserted.second)
	{
		auto& edge0 = vertices[first];
		auto& edge1 = vertices[second];
		auto point = (edge0+edge1).normalized();
		vertices.push_back(point);
	}

	return inserted.first->second;
}

std::vector<Idx3> Subdivide(std::vector<Vec3>& vertices, std::vector<Idx3> triangles)
{
	std::map<std::pair<int32_t, int32_t>, int32_t> lookup;
	std::vector<Idx3> result;

	for (auto&& each : triangles)
	{
		Idx3 mid;
		mid(0) = vertex_for_edge(lookup, vertices, each(0), each(1));
		mid(1) = vertex_for_edge(lookup, vertices, each(1), each(2));
		mid(2) = vertex_for_edge(lookup, vertices, each(2), each(0));
		result.push_back({ each(0), mid(0), mid(2)});
		result.push_back({ each(1), mid(1), mid(0)});
		result.push_back({ each(2), mid(2), mid(1) });
		result.push_back({ mid(0), mid(1), mid(2) });
	}

	return result;
}

std::pair<std::vector<Vec3>, std::vector<Idx3>> ConstructIcosahedron(int subdivisions)
{
	// X,Z,N are coordinates of a icosahedron with edge-length = 2
	// To change radius of the icosahedron
	// multiply all coordinates by r/(2*sin(2*pi/5))
	// Source: https://math.stackexchange.com/questions/441327/coordinates-of-icosahedron-vertices-with-variable-radius
	// subdivisions = 0 means the original icosahedron is returned
	const double X = .525731112119133606;
	const double Z = .850650808352039932;
	const double N = 0.0;
	
	std::vector<Vec3> initVertices =
	{
		{-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
		{N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
		{Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
	};

	std::vector<Idx3> iniTriangles =
	{
		{0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
		{8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
		{7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
		{6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
	};

	std::vector<Vec3> vertices = initVertices;
	std::vector<Idx3> triangles = iniTriangles;

	for (int i = 0; i < subdivisions; ++i)
	{
		triangles = Subdivide(vertices, triangles);
	}

	return{ vertices, triangles };
}

std::vector<Vec3>* GenerateRaysOnIcosahedron(int tessellation, Vec3 origin)
{
	auto icosahedron = ConstructIcosahedron(tessellation);
	std::vector<Vec3> vertices = icosahedron.first;
	std::vector<Idx3> faces = icosahedron.second;

	std::vector<Vec3>* outbound_rays = new std::vector<Vec3>;


	// Launch from the face center
	/*
	for (int i = 0; i < faces.size(); i++)
	{
		//Compute the cente of the triangle
		double center_x = (vertices[faces[i](0)](0) + vertices[faces[i](1)](0) + vertices[faces[i](2)](0)) / 3.0f;
		double center_y = (vertices[faces[i](0)](1) + vertices[faces[i](1)](1) + vertices[faces[i](2)](1)) / 3.0f;
		double center_z = (vertices[faces[i](0)](2) + vertices[faces[i](1)](2) + vertices[faces[i](2)](2)) / 3.0f;
		Vec3 direction = (Vec3(center_x, center_y, center_z)).normalized();
		outbound_rays->push_back(direction);
	}
	*/
	
	
	// Launch from the vertices
	for (int i = 0; i < vertices.size(); i++)
	{
		outbound_rays->push_back(vertices[i].normalized());
	}
	
	return outbound_rays;
}
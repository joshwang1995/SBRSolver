#include "Preprocessor.h"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <iterator>
#include <map>
#include <sstream>
#include <utility>

namespace
{
    // Converts a whole token to a finite number; returns false for anything
    // else (so a textual header line is skipped instead of throwing).
    bool ParseNumber(const std::string& token, double& value)
    {
        if (token.empty()) return false;
        char* end = nullptr;
        value = std::strtod(token.c_str(), &end);
        return end == token.c_str() + token.size() && std::isfinite(value);
    }

    // Parses every token of `line` into `values`; returns false if any token
    // is not a number.
    bool ParseNumericLine(const std::vector<std::string>& tokens, std::vector<double>& values)
    {
        values.clear();
        for (const std::string& token : tokens)
        {
            double v;
            if (!ParseNumber(token, v)) return false;
            values.push_back(v);
        }
        return true;
    }
}

bool Preprocessor::ReadPatternFile(std::string fileName, GainMap& output)
{
    std::string line;
    std::ifstream patternFile(fileName);

    if (patternFile.is_open())
    {
        int lineNo = 0;
        while (getline(patternFile, line))
        {
            lineNo++;
            std::stringstream ss(line);
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> vstrings(begin, end);
            
            if (vstrings.size() != 6) continue; // header or blank line

            std::vector<double> values;
            if (!ParseNumericLine(vstrings, values))
            {
                std::cerr << "Warning: " << fileName << ":" << lineNo << ": skipping non-numeric line" << std::endl;
                continue;
            }
            ThetaPhiAngle angle{ float(values[0]), float(values[1]) };
            GainInfo gain{};
            for (int i = 2; i < 6; i++)
            {
                gain.data[i - 2] = float(values[i]);
            }
            output.insert({ angle,gain });
        }
        patternFile.close();
        return true;
    }
    else
    {
        return false;
    }
    return false;
}

bool Preprocessor::ReadLocationFile(std::string fileName, VecVec3& output)
{
    std::string line;
    std::ifstream locationFile(fileName);

    if (locationFile.is_open())
    {
        int lineNo = 0;
        while (getline(locationFile, line))
        {
            lineNo++;
            std::stringstream ss(line);
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> vstrings(begin, end);

            if (vstrings.size() != 3) continue; // header or blank line

            std::vector<double> values;
            if (!ParseNumericLine(vstrings, values))
            {
                std::cerr << "Warning: " << fileName << ":" << lineNo << ": skipping non-numeric line" << std::endl;
                continue;
            }
            output.push_back(Vec3(values[0], values[1], values[2]));
        }
        locationFile.close();
        return true;
    }
    else
    {
        return false;
    }
    return false;
}

bool Preprocessor::StlToGeometry
(
    const std::vector<StlSource>& sources,
    const std::map<int, int>& surfaceMaterials,
    std::vector<Triangle*>& output, 
    bool saveEdges, 
    bool saveFaces, 
    std::string dataDir
)
{
    // Merged arrays over all files, in stl_reader's layout: coords holds
    // x y z per vertex, tris holds three vertex indices per triangle.
    std::vector<float> coords, normals;
    std::vector<unsigned int> tris;
    std::vector<int> triangleMaterial;
    std::map<std::pair<int, int>, std::vector<int>> edgeMap;

    // Triangles appended by this call are released again if it fails.
    const size_t firstNew = output.size();
    auto fail = [&output, firstNew]()
    {
        for (size_t i = firstNew; i < output.size(); i++) delete output[i];
        output.resize(firstNew);
        return false;
    };

    try
    {
        for (const StlSource& source : sources)
        {
            std::vector<float> fileCoords, fileNormals;
            std::vector<unsigned int> fileTris, solids;
            stl_reader::ReadStlFile(source.fileName.c_str(), fileCoords, fileNormals, fileTris, solids);
            if (fileTris.empty())
            {
                throw std::runtime_error("No triangles in " + source.fileName);
            }
            // stl_reader keeps one normal per facet; a facet without endfacet
            // (truncated file) leaves more normals than triangles. Merging the
            // arrays index by index would then hand the wrong normal to every
            // later triangle, in this file and in every file after it.
            if (fileNormals.size() != fileTris.size())
            {
                throw std::runtime_error(source.fileName + ": " + std::to_string(fileNormals.size() / 3)
                    + " facet normals but " + std::to_string(fileTris.size() / 3)
                    + " triangles; the file is truncated or malformed");
            }

            // Vertex indices of this file are rebased onto the merged array.
            // (Vertices shared between files stay distinct, so coplanar
            // surfaces never span two files.)
            const unsigned int vertexOffset = static_cast<unsigned int>(coords.size() / 3);
            coords.insert(coords.end(), fileCoords.begin(), fileCoords.end());
            normals.insert(normals.end(), fileNormals.begin(), fileNormals.end());
            for (unsigned int vertex : fileTris)
            {
                tris.push_back(vertex + vertexOffset);
            }
            const size_t fileTriangles = fileTris.size() / 3;
            triangleMaterial.insert(triangleMaterial.end(), fileTriangles, source.materialId);
            std::cout << "\tLoaded " << fileTriangles << " triangles from " << source.fileName
                << " (material " << source.materialId << ")" << std::endl;
        }

        const size_t numTris = tris.size() / 3;
        std::vector<std::vector<int>> adjacencyList(numTris);
        for (size_t itri = 0; itri < numTris; ++itri)
        {
            size_t v1_idx = 3 * tris[3 * itri];
            size_t v2_idx = 3 * tris[3 * itri + 1];
            size_t v3_idx = 3 * tris[3 * itri + 2];
            
            Vec3 v1{ coords[v1_idx], coords[v1_idx + 1], coords[v1_idx + 2] };
            Vec3 v2{ coords[v2_idx], coords[v2_idx + 1], coords[v2_idx + 2] };
            Vec3 v3{ coords[v3_idx], coords[v3_idx + 1], coords[v3_idx + 2] };
            Vec3 norm{ normals[3 * itri], normals[3 * itri + 1], normals[3 * itri + 2] };

            Triangle* t = new Triangle(v1, v2, v3, norm, int(itri));
            t->materialId = triangleMaterial[itri];
            output.push_back(t);
            InsertEdgeIntoMap(v1_idx, v2_idx, v3_idx, itri, output, edgeMap, adjacencyList);
        }

        BfsCoplanarSurface(adjacencyList, output);

        std::string surfaceError;
        const bool surfacesOk = ApplySurfaceMaterials(surfaceMaterials, output, surfaceError);

        // Written even when an override is rejected, so the error below can
        // point at a Faces.vtk that describes the geometry just loaded (its
        // MaterialID then shows the per-file materials).
        if (saveEdges) SaveEdgesAsVtk(dataDir + "Edges.vtk", edgeMap, coords);
        if (saveFaces) SaveFacesAsVtk(dataDir + "Faces.vtk", output, coords, tris);

        if (!surfacesOk)
        {
            std::cerr << surfaceError
                << (saveFaces ? " (written to " + dataDir + "Faces.vtk)" : " (set save_faces = true to write it)")
                << std::endl;
            return fail();
        }
        return true;
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return fail();
    }
}

bool Preprocessor::ApplySurfaceMaterials
(
    const std::map<int, int>& surfaceMaterials,
    std::vector<Triangle*>& triangles,
    std::string& error
)
{
    error.clear();
    if (surfaceMaterials.empty()) return true;

    int maxSurfaceId = -1;
    for (const Triangle* t : triangles)
    {
        maxSurfaceId = std::max(maxSurfaceId, t->coplanarId);
    }

    std::string missing;
    for (const auto& kv : surfaceMaterials)
    {
        if (kv.first < 0 || kv.first > maxSurfaceId)
        {
            missing += (missing.empty() ? "" : ", ") + std::to_string(kv.first);
        }
    }
    if (!missing.empty())
    {
        error = "[surface_materials] refers to surface id(s) " + missing + " but the scene has "
            + std::to_string(maxSurfaceId + 1) + " coplanar surfaces (ids 0 to " + std::to_string(maxSurfaceId)
            + "); see the CoplanarID field of Faces.vtk";
        return false;
    }

    for (Triangle* t : triangles)
    {
        std::map<int, int>::const_iterator it = surfaceMaterials.find(t->coplanarId);
        if (it != surfaceMaterials.end())
        {
            t->materialId = it->second;
        }
    }
    return true;
}

void Preprocessor::InsertEdgeIntoMap
(
    const size_t& v1, 
    const size_t& v2, 
    const size_t& v3,
    const size_t& triId,
    const std::vector<Triangle*>& triangles,
    std::map<std::pair<int, int>, std::vector<int>>& edgeMap,
    std::vector<std::vector<int>>& adjacencyList
)
{
    // Always insert the two vertex with the smaller value first. This would avoid situations where <v1,v2> and <v2,v1>
    // are treated as different edges.
    std::pair<int, int> edge1 = v1 < v2 ? std::make_pair(int(v1), int(v2)) : std::make_pair(int(v2), int(v1));
    std::pair<int, int> edge2 = v2 < v3 ? std::make_pair(int(v2), int(v3)) : std::make_pair(int(v3), int(v2));
    std::pair<int, int> edge3 = v1 < v3 ? std::make_pair(int(v1), int(v3)) : std::make_pair(int(v3), int(v1));

    // Find if edge1 is in map, append the triangle ID if it is in the map. Insert the edge and the triangle ID if it is not
    std::map<std::pair<int, int>, std::vector<int>>::iterator it = edgeMap.find(edge1);
    if (it != edgeMap.end())
    {
        for (const int& i : it->second)
        {
            if ((triangles[i]->norm == triangles[triId]->norm) || (triangles[i]->norm == -triangles[triId]->norm))
            {
                adjacencyList[i].push_back(triId);
                adjacencyList[triId].push_back(i);
            }
        }
        it->second.push_back(int(triId));
    }
    else
    {
        edgeMap[edge1] = std::vector<int>();
        edgeMap[edge1].push_back(int(triId));
    }

    it = edgeMap.find(edge2);
    if (it != edgeMap.end())
    {
        for (const int& i : it->second)
        {
            if ((triangles[i]->norm == triangles[triId]->norm) || (triangles[i]->norm == -triangles[triId]->norm))
            {
                adjacencyList[i].push_back(triId);
                adjacencyList[triId].push_back(i);
            }
        }
        it->second.push_back(int(triId));
    }
    else
    {
        edgeMap[edge2] = std::vector<int>();
        edgeMap[edge2].push_back(int(triId));
    }

    it = edgeMap.find(edge3);
    if (it != edgeMap.end())
    {
        for (const int& i : it->second)
        {
            if ((triangles[i]->norm == triangles[triId]->norm) || (triangles[i]->norm == -triangles[triId]->norm))
            {
                adjacencyList[i].push_back(triId);
                adjacencyList[triId].push_back(i);
            }
        }
        it->second.push_back(int(triId));
    }
    else
    {
        edgeMap[edge3] = std::vector<int>();
        edgeMap[edge3].push_back(int(triId));
    }
}

void Preprocessor::BfsCoplanarSurface(const std::vector<std::vector<int>>& adj, const std::vector<Triangle*>& triangles)
{
    int V = adj.size();
    std::vector<bool> visited(V, false);
    int coplanarId = 0;

    for (int u = 0; u < V; u++)
    {
        if (visited[u] == false)
        {
            BFSUtil(u, adj, visited, triangles, coplanarId);
            coplanarId++;
        }
    }
}

void Preprocessor::BFSUtil
(
    int u,
    const std::vector<std::vector<int>>& adj,
    std::vector<bool>& visited,
    const std::vector<Triangle*>& triangles,
    int coplanarId
)
{
    // Create a queue for BFS
    std::queue<int> q;

    // Mark the current node as visited and enqueue it
    visited[u] = true;
    q.push(u);

    // 'i' will be used to get all adjacent vertices 4
    // of a vertex list<int>::iterator i;

    while (!q.empty())
    {
        // Dequeue a vertex from queue and print it
        u = q.front();
        triangles[u]->coplanarId = coplanarId;
        q.pop();

        // Get all adjacent vertices of the dequeued
        // vertex s. If an adjacent has not been visited,
        // then mark it visited and enqueue it
        for (int i = 0; i != adj[u].size(); ++i)
        {
            if (!visited[adj[u][i]])
            {
                visited[adj[u][i]] = true;
                q.push(adj[u][i]);
            }
        }
    }
}


bool Preprocessor::ReadMaterialsFile(std::string fileName, Materials& output)
{
    return false;
}

void Preprocessor::GenerateRxPlane(double xMin, double yMin, double xMax, double yMax, double height, double resolution, VecVec3& output)
{
    int nx = int((xMax - xMin) / resolution);
    int ny = int((yMax - yMin) / resolution);
    
    for (int i = 0; i <= nx; i++)
    {
        double x = xMin + i * resolution;
        for (int j = 0; j <= ny; j++)
        {
            double y = yMin + j * resolution;
            output.emplace_back(Vec3(x, y, height));
        }
    }
}

void Preprocessor::GenerateRXLine(Vec3 start, Vec3 end, double resolution, VecVec3& output)
{
    Vec3 direction = end - start;
    int numPoints = static_cast<int>(std::floor(direction.norm() / resolution));
    direction.normalize();
    
    // Guard in case where start and end is the same point
    if (numPoints <= 0) assert("Number of RX points will be smaller or equal to 0, check start and end points!");
    
    for (int i = 0; i < numPoints; i++)
    {
        Vec3 point = start + i * resolution * direction;
        output.push_back(point);
    }
}

void Preprocessor::GenerateAntennaPattern(std::string, double resolution)
{
}

bool Preprocessor::SaveEdgesAsVtk
(
    std::string fileName, 
    const std::map<std::pair<int, int>, std::vector<int>>& edgeMap, 
    const std::vector<float>& coords
)
{
    using namespace std;

    cout << "[Entering] Preprocessor::SaveEdgesAsVtk ..." << endl;

    ofstream ofs;
    ofs.open(fileName);

    // write header
    ofs << "# vtk DataFile Version 2.0" << endl;
    ofs << "3D model of edge geometry" << endl;
    ofs << "ASCII" << endl;
    ofs << "DATASET POLYDATA" << endl;
    ofs << endl;

    ofs << "POINTS " << int(coords.size()/3) << " float" << endl;
    for (int i = 0; i < coords.size(); i++)
    {
        if (i > 0 && i % 3 == 0) ofs << endl;
        ofs << " " << coords[i];
    }

    map<pair<int, int>, vector<int>>::const_iterator it = edgeMap.begin();
    std::vector<pair<int, int>> edgePoints;
    std::vector<int> numConnections;
    while (it != edgeMap.end())
    {
        edgePoints.push_back(it->first);
        numConnections.push_back(int(it->second.size()));
        it++;
    }

    ofs << "LINES " << edgePoints.size() << " " << 3 * edgePoints.size() << endl;
    for (int i = 0; i < edgePoints.size(); i++)
    {
        ofs << " 2 " << edgePoints[i].first / 3 << " " << edgePoints[i].second / 3 << endl;
    }

    ofs << "CELL_DATA " << edgePoints.size() << endl;
    ofs << "FIELD FieldData 2" << endl;
    ofs << "EdgeID 1 " << edgePoints.size() << " int" << endl;
    for (int i = 0; i < edgePoints.size(); i++)
    {
        ofs << " " << i << endl;
    }

    ofs << "Connection 1 " << numConnections.size() << " int" << endl;
    for (int i = 0; i < numConnections.size(); i++)
    {
        ofs << " " << numConnections[i] << endl;
    }

    ofs.close();
    cout << "\tSaved " << edgePoints.size() << " edges into " << fileName << endl;
    cout << "[Leaving] Preprocessor::SaveEdgesAsVtk" << endl;

    return true;
}


bool Preprocessor::SaveFacesAsVtk
(
    std::string fileName, 
    const std::vector<Triangle*>& triangles,
    const std::vector<float>& coords,
    const std::vector<unsigned int> tris
)
{
    using namespace std;

    cout << "[Entering] Preprocessor::SaveFacesAsVtk ..." << endl;

    ofstream ofs;
    ofs.open(fileName);

    // write header
    ofs << "# vtk DataFile Version 2.0" << endl;
    ofs << "3D model of face geometry" << endl;
    ofs << "ASCII" << endl;
    ofs << "DATASET POLYDATA" << endl;
    ofs << endl;

    ofs << "POINTS " << int(coords.size()/3) << " float" << endl;
    for (int i = 0; i < coords.size(); i++)
    {
        if (i > 0 && i % 3 == 0) ofs << endl;
        ofs << " " << coords[i];
    }
    ofs << endl;

    int numTris = int(tris.size() / 3);
    ofs << "POLYGONS " << numTris << " " << 4 * numTris << endl;
    for (int i = 0; i < numTris; i++)
    {
        ofs << " 3 "
            << tris[3 * i] << " "
            << tris[3 * i + 1] << " "
            << tris[3 * i + 2] << endl;
    }

    ofs << "CELL_DATA " << numTris << endl;
    ofs << "FIELD FieldData 3" << endl;
    ofs << "FaceID 1 " << triangles.size() << " int" << endl;
    for (int i = 0; i < triangles.size(); i++)
    {
        ofs << " " << triangles[i]->triangleId << endl;
    }

    ofs << "CoplanarID 1 " << triangles.size() << " int" << endl;
    for (int i = 0; i < triangles.size(); i++)
    {
        ofs << " " << triangles[i]->coplanarId << endl;
    }

    ofs << "MaterialID 1 " << triangles.size() << " int" << endl;
    for (int i = 0; i < triangles.size(); i++)
    {
        ofs << " " << triangles[i]->materialId << endl;
    }

    ofs.close();
    cout << "\tSaved " << numTris << " faces into " << fileName << endl;
    cout << "[Leaving] Preprocessor::SaveFacesAsVtk" << endl;

    return true;
}
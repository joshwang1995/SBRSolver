#include "Preprocessor.h"

bool Preprocessor::ReadPatternFile(std::string fileName, GainMap& output)
{
    std::string line;
    std::ifstream patternFile(fileName);

    if (patternFile.is_open())
    {
        while (getline(patternFile, line))
        {
            std::stringstream ss(line);
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> vstrings(begin, end);
            
            if (vstrings.size() != 6) continue; // Need to add function to check if vector has only numeric values

            ThetaPhiAngle angle{ std::stof(vstrings[0]), std::stof(vstrings[1]) };
            GainInfo gain{};
            for (int i = 2; i < 6; i++)
            {
                gain.data[i - 2] = std::stof(vstrings[i]);
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
        while (getline(locationFile, line))
        {
            std::stringstream ss(line);
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> vstrings(begin, end);
            Vec3 loc{};

            if (vstrings.size() != 3) continue;
            std::transform(vstrings.begin(), vstrings.end(), loc.begin(), [](const std::string& val)
            {
                return std::stod(val);
            });
            output.push_back(loc);
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
    std::string fileName, 
    std::vector<Triangle*>& output, 
    bool saveEdges, 
    bool saveFaces, 
    std::string dataDir
)
{
    std::vector<float> coords, normals;
    std::vector<unsigned int> tris, solids;
    std::map<std::pair<int, int>, std::vector<int>> edgeMap;

    try
    {
        stl_reader::ReadStlFile(fileName.c_str(), coords, normals, tris, solids);
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
            output.push_back(t);
            InsertEdgeIntoMap(v1_idx, v2_idx, v3_idx, itri, output, edgeMap, adjacencyList);
        }

        Timer timer;
        timer.start();
        BfsCoplanarSurface(adjacencyList, output);
        double endTime = timer.getTime();

        std::cout << "Total Time for BFS algorithm: " << endTime << std::endl;
        std::cout << "\tNumber of Triangles: " << numTris << std::endl;
        
        int numEdges = 0;
        for (std::vector<int> i : adjacencyList)
        {
            numEdges += i.size();
        }

        std::cout << "\tNumber of Edges: " << numEdges << std::endl;

        if (saveEdges) SaveEdgesAsVtk(dataDir + "Edges.vtk", edgeMap, coords);
        // if (saveFaces) SaveFacesAsVtk(dataDir + "Faces.vtk", );
        return true;
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }
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
    std::pair<int, int> edge1 = v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2,v1);
    std::pair<int, int> edge2 = v2 < v3 ? std::make_pair(v2, v3) : std::make_pair(v3, v2);
    std::pair<int, int> edge3 = v1 < v3 ? std::make_pair(v1, v3) : std::make_pair(v3, v1);

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

/*
bool Preprocessor::SaveFacesAsVtk(std::string fileName, const std::vector<float>& )
{
}
*/
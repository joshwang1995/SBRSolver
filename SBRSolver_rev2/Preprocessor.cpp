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

bool Preprocessor::StlToGeometry(std::string fileName, std::vector<Triangle*>& output)
{
    std::vector<float> coords, normals;
    std::vector<unsigned int> tris, solids;

    try
    {
        stl_reader::ReadStlFile(fileName.c_str(), coords, normals, tris, solids);
        const size_t numTris = tris.size() / 3;
        // std::cout << "Total Number of faces: " << (int)numTris << std::endl;
        for (size_t itri = 0; itri < numTris; ++itri)
        {
            size_t v1_idx = 3 * tris[3 * itri];
            size_t v2_idx = 3 * tris[3 * itri + 1];
            size_t v3_idx = 3 * tris[3 * itri + 2];

            Vec3 v1{ coords[v1_idx], coords[v1_idx + 1], coords[v1_idx + 2] };
            Vec3 v2{ coords[v2_idx], coords[v2_idx + 1], coords[v2_idx + 2] };
            Vec3 v3{ coords[v3_idx], coords[v3_idx + 1], coords[v3_idx + 2] };
            Vec3 norm{ normals[3 * itri], normals[3 * itri + 1], normals[3 * itri + 2] };

            Triangle* t = new Triangle(v1, v2, v3, norm);
            output.push_back(t);
            //std::cout << "Triangle " << itri << ": ";
            //std::cout << v1.format(CommaInitFmt) << "," << v2.format(CommaInitFmt) << "," << v3.format(CommaInitFmt) << std::endl;
        }
        return true;
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }
}

bool Preprocessor::ReadMaterialsFile(std::string fileName, Materials& output)
{
    return false;
}

void Preprocessor::GenerateRxPlane(double xMin, double yMin, double xMax, double yMax, double height, double resolution)
{
}

void Preprocessor::GenerateAntennaPattern(std::string, double resolution)
{
}
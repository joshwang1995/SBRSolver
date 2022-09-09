#pragma once

#include "common/VecMatDef.h"
#include "BVH/BVHNode.h"
#include "common/DataStructures.h"
#include "common/Constants.h"

class Triangle
{
public:
    Triangle();
    ~Triangle();

    Triangle(Vec3 a, Vec3 b, Vec3 c, Vec3 d, Mat23 e, Vec3 f)
        :
        v1(a), v2(b), v3(c), norm(d), bbox(e), center(f) {}

    Triangle(Vec3 a, Vec3 b, Vec3 c); //v1, v2, v3
    Triangle(Vec3 a, Vec3 b, Vec3 c, Vec3 d); //v1, v2, v3, normal

    Vec3 v1; // Vertex 1
    Vec3 v2; // Vertex 2
    Vec3 v3; // Vertex 3
    Vec3 norm; // Triangle Normal
    Mat23 bbox; // Bounding box
    Vec3 center; // Bounding box center

    Mat23 findBbox() const;
    Vec3 findCenter() const;
    Vec3 findNormal() const;
    bool RayIntersects
    (
        const Vec3& rayOrg,
        const Vec3& rayDir,
        HitInfo& info
    );

    bool RayIntersects
    (
        const Vec3& rayOrg,
        const Vec3& rayDir,
        double& distance,
        Vec3& pointIntersect
    );
};
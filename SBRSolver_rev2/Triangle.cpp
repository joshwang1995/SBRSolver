#include "Triangle.h"

Triangle::Triangle()
{
    BVHNode<Triangle>::EmptyAABB(bbox);
    v1 = Vec3();
    v2 = Vec3();
    v3 = Vec3();
    norm = Vec3();
    center = Vec3();
}

Triangle::~Triangle()
{
}

Triangle::Triangle(Vec3 a, Vec3 b, Vec3 c, int d)
{
    v1 = a; v2 = b; v3 = c;
    triangleId = d;
    norm = findNormal();
    bbox = findBbox();
    center = findCenter();
}

Triangle::Triangle(Vec3 a, Vec3 b, Vec3 c, Vec3 d, int e)
{
    v1 = a; v2 = b; v3 = c; norm = d;
    triangleId = e;
    bbox = findBbox();
    center = findCenter();
}

Mat23 Triangle::findBbox() const
{
    Mat23 bbox;
    bbox(0, 0) = std::min({ v1(0), v2(0), v3(0) });
    bbox(0, 1) = std::min({ v1(1), v2(1), v3(1) });
    bbox(0, 2) = std::min({ v1(2), v2(2), v3(2) });
    bbox(1, 0) = std::max({ v1(0), v2(0), v3(0) });
    bbox(1, 1) = std::max({ v1(1), v2(1), v3(1) });
    bbox(1, 2) = std::max({ v1(2), v2(2), v3(2) });
    return bbox;
}

Vec3 Triangle::findCenter() const
{
    Vec3 center;
    center(0) = (std::min({ v1(0), v2(0), v3(0) }) + std::max({ v1(0), v2(0), v3(0) })) / 2.0;
    center(1) = (std::min({ v1(1), v2(1), v3(1) }) + std::max({ v1(1), v2(1), v3(1) })) / 2.0;
    center(2) = (std::min({ v1(2), v2(2), v3(2) }) + std::max({ v1(2), v2(2), v3(2) })) / 2.0;
    return center;
}

Vec3 Triangle::findNormal() const
{
    Vec3 edge1 = v2 - v1;
    Vec3 edge2 = v3 - v1;
    Vec3 normal = (edge1.cross(edge2)).normalized();
    return normal;
}

bool Triangle::RayIntersects
(
    const Vec3& rayOrg,
    const Vec3& rayDir,
    HitInfo& info
)
{
    Vec3 edge1, edge2, h, s, q;
    double a, f, u, v;

    edge1 = v2 - v1;
    edge2 = v3 - v1;
    h = rayDir.cross(edge2);
    a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.

    f = 1.0 / a;

    s = rayOrg - v1;
    u = f * s.dot(h);
    if (u < 0.0 || u > 1.0)
        return false;

    q = s.cross(edge1);
    v = f * rayDir.dot(q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot(q);
    if (t > EPSILON) // ray intersection
    {
        info.pointIntersect = rayOrg + rayDir * t;
        info.distance = (info.pointIntersect - rayOrg).norm();
        info.normal = *&norm;
        info.materialId = materialId;
        info.triangelId = triangleId;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
};
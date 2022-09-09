//
// (c) 2020 Takahiro Hashimoto
//
#pragma once
#include "../common/VecMatDef.h"
#include "../common/Constants.h"

template <typename T>
class BVHNode
{
public:
	BVHNode();
	~BVHNode();

	Mat23 bbox; // bounding box
	Idx2 children; // child nodes

	std::vector<T*> objects; // objects in leaf box

	VecVec3 CornerPoints() const;
	bool IntersectsAABB(const Vec3& rayOrigin, const Vec3& rayVector);

	static void EmptyAABB(Mat23& bbox);

};

template<typename T>
BVHNode<T>::BVHNode()
{
	EmptyAABB(bbox);
}

template<typename T>
BVHNode<T>::~BVHNode()
{
}

template<typename T>
VecVec3 BVHNode<T>::CornerPoints() const
{
	VecVec3 varr(8);
	varr[0] = { bbox(0, 0), bbox(0, 1), bbox(0, 2) };
	varr[1] = { bbox(1, 0), bbox(0, 1), bbox(0, 2) };
	varr[2] = { bbox(0, 0), bbox(1, 1), bbox(0, 2) };
	varr[3] = { bbox(1, 0), bbox(1, 1), bbox(0, 2) };
	varr[4] = { bbox(0, 0), bbox(0, 1), bbox(1, 2) };
	varr[5] = { bbox(1, 0), bbox(0, 1), bbox(1, 2) };
	varr[6] = { bbox(0, 0), bbox(1, 1), bbox(1, 2) };
	varr[7] = { bbox(1, 0), bbox(1, 1), bbox(1, 2) };

	return varr;
}

template<typename T>
void BVHNode<T>::EmptyAABB(Mat23& bbox)
{
	for (int axis = 0; axis < 3; axis++)
	{
		bbox(0, axis) = std::numeric_limits<float>::max();
		bbox(1, axis) = -std::numeric_limits<float>::max();
	}
}

template<typename T>
bool BVHNode<T>::IntersectsAABB
(
	const Vec3& rayOrigin,
	const Vec3& rayVector
)
{
	double t_max = INF;
	double t_min = -INF;

	for (int i = 0; i < 3; i++)
	{
		double t1 = (bbox(0, i) - rayOrigin(i)) / rayVector(i);
		double t2 = (bbox(1, i) - rayOrigin(i)) / rayVector(i);
		double t_near = std::min(t1, t2);
		double t_far = std::max(t1, t2);
		t_max = std::min(t_max, t_far);
		t_min = std::max(t_min, t_near);

		if (t_min > t_max) return false;
	}
	return true;
}
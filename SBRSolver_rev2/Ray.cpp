#include "Ray.h"

extern std::atomic<int> RayGlobalId = 0;

PathTreeNode* newPathTreeNode
(
	Vec3 SourcePoint, 
	Vec3 TargetPoint, 
	int PenetrationMaterialId, 
	int ReflectionMaterialId, 
	double PathLength, 
	double AngleFromSurfaceNormal
)
{
	PathTreeNode* result = new PathTreeNode;
	result->ray.id = RayGlobalId++;
	result->ray.sourcePoint = SourcePoint;
	result->ray.targetPoint = TargetPoint;
	result->ray.penetrationMaterialId = PenetrationMaterialId;
	result->ray.reflectionMaterialId = ReflectionMaterialId;
	result->ray.angleFromSurfaceNormal = AngleFromSurfaceNormal;
	result->ray.pathLength = PathLength;
	result->childDirect = nullptr;
	result->childReflect = nullptr;
	return result;
}

TreeNode::TreeNode()
{
}

TreeNode::~TreeNode()
{
	child.clear();
}

Paths::Paths()
{

}

Paths::Paths(TreeNode* rootNode)
{
	if (rootNode == nullptr || rootNode->child.empty())
	{
		return;
	}

	std::queue<std::pair<std::pair<std::vector<Ray>, bool>, TreeNode*>> q;
	std::vector<Ray> newVect;
	newVect.push_back(rootNode->ray);
	q.push(std::make_pair(std::make_pair(newVect, false), rootNode));
	while (!q.empty())
	{
		auto p = q.front();
		q.pop();
		if (p.second->child.size() > 0)
		{
			for (auto i = p.second->child.begin(); i != p.second->child.end(); ++i)
			{
				std::vector<Ray> newVect(p.first.first);
				newVect.push_back(i->ray);
				q.push(std::make_pair(std::make_pair(newVect, p.first.second || i->ray.reflectionMaterialId >= 0), &*i));
			}
		}
		else
		{
			if (p.first.second)
			{
				rayPaths.push_back(p.first.first);
			}
		}
	}
}

Paths::Paths(PathTreeNode* rootNode)
{
	// This function traverses through the pathtreenode and pushes each node's 
	// ray into newVect. 
	if (rootNode == nullptr || (rootNode->childDirect == nullptr && rootNode->childReflect == nullptr))
	{
		return;
	}

	std::queue<std::pair<std::pair<std::vector<Ray>, bool>, PathTreeNode*>> q;
	std::vector<Ray> newVect;
	newVect.push_back(rootNode->ray);
	q.push(std::make_pair(std::make_pair(newVect, false), rootNode));
	while (!q.empty())
	{
		auto p = q.front();
		q.pop();
		if (p.second->childDirect != nullptr || p.second->childReflect != nullptr)
		{
			if (p.second->childDirect != nullptr)
			{
				bool terminate = p.first.second || p.second->childDirect->ray.reflectionMaterialId >= 0 || p.second->childDirect->ray.penetrationMaterialId >= 0;
				std::vector<Ray> newVect(p.first.first);
				newVect.push_back(p.second->childDirect->ray);
				q.push(std::make_pair(std::make_pair(newVect, terminate), p.second->childDirect));
			}

			if (p.second->childReflect != nullptr)
			{
				bool terminate = p.first.second || p.second->childReflect->ray.reflectionMaterialId >= 0 || p.second->childReflect->ray.penetrationMaterialId >= 0;
				std::vector<Ray> newVect(p.first.first);
				newVect.push_back(p.second->childReflect->ray);
				q.push(std::make_pair(std::make_pair(newVect, terminate), p.second->childReflect));
			}
		}
		else
		{
			if (p.first.second)
			{
				rayPaths.push_back(p.first.first);
			}
		}
	}
}


Paths::~Paths()
{
	for (auto i = rayPaths.begin(); i != rayPaths.end(); ++i)
	{
		i->clear();
	}
	rayPaths.clear();
}
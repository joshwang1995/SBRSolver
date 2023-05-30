#include "Ray.h"

extern std::atomic<int> RayGlobalId = 0;

PathTreeNode* newPathTreeNode
(
	Vec3 sourcePoint, 
	Vec3 targetPoint, 
	int penetrationMaterialId, 
	int reflectionMaterialId,
	int surfaceId,
	int coplanarId,
	double pathLength
)
{
	PathTreeNode* result = new PathTreeNode;
	result->ray.id = RayGlobalId++;
	result->ray.sourcePoint = sourcePoint;
	result->ray.targetPoint = targetPoint;
	result->ray.penetrationMaterialId = penetrationMaterialId;
	result->ray.reflectionMaterialId = reflectionMaterialId;
	result->ray.hitSurfaceID = surfaceId;
	result->ray.hitCoplanarId = coplanarId;
	result->ray.pathLength = pathLength;
	result->childTransmit = nullptr;
	result->childReflect = nullptr;
	return result;
}

PathTreeNode* DeleteChildNodes(PathTreeNode* node)
{
	if (node != nullptr)
	{
		if (node->childTransmit != nullptr)
		{
			node->childTransmit = DeleteChildNodes(node->childTransmit);
		}

		if (node->childReflect != nullptr)
		{
			node->childReflect = DeleteChildNodes(node->childReflect);
		}
		delete node;
	}
	return nullptr;
}

void CloneTree(PathTreeNode* orgTree, PathTreeNode* cloneTree)
{
	if (orgTree != nullptr) 
	{
		//Direct ray
		PathTreeNode* newTransmitNode = CloneNode(orgTree->childTransmit);
		cloneTree->childTransmit = newTransmitNode;
		CloneTree(orgTree->childTransmit, cloneTree->childTransmit);

		//Reflect ray
		PathTreeNode* newReflectNode = CloneNode(orgTree->childReflect);
		cloneTree->childReflect = newReflectNode;
		CloneTree(orgTree->childReflect, cloneTree->childReflect);
	}
}

PathTreeNode* CloneNode(PathTreeNode* node)
{
	if (node != nullptr)
	{
		PathTreeNode* newNode = new PathTreeNode;
		newNode->childTransmit = nullptr;
		newNode->childReflect = nullptr;
		newNode->ray.id = node->ray.id;
		newNode->ray.sourcePoint = node->ray.sourcePoint;
		newNode->ray.targetPoint = node->ray.targetPoint;
		newNode->ray.reflectionMaterialId = node->ray.reflectionMaterialId;
		newNode->ray.penetrationMaterialId = node->ray.penetrationMaterialId;
		newNode->ray.pathLength = node->ray.pathLength;
		newNode->ray.hitSurfaceID = node->ray.hitSurfaceID;
		newNode->ray.hitCoplanarId = node->ray.hitCoplanarId;
		newNode->ray.captured = node->ray.captured;
		return newNode;
	}
	return nullptr;
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
#pragma omp parallel
	{
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
#pragma omp critical
					rayPaths.push_back(p.first.first);
				}
			}
		}
	}
}

Paths::Paths(PathTreeNode* rootNode)
{
	// This function traverses through the pathtreenode and pushes each node's 
	// ray into newVect. 
	// if (rootNode == nullptr || (rootNode->childTransmit == nullptr && rootNode->childReflect == nullptr))
	if (rootNode == nullptr)
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
		if (p.second->childTransmit != nullptr || p.second->childReflect != nullptr)
		{
			if (p.second->childTransmit != nullptr)
			{
				//bool terminate = p.first.second || p.second->childTransmit->ray.reflectionMaterialId >= 0 || p.second->childTransmit->ray.penetrationMaterialId >= 0;
				bool terminate = p.first.second || p.second->childTransmit->ray.captured;
				std::vector<Ray> newVect(p.first.first);
				newVect.push_back(p.second->childTransmit->ray);
				q.push(std::make_pair(std::make_pair(newVect, terminate), p.second->childTransmit));
			}

			if (p.second->childReflect != nullptr)
			{
				//bool terminate = p.first.second || p.second->childReflect->ray.reflectionMaterialId >= 0 || p.second->childReflect->ray.penetrationMaterialId >= 0;
				bool terminate = p.first.second || p.second->childReflect->ray.captured;
				std::vector<Ray> newVect(p.first.first);
				newVect.push_back(p.second->childReflect->ray);
				q.push(std::make_pair(std::make_pair(newVect, terminate), p.second->childReflect));
			}
		}
		else
		{
			if (p.first.second || p.second->ray.captured)
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
#include "Ray.h"

std::atomic<int> RayGlobalId{0};

PathTreeNode* newPathTreeNode
(
	Vec3 sourcePoint, 
	Vec3 targetPoint, 
	int penetrationMaterialId, 
	int reflectionMaterialId,
	int surfaceId,
	int coplanarId,
	double pathLength, 
	double angleFromSurfaceNormal
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
	result->ray.angleFromSurfaceNormal = angleFromSurfaceNormal;
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
		newNode->ray.angleFromSurfaceNormal = node->ray.angleFromSurfaceNormal;
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
	// Collects one candidate path per captured node: the rays from the root
	// down to that node, the last one ending at the receiver. A captured node
	// keeps its children, so a ray that passes a receiver, bounces and passes
	// it again yields both the shorter and the longer path. The last ray's
	// surface ids are cleared because it now ends at the receiver, not at the
	// surface it hit; the tree itself is left untouched for the other paths.
	if (rootNode == nullptr)
	{
		return;
	}

	std::queue<std::pair<std::vector<Ray>, PathTreeNode*>> q;
	q.push(std::make_pair(std::vector<Ray>{ rootNode->ray }, rootNode));
	while (!q.empty())
	{
		std::pair<std::vector<Ray>, PathTreeNode*> p = q.front();
		q.pop();
		if (p.second->ray.captured)
		{
			std::vector<Ray> captured(p.first);
			captured.back().hitSurfaceID = -1;
			captured.back().hitCoplanarId = -1;
			rayPaths.push_back(captured);
		}
		for (PathTreeNode* child : { p.second->childTransmit, p.second->childReflect })
		{
			if (child == nullptr) continue;
			std::vector<Ray> extended(p.first);
			extended.push_back(child->ray);
			q.push(std::make_pair(extended, child));
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
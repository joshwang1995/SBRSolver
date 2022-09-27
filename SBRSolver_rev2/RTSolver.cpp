#include "RTSolver.h"

RTSolver::RTSolver()
{
	Cleanup();
}

RTSolver::~RTSolver()
{
	Cleanup();
}

int RTSolver::Init(MaterialProperties materialProperties[], int materialPropertiesCount, BVH<Triangle>& triangles)
{
	Cleanup();

	_materialProperties = new MaterialProperties[materialPropertiesCount];
	std::copy(materialProperties, materialProperties + materialPropertiesCount, _materialProperties);
	_bvh = &triangles;
	_isSceneInitialized = true;
	return 0;
}

void RTSolver::Cleanup()
{
	if (_rayPaths != nullptr)
	{
		DeleteRayPaths();
		_rayPaths = nullptr;
	}

	if (_materialProperties != nullptr)
	{
		delete[] _materialProperties;
		_materialProperties = nullptr;
	}

	if (_shootRayList != nullptr)
	{
		delete _shootRayList;
		_shootRayList = nullptr;
	}
	_isSceneInitialized = false;
}

int RTSolver::ExecuteRayTracing
(
	Vec3 sourcePoint,
	int maxReflectionCount,
	int maxTransmissionCount,
	VecVec3 receivers,
	int txTesslation
)
{
	_maxReflectionCount = maxReflectionCount;
	_maxTransmissionCount = maxTransmissionCount;
	_txTesslation = txTesslation + 1; // Since user can enter 0 as txTesslation
	RayGlobalId = 0;
	_shootRayList = GenerateRaysOnIcosahedron(txTesslation, sourcePoint);
	_receiverCount = int(receivers.size());
	_pathsCount = int(_shootRayList->size());
	InitRayPaths();

	//Timer timer;
	//timer.start();
	// #pragma omp parallel for
	for (int i = 0; i < _pathsCount; i++)
	{
		PathTreeNode rootNode;
		rootNode.ray.id = RayGlobalId++;
		rootNode.childTransmit = nullptr;
		rootNode.childReflect = nullptr;
		RayLaunch(&rootNode, sourcePoint, _shootRayList->at(i), -1, -1, 0, 0, 0, true, FLT_MAX);
		for (int j = 0; j < _receiverCount; j++)
		{
			PathTreeNode* cloneRoot = CloneNode(&rootNode);
			CloneTree(&rootNode, cloneRoot);
			RayCapture(cloneRoot, receivers[j], 0.0);
			_rayPaths[i][j] = Paths(cloneRoot);
			cloneRoot = DeleteChildNodes(cloneRoot);
			delete cloneRoot;
		}
	}
	//std::cout << "\tTotal Time in for loop -> " << timer.getTime() << std::endl;
	return _pathsCount * _receiverCount;
}

void RTSolver::RayLaunch
(
	PathTreeNode* rayTreeNode, 
	Vec3& sourcePoint, 
	Vec3& directionPoint, 
	int transMaterialID, 
	int refMaterialID, 
	int reflectionCnt, 
	int transmissionCnt, 
	double totalPathLength,
	bool isRoot, 
	double lastAnglefromN
)
{
	if (reflectionCnt >= _maxReflectionCount || transmissionCnt >= _maxTransmissionCount)
	{
		return;
	}

	HitInfo hitResult;
	bool hasHit = _bvh->RayIntersects(0, sourcePoint, directionPoint, hitResult);

	if (hasHit)
	{
		double anglefromN = AngleBetween(hitResult.normal, directionPoint);
		PathTreeNode* newRayTreeNode = rayTreeNode;
		if (!isRoot)
		{
			PathTreeNode* newNode = newPathTreeNode
			(
				sourcePoint,
				hitResult.pointIntersect,
				transMaterialID,
				refMaterialID,
				hitResult.triangelId,
				totalPathLength + hitResult.distance,
				lastAnglefromN
			);

			if (refMaterialID >= 0)
				rayTreeNode->childReflect = newNode;
			else
				rayTreeNode->childTransmit = newNode;
			newRayTreeNode = newNode;
		}
		else //it is the root
		{
			rayTreeNode->ray.sourcePoint = sourcePoint;
			rayTreeNode->ray.targetPoint = hitResult.pointIntersect;
			rayTreeNode->ray.penetrationMaterialId = transMaterialID;
			rayTreeNode->ray.reflectionMaterialId = refMaterialID;
			rayTreeNode->ray.hitSurfaceID = hitResult.triangelId;
			rayTreeNode->ray.pathLength = hitResult.distance;
			rayTreeNode->ray.angleFromSurfaceNormal = lastAnglefromN;
		}

		if (reflectionCnt < _maxReflectionCount)
		{
			Vec3 newDirectionPoint = Reflect(directionPoint, hitResult.normal).normalized();
			// Rebounce ray
			RayLaunch
			(
				newRayTreeNode, 
				hitResult.pointIntersect, 
				newDirectionPoint, 
				-1, 
				hitResult.materialId, 
				reflectionCnt + 1,
				transmissionCnt,
				totalPathLength + hitResult.distance,
				false, 
				anglefromN
			);
		}

		if (transmissionCnt < _maxTransmissionCount)
		{
			// Forward ray
			RayLaunch
			(
				newRayTreeNode,
				hitResult.pointIntersect,
				directionPoint,
				hitResult.materialId,
				-1,
				reflectionCnt,
				transmissionCnt + 1,
				totalPathLength + hitResult.distance,
				false,
				anglefromN
			);
		}

	}
	else
	{
		if (reflectionCnt > 0 || transmissionCnt > 0)
		{
			PathTreeNode* newNode = newPathTreeNode
			(
				sourcePoint, 
				sourcePoint + (directionPoint * 1E2), 
				transMaterialID,
				refMaterialID,
				-1,
				FLT_MAX,
				lastAnglefromN
			);
			if (refMaterialID >= 0)
				rayTreeNode->childReflect = newNode;
			else if (transMaterialID >= 0)
				rayTreeNode->childTransmit = newNode;
		}

		if (isRoot)
		{
			// Direct Ray
			rayTreeNode->ray.sourcePoint = sourcePoint;
			rayTreeNode->ray.targetPoint = sourcePoint + (directionPoint * 1E2);
			rayTreeNode->ray.reflectionMaterialId = -1;
			rayTreeNode->ray.penetrationMaterialId = -1;
			rayTreeNode->ray.hitSurfaceID = -1;
			rayTreeNode->ray.pathLength = FLT_MAX;
		}
	}
}

void RTSolver::RayCapture(PathTreeNode* rayTreeNode, const Vec3& receiver, double totalPathLength)
{
	Vec3 capturePoint;
	bool captured = HitReceptionSphere
	(
		rayTreeNode->ray.sourcePoint, 
		rayTreeNode->ray.targetPoint, 
		receiver,
		totalPathLength, 
		_txTesslation, 
		capturePoint
	);

	if (rayTreeNode->childTransmit == nullptr && rayTreeNode->childReflect == nullptr && !captured)
	{
		// Delete current node
		// rayTreeNode = DeleteChildNodes(rayTreeNode);
		return;
	}

	if (captured)
	{
		// Check if vector from capture point to RX is occluded
		// Problem here: 
		//	The ray intersection test checks if a ray intersects with the scene. A ray can be infinitely long in one direction.
		//	Obviously, an infintely long ray can intersect with the surrounding geometry. Therefore, we should instead check if the line segment
		//	from the point of capture to the receiver is occulded. To do so, we need to implement an line segment intersection test.
		HitInfo hitResult;
		bool hasHit = _bvh->RayIntersects(0, capturePoint, (receiver - capturePoint).normalized(), hitResult);
		if (hasHit)
		{
			double segmentLength = (receiver - capturePoint).norm();
			hasHit = hitResult.distance < segmentLength ? true : false;
		}

		if (!hasHit)
		{
			rayTreeNode->ray.captured = true;
			rayTreeNode->ray.targetPoint = capturePoint;
			rayTreeNode->ray.pathLength = (capturePoint - rayTreeNode->ray.sourcePoint).norm();
			rayTreeNode->childTransmit = DeleteChildNodes(rayTreeNode->childTransmit);
			rayTreeNode->childReflect = DeleteChildNodes(rayTreeNode->childReflect);
		}
	}
	else
	{
		if (rayTreeNode->childTransmit != nullptr)
		{
			RayCapture(rayTreeNode->childTransmit, receiver, totalPathLength + rayTreeNode->ray.pathLength);
		}

		if (rayTreeNode->childReflect != nullptr)
		{
			RayCapture(rayTreeNode->childReflect, receiver, totalPathLength + rayTreeNode->ray.pathLength);
		}
	}
	
}

void RTSolver::RemoveDuplicatePath(int receiverId)
{
	int level = 0;
	for (int i = 0; i < _pathsCount; i++)
	{
		// _rayPaths[i][receiverId].rayPaths.at(level);
	}
}

void RTSolver::InitRayPaths()
{
	_rayPaths = new Paths* [_pathsCount];
	for (int i = 0; i < _pathsCount; i++)
	{
		_rayPaths[i] = new Paths[_receiverCount];
	}
}

void RTSolver::DeleteRayPaths()
{ 
	for (int i = 0; i < _pathsCount; i++)
	{
		delete[] _rayPaths[i];
	}
	delete[] _rayPaths;
}

void RTSolver::CmdLineDebug()
{
	using namespace std;
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");

	for (int i = 0; i < _pathsCount; i++)
	{
		cout << "Launch ID: " << i << endl;
		for (int j = 0; j < _receiverCount; j++)
		{
			cout << "\tReceiver ID: " << j << endl;
			for (int k = 0; k < _rayPaths[i][j].rayPaths.size(); k++)
			{
				cout << "\t\tPath " << k + 1 << endl;
				for (int l = 0; l < _rayPaths[i][j].rayPaths[k].size(); l++)
				{
					cout << "\t\t\t Ray ID: " << _rayPaths[i][j].rayPaths[k][l].id << endl;
					cout << "\t\t\t\t Source Point: " << _rayPaths[i][j].rayPaths[k][l].sourcePoint.format(CommaInitFmt) << endl;
					cout << "\t\t\t\t Target Point: " << _rayPaths[i][j].rayPaths[k][l].targetPoint.format(CommaInitFmt) << endl;
					cout << "\t\t\t\t Reflection Material: " << _rayPaths[i][j].rayPaths[k][l].reflectionMaterialId << endl;
					cout << "\t\t\t\t Penetration Material: " << _rayPaths[i][j].rayPaths[k][l].penetrationMaterialId << endl;
					cout << "\t\t\t\t Path Length: " << _rayPaths[i][j].rayPaths[k][l].pathLength << endl;
					cout << "\t\t\t\t Hit Triangle ID: " << _rayPaths[i][j].rayPaths[k][l].hitSurfaceID << endl;
					cout << "\t\t\t\t Angle From Surface Normal: " << _rayPaths[i][j].rayPaths[k][l].angleFromSurfaceNormal << endl;
					cout << "\t\t\t\t Captured: " << _rayPaths[i][j].rayPaths[k][l].captured << endl;
				}
			}
		}
	}
}

bool RTSolver::SavePathsAsVtk(std::string fname)
{
	using namespace std;

	cout << "[Entering] RTSolver::SavePathsAsVtk ..." << endl;

	ofstream ofs;
	ofs.open(fname);

	// write header
	ofs << "# vtk DataFile Version 2.0" << endl;
	ofs << "Ray Path Visualization" << endl;
	ofs << "ASCII" << endl;
	ofs << "DATASET POLYDATA" << endl;
	ofs << endl;

	vector<Vec3> vertex;
	vector<vector<int>> lineIdx;
	vector<int> launchIds;
	vector<int> receiverIds;
	int totalLines = 0;

	for (int i = 0; i < _pathsCount; i++)
	{
		for (int j = 0; j < _receiverCount; j++)
		{
			PathsToVector(_rayPaths[i][j].rayPaths, vertex, lineIdx, launchIds, receiverIds, i, j, totalLines);
		}
	}

	ofs << "POINTS " << vertex.size() << " float" << endl;
	for (const Vec3& v : vertex)
	{
		ofs << " "
			<< v.x() << " "
			<< v.y() << " "
			<< v.z() << endl;
	}

	ofs << "LINES " << lineIdx.size() << " " << totalLines << endl;
	for(int i = 0; i < lineIdx.size(); i++)
	{
		ofs << " " << lineIdx[i].size() << " ";
		for (int j = 0; j < lineIdx[i].size(); j++)
		{
			ofs << lineIdx[i][j] << " ";
		}
		ofs << endl;
	}

	ofs << "CELL_DATA " << launchIds.size() << endl;
	ofs << "FIELD FieldData 2" << endl;
	ofs << "LaunchID 1 " << launchIds.size() << " int" << endl;
	for (int i = 0; i < launchIds.size(); i++)
	{
		ofs << " " << launchIds[i] << endl;
	}

	ofs << "ReceiverID 1 " << receiverIds.size() << " int" << endl;
	for (int i = 0; i < receiverIds.size(); i++)
	{
		ofs << " " << receiverIds[i] << endl;
	}

	ofs.close();
	cout << "\tSaved " << lineIdx.size() << " ray paths into " << fname << endl;
	cout << "[Leaving] RTSolver::SavePathsAsVtk" << endl;

	return true;
}

bool RTSolver::SaveIcosahedronAsVtk(std::string fname, Vec3 rayOrg, int tessellation)
{
	using namespace std;

	cout << "[Entering] RTSolver::SaveIcosahedronAsVtk ..." << endl;

	ofstream ofs;
	ofs.open(fname);

	// write header
	ofs << "# vtk DataFile Version 2.0" << endl;
	ofs << "TX Icosahedron Visualization" << endl;
	ofs << "ASCII" << endl;
	ofs << "DATASET POLYDATA" << endl;
	ofs << endl;

	ofs << "POINTS " << _shootRayList->size() + 1 << " float" << endl;
	ofs << " " << rayOrg.x() << " " << rayOrg.y() << " " << rayOrg.z() << endl;
	for (const Vec3& v : *_shootRayList)
	{
		ofs << " "
			<< rayOrg.x() + 1.0 * v.x() << " "
			<< rayOrg.y() + 1.0 * v.y() << " "
			<< rayOrg.z() + 1.0 * v.z() << endl;
	}

	ofs << "LINES " << _shootRayList->size() << " " << 3 * _shootRayList->size() << endl;
	for (int i = 1; i < _shootRayList->size() + 1; i++)
	{
		ofs << " 2"
			<< " 0 "
			<< i << endl;
	}

	ofs << "CELL_DATA " << _shootRayList->size() << endl;
	ofs << "FIELD FieldData 1" << endl;
	ofs << "LaunchID 1 " << _shootRayList->size() << " int" << endl;
	for (int i = 0; i < _shootRayList->size(); i++)
	{
		ofs << " " << i << endl;
	}

	ofs.close();
	cout << "\tSaved " << _shootRayList->size() << " ray paths into " << fname << endl;
	cout << "[Leaving] RTSolver::SaveIcosahedronAsVtk" << endl;

	return true;
}

void PathsToVector
(
	const std::vector<std::vector<Ray>>& rayPaths,
	std::vector<Vec3>& vertices,
	std::vector<std::vector<int>>& lineIndex,
	std::vector<int>& launchIds,
	std::vector<int>& receiverIds,
	int launchId,
	int receiverId,
	int& totalLines
)
{
	for (int i = 0; i < rayPaths.size(); i++)
	{
		std::vector<int> lines;
		for (int j = 0; j < rayPaths[i].size(); j++)
		{
			if (j == 0)
			{
				vertices.push_back(rayPaths[i][j].sourcePoint);
				lines.push_back(int(vertices.size() - 1));
				vertices.push_back(rayPaths[i][j].targetPoint);
				lines.push_back(int(vertices.size() - 1));
			}
			else
			{
				vertices.push_back(rayPaths[i][j].targetPoint);
				lines.push_back(int(vertices.size() - 1));
			}
		}
		totalLines += int(lines.size()) + 1;
		lineIndex.push_back(lines);
		launchIds.push_back(launchId);
		receiverIds.push_back(receiverId);
	}
}

bool HitReceptionSphere
(
	const Vec3& sourcePoint,
	const Vec3& targetPoint,
	const Vec3& sphereCenter,
	double pathLength,
	int txTesslation,
	Vec3& capturePoint
)
{
	// Find the closest distance between the ray and the sphere center
	// Reference: https://www.geometrictools.com/Documentation/DistancePointLine.pdf
	Vec3 rayDir = targetPoint - sourcePoint; // Do not normalize, this is a line segment
	double t0 = rayDir.dot(sphereCenter - sourcePoint) / rayDir.squaredNorm();

	if (t0 < 0)
	{
		// Projection is negative, receiver is in the opposite direction of the ray direction
		return false;
	}
	else if (t0 >= 1)
	{
		t0 = 1; // Clamp t0 to 1 because the closest point is the target point
	}

	// Find the optimal reception sphere radius
	// Reference: 191115_Ray_Launching.pdf
	// Note: the unfolded path length is the path length from previous ray segments plus the distance
	//  the point of intersection and the receiver

	// Takahiro's method
	//double beta = 1.0 / (3.0 * txTesslation) * acos(-1.0 / sqrt(5.0));
	//double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	//double sphereRadius = unfoldedLength * tan(beta);

	// Rappaport's approximation
	// double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	// double alpha = (69.0 * PI / 180.0) / txTesslation;
	// double sphereRadius = alpha * unfoldedLength / sqrt(3.0);

	double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	double phi = (sqrt(5) + 1) / 2;
	double alpha = acos(phi / (phi * phi + 1)) / (pow(2, txTesslation - 1)) / sqrt(3);
	double sphereRadius = alpha * unfoldedLength;


	double distance = (sphereCenter - (sourcePoint + t0 * rayDir)).norm();
	if (distance <= sphereRadius)
	{
		// Distance from closest point on the ray to sphere center is smaller
		// than the sphere radius, therefore it is not captured
		capturePoint = sourcePoint + t0 * rayDir;
		return true;
	}
	return false;
}

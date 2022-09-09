#include "RTSolver.h"

RTSolver::RTSolver()
{
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
		delete[] _rayPaths;
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

int RTSolver::ExecuteRayTracing(Vec3 sourcePoint, int maxBounceCount, double maxPathLoss, int txTesslation)
{
	if (_rayPaths != nullptr)
	{
		delete[] _rayPaths;
		_rayPaths = nullptr;
	}
	_maxBounceCount = maxBounceCount;
	_maxPathLoss = maxPathLoss;
	RayGlobalId = 0;
	_shootRayList = GenerateRaysOnIcosahedron(txTesslation, sourcePoint);
	_rayPaths = new Paths[_shootRayList->size()];

	Timer timer;
	timer.start();
	#pragma omp parallel for
	for (int i = 0; i < _shootRayList->size(); i++)
	{
		PathTreeNode rootNode;
		rootNode.ray.id = RayGlobalId++;
		rootNode.childDirect = nullptr;
		rootNode.childReflect = nullptr;
		RayLaunch(&rootNode, sourcePoint, _shootRayList->at(i), -1, -1, 0, 0, 0, true, FLT_MAX);
		_rayPaths[i] = Paths(&rootNode);
	}
	std::cout << "\tTotal Time in for loop -> " << timer.getTime() << std::endl;
	_pathsCount = static_cast<int>(_shootRayList->size());
	return _pathsCount;
}


void RTSolver::RayLaunch(PathTreeNode* rayTreeNode, Vec3& sourcePoint, Vec3& directionPoint, int transmissionLossId, int reflectionMaterialId, int bounceCnt, int penetrationCnt, double totalLoss, bool isRoot, float lastAnglefromN)
{
	if (bounceCnt >= _maxBounceCount || penetrationCnt >= 2)
	{
		return;
	}

	if (isRoot)
	{
		directionPoint = directionPoint;
	}

	HitInfo hitResult;
	bool hasHit = _bvh->RayIntersects(0, sourcePoint, directionPoint, hitResult);

	if (hasHit)
	{
		float anglefromN = AngleBetween(hitResult.normal, directionPoint);
		PathTreeNode* newRayTreeNode = rayTreeNode;
		if (!isRoot)
		{
			PathTreeNode* newNode = newPathTreeNode(sourcePoint, hitResult.pointIntersect, transmissionLossId, reflectionMaterialId, lastAnglefromN);
			if (reflectionMaterialId >= 0)
				rayTreeNode->childReflect = newNode;
			else
				rayTreeNode->childDirect = newNode;
			newRayTreeNode = newNode;
		}
		else //it is the root
		{
			rayTreeNode->ray.sourcePoint = sourcePoint;
			rayTreeNode->ray.targetPoint = hitResult.pointIntersect;
			rayTreeNode->ray.penetrationMaterialId = transmissionLossId;
			rayTreeNode->ray.reflectionMaterialId = reflectionMaterialId;
			rayTreeNode->ray.angleFromSurfaceNormal = lastAnglefromN;
		}

		if ((bounceCnt < _maxBounceCount && penetrationCnt == 0) || (bounceCnt < _maxBounceCount && penetrationCnt > 0))
		{
			Vec3 inRay = directionPoint;
			Vec3 wallNormal = hitResult.normal;
			Vec3 newDirectionPoint = Reflect(inRay, wallNormal).normalized();
			// Rebounce ray
			RayLaunch(newRayTreeNode, hitResult.pointIntersect, newDirectionPoint, -1, 0, bounceCnt + 1, penetrationCnt, totalLoss, false, anglefromN);
		}

		// Forward ray
		//Vec3 pointIntersectOffset = hitResult.pointIntersect + (directionPoint * EPSILON);
		RayLaunch(newRayTreeNode, hitResult.pointIntersect, directionPoint, 0, -1, bounceCnt, penetrationCnt + 1, totalLoss, false, anglefromN);
	}
	else
	{
		if (bounceCnt > 0)
		{
			PathTreeNode* newNode = newPathTreeNode(sourcePoint, sourcePoint + (directionPoint * 1E2), transmissionLossId, reflectionMaterialId, lastAnglefromN);
			if (reflectionMaterialId >= 0)
				rayTreeNode->childReflect = newNode;
			else
				rayTreeNode->childDirect = newNode;
		}
	}
}

bool RTSolver::RayLaunch2(Vec3& rayOrig, Vec3& rayDir, double totalPathLen, int bounceCnt, int penetrationCnt)
{

	HitInfo hitResult;
	bool hasHit = _bvh->RayIntersects(0, rayOrig, rayDir, hitResult);

	// Either a direct ray or ray with no interation
	if (!hasHit)
	{

	}
	return true;
}

void RTSolver::DebugFunc()
{
	using namespace std;
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");

	for (int i = 0; i < _pathsCount; i++)
	{
		cout << "Path " << i << ":" << endl;
		for (int j = 0; j < _rayPaths[i].rayPaths.size(); j++)
		{
			cout << "\t Vector " << j << ":" << endl;
			for (int k = 0; k < _rayPaths[i].rayPaths[j].size(); k++)
			{
				Ray r = _rayPaths[i].rayPaths[j][k];
				cout << "\t\t Ray " << r.id << ": ";
				cout << r.sourcePoint.format(CommaInitFmt) << " -> " << r.targetPoint.format(CommaInitFmt) << endl;
			}
		}
		cout << endl;
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

	vector<Vec3> points;
	vector<Idx2> lineIdx;

	for (int i = 0; i < _pathsCount; i++)
	{
		for (int j = 0; j < _rayPaths[i].rayPaths.size(); j++)
		{
			for (int k = 0; k < _rayPaths[i].rayPaths[j].size(); k++)
			{
				points.push_back(_rayPaths[i].rayPaths[j][k].sourcePoint);
				points.push_back(_rayPaths[i].rayPaths[j][k].targetPoint);
				lineIdx.emplace_back(Idx2(points.size() - 2, points.size() - 1));
			}
		}
	}

	ofs << "POINTS " << points.size() << " float" << endl;
	for (const Vec3& v : points)
	{
		ofs << " "
			<< v.x() << " "
			<< v.y() << " "
			<< v.z() << endl;
	}

	ofs << "LINES " << lineIdx.size() << " " << 3 * lineIdx.size() << endl;
	for (const Idx2& l : lineIdx)
	{
		ofs << " 2 "
			<< l.x() << " "
			<< l.y() << endl;
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
			<< rayOrg.x() + 2.0 * v.x() << " "
			<< rayOrg.y() + 2.0 * v.y() << " "
			<< rayOrg.z() + 2.0 * v.z() << endl;
	}

	ofs << "LINES " << _shootRayList->size() << " " << 3 * _shootRayList->size() << endl;
	for (int i = 0; i < _shootRayList->size(); i++)
	{
		ofs << " 2"
			<< " 0 "
			<< i << endl;
	}

	ofs.close();
	cout << "\tSaved " << _shootRayList->size() << " ray paths into " << fname << endl;
	cout << "[Leaving] RTSolver::SaveIcosahedronAsVtk" << endl;

	return true;
}
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

	if (_efield != nullptr)
	{
		delete _efield;
		_efield = nullptr;
	}
	_isSceneInitialized = false;
}

int RTSolver::ExecuteRayTracing
(
	Vec3 sourcePoint,
	const Mat3& txCoordSys,
	int maxReflectionCount,
	int maxTransmissionCount,
	double frequency,
	double pt,
	bool useFresnelCoeff,
	VecVec3& receivers,
	std::vector<Triangle*>& triangleMesh,
	int txTesslation
)
{
	RayGlobalId = 0;

	_maxReflectionCount = maxReflectionCount;
	_maxTransmissionCount = maxTransmissionCount;
	_txTesslation = txTesslation + 1; // user can enter 0 as txTesslation
	_shootRayList = GenerateRaysOnIcosahedron(txTesslation, sourcePoint);
	_receivers = &receivers;
	_receiverCount = int(receivers.size());
	_pathsCount = int(_shootRayList->size());
	_efield = new std::vector<Vec3c>;

	InitRayPaths();

	Timer timer;
	timer.start();
//#pragma omp parallel for collapse(2)
	for (int i = 0; i < _pathsCount; i++)
	{
		PathTreeNode rootNode;
		rootNode.ray.id = RayGlobalId++;
		rootNode.childTransmit = nullptr;
		rootNode.childReflect = nullptr;
		RayLaunch(&rootNode, sourcePoint, _shootRayList->at(i), -1, -1, 0, 0, 0, true, INF);
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

	FieldCompute* fieldCore = new FieldCompute
	(
		_rayPaths, 
		triangleMesh, 
		receivers, 
		_pathsCount, 
		_materialProperties, 
		txCoordSys, 
		frequency, 
		pt,
		useFresnelCoeff
	);

	Vec3c field = { cdouble(0,0), cdouble(0,0), cdouble(0,0) };

//#pragma omp parallel for
	for (int k = 0; k < _receiverCount; k++)
	{
		ImagePathCorrection(k, sourcePoint, triangleMesh);
		RemoveDuplicatePath(k);
		_efield->push_back(fieldCore->FieldAtReceiver(k));
	}
	std::cout << "\tTotal Time in for loop -> " << timer.getTime() << std::endl;
	
	delete fieldCore;
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
	if (reflectionCnt > _maxReflectionCount && transmissionCnt > _maxTransmissionCount)
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
				hitResult.coplanarId,
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
			rayTreeNode->ray.hitCoplanarId = hitResult.coplanarId;
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
				-1,
				INF,
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
			rayTreeNode->ray.hitCoplanarId = -1;
			rayTreeNode->ray.pathLength = INF;
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
		return;
	}

	if (captured)
	{
		rayTreeNode->ray.captured = true;
		rayTreeNode->ray.hitSurfaceID = -1;
		rayTreeNode->ray.hitCoplanarId = -1;
		rayTreeNode->childTransmit = DeleteChildNodes(rayTreeNode->childTransmit);
		rayTreeNode->childReflect = DeleteChildNodes(rayTreeNode->childReflect);
	}

	if (rayTreeNode->childTransmit != nullptr)
	{
		RayCapture(rayTreeNode->childTransmit, receiver, totalPathLength + rayTreeNode->ray.pathLength);
	}

	if (rayTreeNode->childReflect != nullptr)
	{
		RayCapture(rayTreeNode->childReflect, receiver, totalPathLength + rayTreeNode->ray.pathLength);
	}

	
}

void RTSolver::ImagePathCorrection(int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh)
{
	for (int i = 0; i < _pathsCount; i++)
	{
		auto it = _rayPaths[i][receiverId].rayPaths.begin();
		bool pathCorrected = false;
		while (it != _rayPaths[i][receiverId].rayPaths.end())
		{
			if (it->size() == 1)
			{
				pathCorrected = DirectPathCorrection(*it, receiverId, txPoint);
			}
			else
			{
				pathCorrected = MultiPathCorrection(*it, receiverId, txPoint, triangleMesh);
			}
			
			if (!pathCorrected)
			{
				// Delete iterator and retreat the next item
				it = _rayPaths[i][receiverId].rayPaths.erase(it);
			}
			else
			{
				it++;
			}
		}
	}
	return;
}

bool RTSolver::MultiPathCorrection(std::vector<Ray>& multiPath, int receiverId, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh)
{
	std::vector<std::pair<int, Vec3>> imgSources = GetImageSources(multiPath, txPoint, triangleMesh);
	if(imgSources.size() != multiPath.size()) 
	{ 
		assert("The GetImageSources function returned a list that is smaller than the path size");
	}

	HitInfo result;
	Vec3 targetPoint = _receivers->at(receiverId);

	for (int i = (multiPath.size() - 1); i >= 0; i--)
	{
		int hitSurfId = imgSources.at(i).first;
		Vec3 imgSrc = imgSources.at(i).second;
		Vec3 imgRayDir = targetPoint - imgSrc;
		
		if (hitSurfId < 0)
		{
			multiPath.at(i).sourcePoint = imgSrc;
			multiPath.at(i).targetPoint = targetPoint;
			multiPath.at(i).pathLength = (targetPoint - result.pointIntersect).norm();
			multiPath.at(i).angleFromSurfaceNormal = INF;
		}
		else
		{
			bool hasHit = triangleMesh[hitSurfId]->RayIntersects(imgSrc, imgRayDir, result);

			if (hasHit)
			{
				multiPath.at(i).sourcePoint = result.pointIntersect;
				multiPath.at(i).targetPoint = targetPoint;
				multiPath.at(i).pathLength = (targetPoint - result.pointIntersect).norm();
				multiPath.at(i).angleFromSurfaceNormal = AngleBetween(triangleMesh[hitSurfId]->norm, targetPoint - result.pointIntersect);
				targetPoint = result.pointIntersect;
			}
			else
			{
				return false;
			}
		}
	}
	return true;
}

bool RTSolver::DirectPathCorrection(std::vector<Ray>& directPath, int receiverId, const Vec3& txPoint)
{
	if (directPath.size() != 1) { return false; }

	// Direct path correction
	HitInfo result;
	Vec3 directRayDirection = _receivers->at(receiverId) - txPoint;
	bool hasHit = _bvh->RayIntersects(0, txPoint, directRayDirection, result);
	if (!hasHit || directRayDirection.norm() < result.distance)
	{
		directPath.at(0).sourcePoint = txPoint;
		directPath.at(0).targetPoint = _receivers->at(receiverId);
		directPath.at(0).pathLength = (directPath.at(0).targetPoint - directPath.at(0).sourcePoint).norm();
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<std::pair<int, Vec3>> GetImageSources(const std::vector<Ray>& path, const Vec3& txPoint, const std::vector<Triangle*>& triangleMesh)
{
	Vec3 currentImgSrc = txPoint;
	std::vector<std::pair<int, Vec3>> result {std::make_pair(-1, txPoint) };
	
	int surfId = -1;
	for (int i = 0; i < path.size() - 1; i++)
	{
		surfId = path[i].hitSurfaceID;

		if(surfId < 0)
		{
			// Shouldn't happen. This implies a non-direct ray in the path doesn't hit a surface
			return result;
			assert("The ray path is not valid");
		}
		if (path[i+1].reflectionMaterialId >= 0) // The interaction type is stored in the next ray
		{
			// if the interaction is reflection, then find image source
			currentImgSrc = ImagePoint(*triangleMesh[surfId], currentImgSrc);
			result.push_back(std::make_pair(surfId, currentImgSrc));
		}
		else if (path[i + 1].penetrationMaterialId >= 0)
		{
			// if the interaction is transmission, then keep the old image source
			result.push_back(std::make_pair(surfId, currentImgSrc));
		}
		else
		{
			// Shouldn't happen
			// When both reflectionMaterialId and penetrationMaterialId is < 1, this implies the current ray
			// is fired directly from the TX. but we are indexing path[i+1], so the zero case is not included here
			return result;
			assert("The ray path is not valid");
		}
	}
	return result;
}

void RTSolver::RemoveDuplicatePath(int receiverId)
{
	std::map<std::vector<int>, int> surfaceIdMap;

	// std::unordered_set<std::vector<int>> hitIdSet;

	for (int i = 0; i < _pathsCount; i++)
	{
		auto it = _rayPaths[i][receiverId].rayPaths.begin();
		while (it != _rayPaths[i][receiverId].rayPaths.end())
		{
			// std::vector<int> key = GetHitSurfaceIds(*it);
			std::vector<int> key = GetHitCoplanarIds(*it);
			if (surfaceIdMap.count(key))
			{
				it = _rayPaths[i][receiverId].rayPaths.erase(it);
			}
			else
			{
				surfaceIdMap.insert(std::make_pair(key,0));
				it++;
			}
		}
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
	ofs << "FIELD FieldData 5" << endl;
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

	ofs << "Eabs 1 " << receiverIds.size() << " double" << endl;
	for (int i = 0; i < receiverIds.size(); i++)
	{
		double eabs = std::isfinite(_efield->at(receiverIds[i]).norm()) ? _efield->at(receiverIds[i]).norm() : 0.0;
		ofs << " " << eabs << endl;
	}

	ofs << "EdB 1 " << receiverIds.size() << " double" << endl;
	for (int i = 0; i < receiverIds.size(); i++)
	{
		double eabs = std::isfinite(_efield->at(receiverIds[i]).norm()) ? 20 * log10(_efield->at(receiverIds[i]).norm() * 1e6) : -1.0;
		ofs << " " << eabs << endl;
	}

	ofs << "Phase 1 " << receiverIds.size() << " double" << endl;
	for (int i = 0; i < receiverIds.size(); i++)
	{
		double totalPhase = std::arg(_efield->at(receiverIds[i]).x()) + std::arg(_efield->at(receiverIds[i]).y()) + std::arg(_efield->at(receiverIds[i]).z());
		totalPhase = std::isfinite(totalPhase) ? totalPhase : 0.0;
		if (std::isfinite(totalPhase))
		{
			totalPhase >= TWOPI ? (totalPhase - TWOPI) * (180.0 / PI) : totalPhase * (180.0 / PI);
		}
		else
		{
			totalPhase = 361.0;
		}
		ofs << " " << totalPhase << endl;
	}

	ofs.close();
	cout << "\tSaved " << lineIdx.size() << " ray paths into " << fname << endl;
	cout << "[Leaving] RTSolver::SavePathsAsVtk" << endl;

	return true;
}

bool RTSolver::SaveReceiversAsVtk(std::string fname)
{
	using namespace std;

	cout << "[Entering] RTSolver::SaveReceiversAsVtk ..." << endl;

	ofstream ofs;
	ofs.open(fname);

	// write header
	ofs << "# vtk DataFile Version 2.0" << endl;
	ofs << "Ray Path Visualization" << endl;
	ofs << "ASCII" << endl;
	ofs << "DATASET POLYDATA" << endl;
	ofs << endl;

	ofs << "POINTS " << _receiverCount << " float" << endl;
	for(int i = 0; i < _receiverCount; i++)
	{
		ofs << " "
			<< _receivers->at(i).x() << " "
			<< _receivers->at(i).y() << " "
			<< _receivers->at(i).z() << endl;
	}

	ofs << "VERTICES 1 " << _receiverCount + 1 << endl;
	ofs << _receiverCount << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		((i % 5) == 4) ? ofs << " " << i << endl : ofs << " " << i;
	}
	ofs << endl;

	ofs << "POINT_DATA " << _receiverCount << endl;
	ofs << "FIELD FieldData 4" << endl;
	ofs << "ReceiverID 1 " << _receiverCount << " int" << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		ofs << " " << i << endl;
	}

	ofs << "Eabs 1 " << _receiverCount << " double" << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		double eabs = std::isfinite(_efield->at(i).norm()) ? _efield->at(i).norm() : 0.0;
		ofs << " " << eabs << endl;
	}

	ofs << "EdB 1 " << _receiverCount << " double" << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		double eabs = std::isfinite(_efield->at(i).norm()) ? _efield->at(i).norm() : 0.0;

		double edb = eabs > 0.0 ? 20 * log10(_efield->at(i).norm() * 1e6) : -300.0;
		ofs << " " << edb << endl;
	}

	ofs << "Phase 1 " << _receiverCount << " double" << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		double totalPhase = std::arg(_efield->at(i).x()) + std::arg(_efield->at(i).y()) + std::arg(_efield->at(i).z());
		totalPhase = std::isfinite(totalPhase) ? totalPhase : 0.0;
		if (std::isfinite(totalPhase))
		{
			totalPhase >= TWOPI ? (totalPhase - TWOPI) * (180.0 / PI) : totalPhase * (180.0 / PI);
		}
		else
		{
			totalPhase = 361.0;
		}
		ofs << " " << totalPhase << endl;
	}

	ofs.close();
	cout << "\tSaved " << _receiverCount << " receiver info into " << fname << endl;
	cout << "[Leaving] RTSolver::SaveReceiversAsVtk" << endl;

	return true;
}

bool RTSolver::SaveFieldAsCsv(std::string fname)
{
	using namespace std;
	cout << "[Entering] RTSolver::SaveFieldAsCsv ..." << endl;

	ofstream ofs;
	ofs.open(fname);

	// write header
	ofs << "X,Y,Z,Ex,Ey,Ez" << endl;
	for (int i = 0; i < _receiverCount; i++)
	{
		ofs << _receivers->at(i).x() << ","
			<< _receivers->at(i).y() << ","
			<< _receivers->at(i).z() << ",";

		ofs << showpos
			<< _efield->at(i).x().real() << _efield->at(i).x().imag() << "i" << ","
			<< _efield->at(i).y().real() << _efield->at(i).y().imag() << "i" << ","
			<< _efield->at(i).z().real() << _efield->at(i).z().imag() << "i" << endl;
	}

	cout << "\tSaved " << _receiverCount << " received field into " << fname << endl;
	cout << "[Leaving] RTSolver::SaveFieldAsCsv ..." << endl;
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
		// We can choose to clamp t0 to 1 here. This means that the closest point is the target point
		// However, I choose to ignore it because the reception sphere algorithm checks to see if the capture\
		// point is in LOS of the receiver.
		t0 = 1; // Clamp t0 to 1 because the closest point is the target point
	}

	// Find the optimal reception sphere radius
	// Reference: 191115_Ray_Launching.pdf
	// Note: the unfolded path length is the path length from previous ray segments plus the distance
	//  the point of intersection and the receiver

	// Takahiro's method
	double beta = 1.0 / (3.0 * txTesslation) * acos(-1.0 / sqrt(5.0));
	double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	double sphereRadius = unfoldedLength * tan(beta);

	// Rappaport's approximation
	// double unfoldedLength = pathLength + (sourcePoint + t0 * rayDir).norm();
	// double alpha = (69.0 * PI / 180.0) / txTesslation;
	// double sphereRadius = alpha * unfoldedLength / sqrt(3.0);

	// double unfoldedLength = pathLength + (t0 * rayDir).norm();
	// double phi = (sqrt(5) + 1) / 2;
	// double alpha = acos(phi / (phi * phi + 1)) / (pow(2, txTesslation - 1)) / sqrt(3);
	// double sphereRadius = alpha * unfoldedLength;


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

std::vector<int> GetHitSurfaceIds(const std::vector<Ray>& rayPaths)
{
	std::vector<int> result;
	for (int i = 0; i < rayPaths.size(); i++)
	{
		result.push_back(rayPaths[i].hitSurfaceID);
	}
	return result;
}

std::vector<int> GetHitCoplanarIds(const std::vector<Ray>& rayPaths)
{
	std::vector<int> result;
	for (int i = 0; i < rayPaths.size(); i++)
	{
		result.push_back(rayPaths[i].hitCoplanarId);
	}
	return result;
}

double GetTotalRayLength(const std::vector<Ray>& rayPaths)
{
	if (rayPaths.empty())
	{
		return INF;
	}

	double totalRayLength = 0.0;
	for (int i = 0; i < rayPaths.size(); i++)
	{
		totalRayLength += rayPaths[i].pathLength;
	}
	return totalRayLength;
}

double DistanceToReceiver(const std::vector<Ray>& rayPaths, const Vec3& receiver)
{
	if (rayPaths.empty())
	{
		return INF;
	}

	return (rayPaths.back().targetPoint - receiver).norm();
}

Vec3 ImagePoint(const Triangle& triangle, const Vec3& src)
{
	// We will treat triangle as an infinite plane Ax + By + Cz + d = 0
	// where n = <A,B,C> and d is calcuated by based on a known point on the plane
	// in this case the point is one of the triangle's vertices

	// Algorithm: calculate perpendicular distance between src point and the plane
	// the image point is [src - 2(normal * dist)] 

	// Assumption is that the normal vector is normalized

	double dist = triangle.norm.dot(src - triangle.v1);

	return src - (2*dist*triangle.norm);
}
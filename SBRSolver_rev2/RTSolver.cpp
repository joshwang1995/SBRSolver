#include "RTSolver.h"
#include "CpuTracer.h"
#include "common/ParallelFor.h"
#include "common/Profiling.h"

#include <algorithm>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <map>
#include <memory>
#include <utility>

namespace
{
	// Segments shorter than this (metres) carry no field: 1/L would overflow.
	const double kMinSegmentLength = 1e-9;
}

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
	const std::vector<TransmitterSource>& transmitters,
	int maxReflectionCount,
	int maxTransmissionCount,
	double frequency,
	bool useFresnelCoeff,
	bool powerSum,
	VecVec3& receivers,
	std::vector<Triangle*>& triangleMesh
)
{
	SBR_ZONE_NAMED("ExecuteRayTracing");

	_maxReflectionCount = maxReflectionCount;
	_maxTransmissionCount = maxTransmissionCount;
	// One launch list for all transmitters: every direction remembers the
	// transmitter it starts from, so the launches of all of them run
	// together in the loop below.
	_transmitters = transmitters;
	_shootRayList = new std::vector<Vec3>;
	_launchTx.clear();
	for (size_t t = 0; t < _transmitters.size(); t++)
	{
		std::unique_ptr<std::vector<Vec3>> directions(GenerateRaysOnIcosahedron(_transmitters[t].tessellation, _transmitters[t].position));
		_shootRayList->insert(_shootRayList->end(), directions->begin(), directions->end());
		_launchTx.insert(_launchTx.end(), directions->size(), int(t));
	}
	_receivers = &receivers;
	_receiverCount = int(receivers.size());
	_pathsCount = int(_shootRayList->size());
	_efield = new std::vector<Vec3c>;

	_coplanarGroups.clear();
	for (Triangle* t : triangleMesh)
	{
		if (t->coplanarId < 0) continue;
		if (int(_coplanarGroups.size()) <= t->coplanarId) _coplanarGroups.resize(size_t(t->coplanarId) + 1);
		_coplanarGroups[t->coplanarId].push_back(t);
	}

	InitRayPaths();

	const int threads = (_threadCount > 0) ? _threadCount : DefaultThreadCount();
	std::cout << "\tThreads                  -> " << threads << std::endl;

	Timer timer;
	timer.start();
	// The tracer finds every segment of every launch and the segments that
	// cross each receiver's reception sphere. A candidate path is the chain
	// of segments from the root of its launch to the captured one.
	TracerInput input;
	input.bvh = _bvh;
	input.transmitters = &_transmitters;
	input.launches = _shootRayList;
	input.launchTx = &_launchTx;
	input.receivers = _receivers;
	for (const TransmitterSource& tx : _transmitters)
	{
		input.tanBeta.push_back(ReceptionSlope(tx.tessellation + 1)); // a user may give level 0
	}
	input.maxReflections = _maxReflectionCount;
	input.maxTransmissions = _maxTransmissionCount;
	input.threads = threads;
	TraceResult trace;
	TraceAndCapture(input, trace);
	for (const std::pair<int, int>& candidate : trace.candidates)
	{
		std::vector<Ray> path;
		for (int s = candidate.second; s >= 0; s = trace.segments[s].parent)
		{
			path.push_back(trace.segments[s].ray);
		}
		std::reverse(path.begin(), path.end());
		// The last ray ends at the receiver now, not at the surface it hit.
		path.back().hitSurfaceID = -1;
		path.back().hitCoplanarId = -1;
		_rayPaths[trace.segments[candidate.second].launch][candidate.first].rayPaths.push_back(std::move(path));
	}
	std::cout << "\tSegments                 -> " << trace.segments.size() << std::endl;
	std::cout << "\tCandidates               -> " << trace.candidates.size() << std::endl;
	std::cout << "\tTotal Time in ray launching -> " << timer.getTime() << std::endl;

	// The source field differs per transmitter (power, antenna, frame); the
	// propagation does not.
	std::vector<std::unique_ptr<FieldCompute>> fieldCores;
	for (const TransmitterSource& tx : _transmitters)
	{
		fieldCores.emplace_back(new FieldCompute
		(
			_rayPaths,
			triangleMesh,
			receivers,
			_pathsCount,
			_materialProperties,
			tx.coordSys,
			frequency,
			tx.powerWatts,
			useFresnelCoeff,
			tx.antennaType
		));
	}

	// Each receiver only touches its own column of _rayPaths and its own
	// entry of _efield, so the receivers run independently as well.
	_efield->assign(_receiverCount, Vec3c{ cdouble(0,0), cdouble(0,0), cdouble(0,0) });
	_pathFields.assign(_receiverCount, std::vector<Vec3c>());
	timer.restart();
	ParallelFor(_receiverCount, threads, [&](int k)
	{
		SBR_ZONE_NAMED("ReceiverField");
		{
			SBR_ZONE_NAMED("ImagePathCorrection");
			ImagePathCorrection(k, triangleMesh);
			RemoveDuplicatePath(k);
		}
		// The field of every path is kept (the viewer shows it). The paths
		// of one transmitter add as phasors. Between transmitters the fields
		// add as phasors too, or, for independent sources, the powers add
		// and the field written is the root of the summed squared
		// components, without phase.
		std::vector<Vec3c>& fields = _pathFields[k];
		fields.clear();
		std::vector<Vec3c> perTransmitter(_transmitters.size(), Vec3c::Zero());
		for (int i = 0; i < _pathsCount; i++)
		{
			const FieldCompute& fieldCore = *fieldCores[_launchTx[i]];
			for (const std::vector<Ray>& path : _rayPaths[i][k].rayPaths)
			{
				const Vec3c field = fieldCore.FieldForPath(path);
				fields.push_back(field);
				perTransmitter[_launchTx[i]] += field;
			}
		}
		Vec3c total = Vec3c::Zero();
		if (powerSum)
		{
			Vec3 power = Vec3::Zero();
			for (const Vec3c& field : perTransmitter) power += field.cwiseAbs2();
			total = power.cwiseSqrt().cast<cdouble>();
		}
		else
		{
			for (const Vec3c& field : perTransmitter) total += field;
		}
		(*_efield)[k] = total;
	});
	std::cout << "\tTotal Time in field computation -> " << timer.getTime() << std::endl;

	return _pathsCount * _receiverCount;
}

const std::vector<Vec3c>& RTSolver::Fields() const
{
	static const std::vector<Vec3c> empty;
	return _efield ? *_efield : empty;
}

std::vector<sbr::viewer::ScenePath> RTSolver::CollectPaths() const
{
	using sbr::viewer::InteractionType;
	using sbr::viewer::ScenePath;
	std::vector<ScenePath> result;
	if (_rayPaths == nullptr) return result;
	// _pathFields[j] holds the fields in (launch, path) order: the order of
	// the loops below for a fixed receiver j.
	std::vector<size_t> nextField(_receiverCount, 0);
	for (int i = 0; i < _pathsCount; i++)
	{
		for (int j = 0; j < _receiverCount; j++)
		{
			const Paths& paths = _rayPaths[i][j];
			for (size_t p = 0; p < paths.rayPaths.size(); p++)
			{
				const std::vector<Ray>& rays = paths.rayPaths[p];
				const size_t fieldIndex = nextField[j]++;
				if (rays.empty()) continue;
				ScenePath sp;
				sp.receiverId = j;
				sp.launchId = i;
				sp.transmitterId = _launchTx[i];
				sp.points.reserve(rays.size() + 1);
				sp.points.push_back(rays[0].sourcePoint);
				for (size_t r = 0; r < rays.size(); r++)
				{
					sp.points.push_back(rays[r].targetPoint);
					sp.length += (rays[r].targetPoint - rays[r].sourcePoint).norm();
					if (r > 0)
					{
						sp.interactions.push_back(rays[r].reflectionMaterialId >= 0 ? InteractionType::Reflection : InteractionType::Transmission);
					}
				}
				if (j < int(_pathFields.size()) && fieldIndex < _pathFields[j].size()) sp.field = _pathFields[j][fieldIndex];
				result.push_back(sp);
			}
		}
	}
	return result;
}

void RTSolver::ImagePathCorrection(int receiverId, const std::vector<Triangle*>& triangleMesh)
{
	for (int i = 0; i < _pathsCount; i++)
	{
		const Vec3& txPoint = _transmitters[_launchTx[i]].position;
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
	// imgSources[i] is the surface at which segment i starts (-1 for the
	// transmitter) and the image of the transmitter that segment i appears
	// to come from.
	std::vector<std::pair<int, Vec3>> imgSources = GetImageSources(multiPath, txPoint, triangleMesh);
	if (imgSources.size() != multiPath.size())
	{
		return false;
	}

	// Walk back from the receiver: each interaction point is where the line
	// from the segment's image source to the next point crosses the whole
	// coplanar surface (the launched ray may have hit another triangle of it).
	Vec3 targetPoint = _receivers->at(receiverId);
	for (int i = int(multiPath.size()) - 1; i >= 1; i--)
	{
		const int surfaceId = imgSources[i].first;
		if (surfaceId < 0 || surfaceId >= int(triangleMesh.size())) return false;
		Vec3 point;
		Triangle* hit = IntersectSurface(triangleMesh[surfaceId]->coplanarId, imgSources[i].second, targetPoint, point);
		if (hit == nullptr) return false;

		Ray& ray = multiPath[i];
		ray.sourcePoint = point;
		ray.targetPoint = targetPoint;
		ray.pathLength = (targetPoint - point).norm();
		Ray& previous = multiPath[i - 1];
		previous.hitSurfaceID = hit->triangleId;
		previous.hitCoplanarId = hit->coplanarId;
		previous.angleFromSurfaceNormal = AngleBetween(hit->norm, targetPoint - point);
		targetPoint = point;
	}
	multiPath[0].sourcePoint = txPoint;
	multiPath[0].targetPoint = targetPoint;
	multiPath[0].pathLength = (targetPoint - txPoint).norm();

	// The exact path exists only if nothing else stands in the way of any of
	// its segments. The surfaces at a segment's ends are not obstacles. A
	// segment of zero length (receiver on the transmitter or on an
	// interaction point) has no field and is dropped, as the reference does.
	for (size_t i = 0; i < multiPath.size(); i++)
	{
		if (multiPath[i].pathLength < kMinSegmentLength) return false;
		const int startSurface = (i == 0) ? -1 : multiPath[i - 1].hitCoplanarId;
		const int endSurface = (i + 1 == multiPath.size()) ? -1 : multiPath[i].hitCoplanarId;
		if (SegmentBlocked(multiPath[i].sourcePoint, multiPath[i].targetPoint, startSurface, endSurface))
		{
			return false;
		}
	}
	return true;
}

bool RTSolver::DirectPathCorrection(std::vector<Ray>& directPath, int receiverId, const Vec3& txPoint)
{
	if (directPath.size() != 1) { return false; }

	const Vec3& receiver = _receivers->at(receiverId);
	if ((receiver - txPoint).norm() < kMinSegmentLength || SegmentBlocked(txPoint, receiver, -1, -1))
	{
		return false;
	}
	directPath[0].sourcePoint = txPoint;
	directPath[0].targetPoint = receiver;
	directPath[0].pathLength = (receiver - txPoint).norm();
	return true;
}

bool RTSolver::SegmentBlocked(const Vec3& from, const Vec3& to, int skipCoplanarA, int skipCoplanarB) const
{
	const Vec3 segment = to - from;
	const double length = segment.norm();
	if (length <= 0.0) return false;
	// Hits closer than a hair's breadth to either end are the end surfaces
	// themselves (or the receiver sitting on a surface), not obstacles.
	const double margin = std::max(1e-9 * length, 1e-9);
	return _bvh->RayOccluded(0, from, segment, length - margin, margin, skipCoplanarA, skipCoplanarB);
}

Triangle* RTSolver::IntersectSurface(int coplanarId, const Vec3& origin, const Vec3& target, Vec3& point) const
{
	if (coplanarId < 0 || coplanarId >= int(_coplanarGroups.size())) return nullptr;
	const Vec3 direction = target - origin;
	const double length = direction.norm();
	for (Triangle* t : _coplanarGroups[coplanarId])
	{
		HitInfo info;
		if (t->RayIntersects(origin, direction, info) && info.distance < length)
		{
			point = info.pointIntersect;
			return t;
		}
	}
	return nullptr;
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
	// A path is its transmitter and its sequence of surfaces: the same
	// sequence from two transmitters is two paths.
	std::map<std::vector<int>, int> surfaceIdMap;

	// std::unordered_set<std::vector<int>> hitIdSet;

	for (int i = 0; i < _pathsCount; i++)
	{
		auto it = _rayPaths[i][receiverId].rayPaths.begin();
		while (it != _rayPaths[i][receiverId].rayPaths.end())
		{
			// std::vector<int> key = GetHitSurfaceIds(*it);
			std::vector<int> key = GetHitCoplanarIds(*it);
			key.insert(key.begin(), _launchTx[i]);
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
	ofs << "FIELD FieldData 6" << endl;
	ofs << "LaunchID 1 " << launchIds.size() << " int" << endl;
	for (int i = 0; i < launchIds.size(); i++)
	{
		ofs << " " << launchIds[i] << endl;
	}

	ofs << "TransmitterID 1 " << launchIds.size() << " int" << endl;
	for (int i = 0; i < launchIds.size(); i++)
	{
		ofs << " " << _launchTx[launchIds[i]] << endl;
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
		double edb = std::isfinite(_efield->at(receiverIds[i]).norm()) ? 20 * log10(_efield->at(receiverIds[i]).norm() * 1e6) : -300.0;
		edb = std::isfinite(edb) ? edb : -300.0;
		ofs << " " << edb << endl;
	}

	ofs << "Phase 1 " << receiverIds.size() << " double" << endl;
	for (int i = 0; i < receiverIds.size(); i++)
	{
		double totalPhase = std::arg(_efield->at(receiverIds[i]).x()) + std::arg(_efield->at(receiverIds[i]).y()) + std::arg(_efield->at(receiverIds[i]).z());
		// Written in radians, unwrapped (non-finite values are clamped to 0).
		totalPhase = std::isfinite(totalPhase) ? totalPhase : 0.0;
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
		// Written in radians, unwrapped (non-finite values are clamped to 0).
		totalPhase = std::isfinite(totalPhase) ? totalPhase : 0.0;
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
	// 12 significant digits: enough to compare with the reference solver's 8
	// without the printing itself adding a visible difference.
	ofs << std::setprecision(12);
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

bool RTSolver::SaveIcosahedronAsVtk(std::string fname)
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

	// The first points are the transmitters, then one point per launch
	// direction at unit distance from its transmitter.
	const size_t launches = _shootRayList->size();
	const size_t origins = _transmitters.size();
	ofs << "POINTS " << launches + origins << " float" << endl;
	for (const TransmitterSource& tx : _transmitters)
	{
		ofs << " " << tx.position.x() << " " << tx.position.y() << " " << tx.position.z() << endl;
	}
	for (size_t i = 0; i < launches; i++)
	{
		const Vec3& origin = _transmitters[_launchTx[i]].position;
		const Vec3& v = _shootRayList->at(i);
		ofs << " "
			<< origin.x() + 1.0 * v.x() << " "
			<< origin.y() + 1.0 * v.y() << " "
			<< origin.z() + 1.0 * v.z() << endl;
	}

	ofs << "LINES " << launches << " " << 3 * launches << endl;
	for (size_t i = 0; i < launches; i++)
	{
		ofs << " 2 " << _launchTx[i] << " " << origins + i << endl;
	}

	ofs << "CELL_DATA " << launches << endl;
	ofs << "FIELD FieldData 2" << endl;
	ofs << "LaunchID 1 " << launches << " int" << endl;
	for (size_t i = 0; i < launches; i++)
	{
		ofs << " " << i << endl;
	}
	ofs << "TransmitterID 1 " << launches << " int" << endl;
	for (size_t i = 0; i < launches; i++)
	{
		ofs << " " << _launchTx[i] << endl;
	}

	ofs.close();
	cout << "\tSaved " << launches << " launch directions into " << fname << endl;
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
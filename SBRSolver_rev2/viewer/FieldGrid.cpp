#include "viewer/FieldGrid.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

namespace sbr
{
namespace viewer
{

namespace
{
	// A query closer than this fraction of the sample spacing to a sample takes
	// that sample's value unchanged (the inverse-distance weight would blow up).
	const double kCoincidenceFraction = 1e-12;

	const double kInfinity = std::numeric_limits<double>::infinity();

	// Balanced k-d tree over the sample positions, stored implicitly: the
	// median of every index range [lo, hi) of `_order` is the node, the ranges
	// before and after it are its children. Build O(n log n), a k-nearest
	// query O(log n + k) on average whatever the distance to the samples, so
	// grid nodes far from a sparse or tilted receiver layout (which still need
	// the nearest sample) cost no more than nodes among the samples.
	class KdTree
	{
	public:
		typedef std::pair<double, int> Hit;   // squared distance, sample index

		// The coordinates are copied into a flat array: the search loop then
		// runs on plain doubles (fast in unoptimised builds too) and the tree
		// does not refer to the caller's vector.
		explicit KdTree(const std::vector<Vec3>& points)
			: _coords(points.size() * 3), _order(points.size()), _axis(points.size(), 0)
		{
			for (size_t i = 0; i < points.size(); i++)
			{
				_coords[3 * i] = points[i].x();
				_coords[3 * i + 1] = points[i].y();
				_coords[3 * i + 2] = points[i].z();
				_order[i] = int(i);
			}
			Build(0, int(_order.size()));
		}

		// The (at most) k samples nearest to q whose squared distance does not
		// exceed maxDist2, nearest first; sample `exclude` (-1: none) is skipped.
		void Nearest(const Vec3& q, int k, double maxDist2, int exclude, std::vector<Hit>& hits) const
		{
			hits.clear();
			if (k <= 0 || _order.empty()) return;
			const double query[3] = {q.x(), q.y(), q.z()};
			Search(query, 0, int(_order.size()), k, maxDist2, exclude, hits);
		}

	private:
		double Coord(int index, int axis) const
		{
			return _coords[size_t(index) * 3 + size_t(axis)];
		}

		void Build(int lo, int hi)
		{
			while (hi - lo > 1)
			{
				// Split along the axis of largest spread, so flat and linear
				// layouts are not split along a direction without extent.
				double lower[3], upper[3];
				for (int a = 0; a < 3; a++) lower[a] = upper[a] = Coord(_order[lo], a);
				for (int i = lo + 1; i < hi; i++)
				{
					for (int a = 0; a < 3; a++)
					{
						const double c = Coord(_order[i], a);
						lower[a] = std::min(lower[a], c);
						upper[a] = std::max(upper[a], c);
					}
				}
				int axis = 0;
				for (int a = 1; a < 3; a++)
				{
					if (upper[a] - lower[a] > upper[axis] - lower[axis]) axis = a;
				}

				const int mid = lo + (hi - lo) / 2;
				std::nth_element(_order.begin() + lo, _order.begin() + mid, _order.begin() + hi,
				                 [&](int a, int b) { return Coord(a, axis) < Coord(b, axis); });
				_axis[mid] = static_cast<unsigned char>(axis);

				// Recurse into the smaller half, loop on the larger one.
				if (mid - lo < hi - (mid + 1))
				{
					Build(lo, mid);
					lo = mid + 1;
				}
				else
				{
					Build(mid + 1, hi);
					hi = mid;
				}
			}
		}

		void Search(const double* q, int lo, int hi, int k, double maxDist2, int exclude, std::vector<Hit>& hits) const
		{
			while (lo < hi)
			{
				const int mid = lo + (hi - lo) / 2;
				const int index = _order[mid];
				const double* p = &_coords[size_t(index) * 3];
				const double dx = p[0] - q[0];
				const double dy = p[1] - q[1];
				const double dz = p[2] - q[2];
				const double dist2 = dx * dx + dy * dy + dz * dz;
				if (index != exclude && dist2 <= maxDist2 && (int(hits.size()) < k || dist2 < hits.back().first))
				{
					const Hit hit(dist2, index);
					hits.insert(std::upper_bound(hits.begin(), hits.end(), hit,
					                             [](const Hit& a, const Hit& b) { return a.first < b.first; }),
					            hit);
					if (int(hits.size()) > k) hits.pop_back();
				}

				const int axis = _axis[mid];
				const double diff = q[axis] - p[axis];
				int nearLo, nearHi, farLo, farHi;
				if (diff < 0.0)
				{
					nearLo = lo;
					nearHi = mid;
					farLo = mid + 1;
					farHi = hi;
				}
				else
				{
					nearLo = mid + 1;
					nearHi = hi;
					farLo = lo;
					farHi = mid;
				}
				Search(q, nearLo, nearHi, k, maxDist2, exclude, hits);

				// Every point on the far side is at least |diff| away.
				const double bound = (int(hits.size()) < k) ? maxDist2 : hits.back().first;
				if (diff * diff > bound) return;
				lo = farLo;
				hi = farHi;
			}
		}

		std::vector<double> _coords;        // x, y, z of every sample
		std::vector<int> _order;            // sample indices, partitioned by the tree
		std::vector<unsigned char> _axis;   // split axis of the node at that position of _order
	};

	double MedianNearestNeighbourDistance(const KdTree& tree, const std::vector<Vec3>& points)
	{
		const size_t n = points.size();
		if (n < 2) return 0.0;

		std::vector<double> distances(n, 0.0);
		std::vector<KdTree::Hit> hits;
		for (size_t i = 0; i < n; i++)
		{
			tree.Nearest(points[i], 1, kInfinity, int(i), hits);
			if (!hits.empty()) distances[i] = std::sqrt(hits[0].first);
		}

		const size_t mid = n / 2;
		std::nth_element(distances.begin(), distances.begin() + mid, distances.end());
		double median = distances[mid];
		if (n % 2 == 0)
		{
			median = 0.5 * (median + *std::max_element(distances.begin(), distances.begin() + mid));
		}
		return median;
	}

	// Interpolation with an existing tree and sample spacing (see the header).
	void InterpolateWithTree(const KdTree& tree, double spacing, const std::vector<Vec3c>& values,
	                         const std::vector<Vec3>& queries, const InterpolationSettings& settings,
	                         std::vector<Vec3c>& out, std::vector<unsigned char>& valid)
	{
		out.assign(queries.size(), Vec3c::Zero());
		valid.assign(queries.size(), 0);
		if (values.empty()) return;

		// Without a length scale (one sample, or coincident samples) the search
		// radius is unlimited, so every query gets the sample.
		const bool limited = settings.radiusFactor > 0.0 && spacing > 0.0;
		const double radius2 = limited ? (settings.radiusFactor * spacing) * (settings.radiusFactor * spacing) : kInfinity;
		const double coincidence2 = (kCoincidenceFraction * spacing) * (kCoincidenceFraction * spacing);
		const bool nearestOnly = settings.method == InterpolationMethod::Nearest;
		const int neighbours = nearestOnly ? 1 : std::max(1, settings.neighbours);

		std::vector<KdTree::Hit> hits;
		for (size_t i = 0; i < queries.size(); i++)
		{
			tree.Nearest(queries[i], neighbours, radius2, -1, hits);
			if (hits.empty())
			{
				// Nothing within the radius: flag the node, show the nearest sample.
				tree.Nearest(queries[i], 1, kInfinity, -1, hits);
				out[i] = values[hits[0].second];
				continue;
			}
			valid[i] = 1;
			if (nearestOnly || hits[0].first <= coincidence2)
			{
				out[i] = values[hits[0].second];
				continue;
			}

			// Shepard: w = 1 / d^power, scaled by the nearest distance so the
			// weights stay in (0, 1] whatever the power.
			Vec3c sum = Vec3c::Zero();
			double weightSum = 0.0;
			for (const KdTree::Hit& hit : hits)
			{
				const double w = std::pow(hits[0].first / hit.first, 0.5 * settings.power);
				sum += w * values[hit.second];
				weightSum += w;
			}
			out[i] = sum / weightSum;
		}
	}

	// Nodes on a side of length `length` when the longest side (`longest`) has
	// `resolution` nodes: the same node spacing on every side, at least `minimum`.
	int ProportionalNodeCount(double length, double longest, int resolution, int minimum)
	{
		if (longest <= 0.0) return minimum;
		const double step = longest / double(resolution - 1);
		const int count = int(std::lround(length / step)) + 1;
		return std::max(minimum, count);
	}

	double PaddingFor(double spacing)
	{
		return spacing > 0.0 ? 0.5 * spacing : 0.5;
	}

	// Flip a unit vector so that its component of largest magnitude is
	// positive: the eigenvector signs are arbitrary, this keeps the frame
	// stable between runs and aligned with the global axes when possible.
	Vec3 OrientAxis(const Vec3& axis)
	{
		int largest = 0;
		axis.cwiseAbs().maxCoeff(&largest);
		return axis[largest] < 0.0 ? Vec3(-axis) : axis;
	}
}

double SampleSpacing(const std::vector<Vec3>& samples)
{
	if (samples.size() < 2) return 0.0;
	const KdTree tree(samples);
	return MedianNearestNeighbourDistance(tree, samples);
}

PlaneFit FitPlane(const std::vector<Vec3>& samples)
{
	PlaneFit fit;
	if (samples.empty()) return fit;

	Vec3 centroid = Vec3::Zero();
	for (const Vec3& p : samples) centroid += p;
	centroid /= double(samples.size());
	fit.origin = centroid;
	if (samples.size() < 2) return fit;

	Mat3 covariance = Mat3::Zero();
	for (const Vec3& p : samples)
	{
		const Vec3 d = p - centroid;
		covariance += d * d.transpose();
	}
	covariance /= double(samples.size());

	// Eigenvalues ascending: the principal axes are the last two columns. For
	// collinear samples the second eigenvalue is 0 and v is an arbitrary
	// direction perpendicular to u (extentV = 0): a degenerate but usable plane.
	Eigen::SelfAdjointEigenSolver<Mat3> solver(covariance);
	if (solver.info() != Eigen::Success) return fit;
	const Vec3 u = OrientAxis(solver.eigenvectors().col(2));
	const Vec3 v = OrientAxis(solver.eigenvectors().col(1));
	const Vec3 normal = u.cross(v).normalized();

	double sMin = kInfinity, sMax = -kInfinity, tMin = kInfinity, tMax = -kInfinity, thickness = 0.0;
	for (const Vec3& p : samples)
	{
		const Vec3 d = p - centroid;
		const double s = d.dot(u);
		const double t = d.dot(v);
		sMin = std::min(sMin, s);
		sMax = std::max(sMax, s);
		tMin = std::min(tMin, t);
		tMax = std::max(tMax, t);
		thickness = std::max(thickness, std::abs(d.dot(normal)));
	}

	// Coincident samples span nothing: no plane, only the centroid.
	if (!(sMax - sMin > 0.0)) return fit;
	fit.ok = true;
	fit.u = u;
	fit.v = v;
	fit.normal = normal;
	fit.extentU = sMax - sMin;
	fit.extentV = tMax - tMin;
	fit.thickness = thickness;
	return fit;
}

void Interpolate(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, const std::vector<Vec3>& queries,
                 const InterpolationSettings& settings, std::vector<Vec3c>& out, std::vector<unsigned char>& valid)
{
	assert(values.size() == samples.size());
	out.clear();
	valid.clear();
	if (samples.empty()) return;

	const KdTree tree(samples);
	const double spacing = MedianNearestNeighbourDistance(tree, samples);
	InterpolateWithTree(tree, spacing, values, queries, settings, out, valid);
}

PlaneGrid InterpolateOnPlane(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, const PlaneFit& plane,
                             int resolution, const InterpolationSettings& settings)
{
	assert(values.size() == samples.size());
	PlaneGrid grid;
	grid.plane = plane;
	if (samples.empty()) return grid;

	const KdTree tree(samples);
	const double spacing = MedianNearestNeighbourDistance(tree, samples);
	const double pad = PaddingFor(spacing);

	// Sample extents in the plane's frame, padded by half a spacing on every
	// side; a side without extent (a receiver line) becomes one spacing wide.
	double sMin = kInfinity, sMax = -kInfinity, tMin = kInfinity, tMax = -kInfinity;
	for (const Vec3& p : samples)
	{
		const Vec3 d = p - plane.origin;
		const double s = d.dot(plane.u);
		const double t = d.dot(plane.v);
		sMin = std::min(sMin, s);
		sMax = std::max(sMax, s);
		tMin = std::min(tMin, t);
		tMax = std::max(tMax, t);
	}
	sMin -= pad;
	sMax += pad;
	tMin -= pad;
	tMax += pad;
	const double lengthS = sMax - sMin;
	const double lengthT = tMax - tMin;
	const double longest = std::max(lengthS, lengthT);

	resolution = std::max(2, resolution);
	grid.nu = ProportionalNodeCount(lengthS, longest, resolution, 2);
	grid.nv = ProportionalNodeCount(lengthT, longest, resolution, 2);

	grid.nodes.resize(size_t(grid.nu) * size_t(grid.nv));
	for (int j = 0; j < grid.nv; j++)
	{
		const double t = tMin + lengthT * double(j) / double(grid.nv - 1);
		for (int i = 0; i < grid.nu; i++)
		{
			const double s = sMin + lengthS * double(i) / double(grid.nu - 1);
			grid.nodes[grid.Index(i, j)] = plane.origin + s * plane.u + t * plane.v;
		}
	}

	InterpolateWithTree(tree, spacing, values, grid.nodes, settings, grid.values, grid.valid);
	return grid;
}

VolumeGrid InterpolateInVolume(const std::vector<Vec3>& samples, const std::vector<Vec3c>& values, int resolution,
                               const InterpolationSettings& settings)
{
	assert(values.size() == samples.size());
	VolumeGrid grid;
	if (samples.empty()) return grid;

	const KdTree tree(samples);
	const double spacing = MedianNearestNeighbourDistance(tree, samples);
	const double pad = PaddingFor(spacing);

	Vec3 lower = samples[0];
	Vec3 upper = samples[0];
	for (const Vec3& p : samples)
	{
		lower = lower.cwiseMin(p);
		upper = upper.cwiseMax(p);
	}
	const Vec3 extent = upper - lower;   // before padding
	grid.boundMin = lower - Vec3::Constant(pad);
	grid.boundMax = upper + Vec3::Constant(pad);
	const Vec3 length = grid.boundMax - grid.boundMin;

	// The longest side gets `resolution` nodes; a side thinner than one sample
	// spacing (a flat or linear layout) gets a single node at its centre, the
	// others the same node spacing as the longest side.
	int longestAxis = 0;
	length.maxCoeff(&longestAxis);
	resolution = std::max(2, resolution);
	int counts[3];
	for (int a = 0; a < 3; a++)
	{
		const bool thin = spacing > 0.0 ? extent[a] < spacing * (1.0 - 1e-9) : extent[a] <= 0.0;
		if (a == longestAxis)
		{
			counts[a] = resolution;
		}
		else if (thin)
		{
			counts[a] = 1;
		}
		else
		{
			counts[a] = ProportionalNodeCount(length[a], length[longestAxis], resolution, 2);
		}
	}
	grid.nx = counts[0];
	grid.ny = counts[1];
	grid.nz = counts[2];

	// Node coordinates along each axis.
	std::vector<double> coords[3];
	for (int a = 0; a < 3; a++)
	{
		coords[a].resize(size_t(counts[a]));
		if (counts[a] == 1)
		{
			coords[a][0] = 0.5 * (grid.boundMin[a] + grid.boundMax[a]);
		}
		else
		{
			for (int i = 0; i < counts[a]; i++)
			{
				coords[a][size_t(i)] = grid.boundMin[a] + length[a] * double(i) / double(counts[a] - 1);
			}
		}
	}

	grid.nodes.resize(size_t(grid.nx) * size_t(grid.ny) * size_t(grid.nz));
	for (int k = 0; k < grid.nz; k++)
	{
		for (int j = 0; j < grid.ny; j++)
		{
			for (int i = 0; i < grid.nx; i++)
			{
				grid.nodes[grid.Index(i, j, k)] = Vec3(coords[0][size_t(i)], coords[1][size_t(j)], coords[2][size_t(k)]);
			}
		}
	}

	InterpolateWithTree(tree, spacing, values, grid.nodes, settings, grid.values, grid.valid);
	return grid;
}

} // namespace viewer
} // namespace sbr

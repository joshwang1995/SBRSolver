#include "ReceiverGrid.h"

#include <algorithm>
#include <cmath>

namespace
{
	const double kReceiversPerCell = 8.0;
	// The cell table may not grow past this many cells per receiver, whatever
	// the layout's shape.
	const double kMaxCellsPerReceiver = 4.0;
}

ReceiverGrid::ReceiverGrid(const VecVec3& receivers)
{
	const int count = int(receivers.size());
	if (count == 0)
	{
		_first.assign(1, 0);
		_count.assign(1, 0);
		return;
	}

	Vec3 low = receivers[0];
	Vec3 high = receivers[0];
	for (const Vec3& r : receivers)
	{
		low = low.cwiseMin(r);
		high = high.cwiseMax(r);
	}
	_origin = low;
	_extent = high - low;

	// Cell size from the axes with an extent only: a plane or a line of
	// receivers has a zero extent on the other axes.
	int axes = 0;
	double product = 1.0;
	for (int axis = 0; axis < 3; axis++)
	{
		if (_extent[axis] > 0.0)
		{
			axes++;
			product *= _extent[axis];
		}
	}
	if (axes > 0)
	{
		double cell = std::pow(product * kReceiversPerCell / count, 1.0 / axes);
		for (;;)
		{
			double total = 1.0;
			for (int axis = 0; axis < 3; axis++)
			{
				_cells[axis] = (_extent[axis] > 0.0) ? std::max(1, int(std::ceil(_extent[axis] / cell))) : 1;
				_cellSize[axis] = (_extent[axis] > 0.0) ? cell : 0.0;
				total *= _cells[axis];
			}
			if (total <= kMaxCellsPerReceiver * count) break;
			cell *= 1.25;
		}
	}

	const int cells = _cells[0] * _cells[1] * _cells[2];
	std::vector<int> cellOf(count);
	_count.assign(cells, 0);
	for (int i = 0; i < count; i++)
	{
		cellOf[i] = CellIndex(CellCoordinate(0, receivers[i].x()), CellCoordinate(1, receivers[i].y()), CellCoordinate(2, receivers[i].z()));
		_count[cellOf[i]]++;
	}
	_first.assign(cells, 0);
	for (int c = 1; c < cells; c++) _first[c] = _first[c - 1] + _count[c - 1];
	std::vector<int> next(_first);
	_sorted.resize(count);
	for (int i = 0; i < count; i++) _sorted[next[cellOf[i]]++] = i;
}

int ReceiverGrid::CellCoordinate(int axis, double value) const
{
	if (_cells[axis] == 1) return 0;
	const int c = int(std::floor((value - _origin[axis]) / _cellSize[axis]));
	return std::min(_cells[axis] - 1, std::max(0, c));
}

void ReceiverGrid::Query(const Mat23& box, std::vector<int>& out) const
{
	if (_sorted.empty()) return;
	int low[3];
	int high[3];
	for (int axis = 0; axis < 3; axis++)
	{
		if (box(1, axis) < _origin[axis] || box(0, axis) > _origin[axis] + _extent[axis]) return;
		low[axis] = CellCoordinate(axis, box(0, axis));
		high[axis] = CellCoordinate(axis, box(1, axis));
	}
	for (int z = low[2]; z <= high[2]; z++)
	{
		for (int y = low[1]; y <= high[1]; y++)
		{
			for (int x = low[0]; x <= high[0]; x++)
			{
				const int cell = CellIndex(x, y, z);
				out.insert(out.end(), _sorted.begin() + _first[cell], _sorted.begin() + _first[cell] + _count[cell]);
			}
		}
	}
}

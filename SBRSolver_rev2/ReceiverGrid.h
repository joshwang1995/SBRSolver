#pragma once

#include <vector>

#include "common/VecMatDef.h"

// Uniform grid over the receivers' bounding box, so that a capture test
// visits only the receivers near a segment. A cell holds about eight
// receivers along the axes that have an extent. An axis without extent (a
// plane or a line of receivers) has one cell.
class ReceiverGrid
{
public:
	explicit ReceiverGrid(const VecVec3& receivers);

	/**
	 * @brief Appends the indices of the receivers in every cell that `box`
	 * overlaps. Row 0 of `box` is the minimum corner, row 1 the maximum.
	 * @param box The query box in metres.
	 * @param out Receives the indices; it is not cleared first.
	 */
	void Query(const Mat23& box, std::vector<int>& out) const;

private:
	int CellCoordinate(int axis, double value) const;
	int CellIndex(int x, int y, int z) const { return (z * _cells[1] + y) * _cells[0] + x; }

	Vec3 _origin = Vec3::Zero();   // minimum corner of the grid
	Vec3 _extent = Vec3::Zero();   // size of the grid along every axis
	Vec3 _cellSize = Vec3::Zero(); // zero on an axis with one cell
	int _cells[3] = { 1, 1, 1 };
	std::vector<int> _first;       // per cell: first entry in _sorted
	std::vector<int> _count;       // per cell: number of entries
	std::vector<int> _sorted;      // receiver indices grouped by cell
};

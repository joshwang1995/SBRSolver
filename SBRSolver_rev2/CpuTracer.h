#pragma once

#include "Tracer.h"
#include "common/VecMatDef.h"

// The CPU ray tracer: the reference implementation of the trace and capture
// stages. It launches every ray of every transmitter, follows the reflected
// and transmitted rays up to the bounce limits, and finds the segments that
// pass through each receiver's reception sphere.

/** @brief Traces every launch into `result.segments` and fills `result.launchFirst`. */
void TraceSegments(const TracerInput& input, TraceResult& result);

/** @brief Fills `result.candidates` from `result.segments` and the receivers. */
void CaptureCandidates(const TracerInput& input, TraceResult& result);

/** @brief Runs TraceSegments and CaptureCandidates. */
void TraceAndCapture(const TracerInput& input, TraceResult& result);

/**
 * @brief Slope tan(beta) of a transmitter's reception sphere: the sphere
 * radius is the unfolded path length times this value.
 * @param tessellation The transmitter's tessellation level plus one.
 */
double ReceptionSlope(int tessellation);

/**
 * @brief True when the segment from `sourcePoint` to `targetPoint` passes
 * through the reception sphere at `sphereCenter`.
 * @param pathLength Unfolded length from the transmitter to `sourcePoint`.
 * @param tanBeta Slope of the reception sphere (ReceptionSlope).
 * @param capturePoint Receives the point of the segment closest to the centre.
 */
bool HitReceptionSphere
(
	const Vec3& sourcePoint,
	const Vec3& targetPoint,
	const Vec3& sphereCenter,
	double pathLength,
	double tanBeta,
	Vec3& capturePoint
);

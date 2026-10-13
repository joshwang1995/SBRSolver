// Copyright NVIDIA Corporation 2002-2005
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of NVIDIA CORPORATION nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Derived from the Timer class in the NVIDIA samples; ported to
// std::chrono::steady_clock so it needs no platform-specific code.

#pragma once

#include <chrono>

/*! \brief A simple, portable timer class.
 *
 * Measures wall-clock time intervals in seconds using std::chrono, so it
 * behaves identically on Windows, Linux and macOS. The timer can be started
 * and stopped several times and accumulates the time elapsed between the
 * start() and stop() calls.
 */
class Timer
{
public:
	//! Constructs a Timer, but does not start it yet.
	Timer();

	//! Starts the timer (no-op if already running).
	void start();

	//! Stops the timer and accumulates the elapsed time.
	void stop();

	//! Resets the accumulated time and stops the timer.
	void reset();

	//! Resets the timer and starts it.
	void restart();

	//! Returns the accumulated time in seconds (including the running interval).
	double getTime() const;

	//! Returns whether the timer is currently running.
	bool isRunning() const { return m_running; }

private:
	using Clock = std::chrono::steady_clock;

	Clock::time_point m_begin;
	bool m_running;
	double m_seconds;
};

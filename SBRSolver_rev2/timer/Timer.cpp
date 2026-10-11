// Derived from the NVIDIA nvpro-pipeline Timer (https://github.com/nvpro-pipeline/pipeline),
// see the licence notice in Timer.h. Ported to std::chrono::steady_clock.
#include "Timer.h"

Timer::Timer()
	: m_begin()
	, m_running(false)
	, m_seconds(0.0)
{
}

void Timer::start()
{
	if (!m_running)
	{
		m_running = true;
		// starting a timer: store the starting time last
		m_begin = Clock::now();
	}
}

void Timer::stop()
{
	// stopping a timer: store the stopping time first
	const Clock::time_point now = Clock::now();
	if (m_running)
	{
		m_seconds += std::chrono::duration<double>(now - m_begin).count();
		m_running = false;
	}
}

void Timer::reset()
{
	m_running = false;
	m_seconds = 0.0;
}

void Timer::restart()
{
	reset();
	start();
}

double Timer::getTime() const
{
	if (m_running)
	{
		return m_seconds + std::chrono::duration<double>(Clock::now() - m_begin).count();
	}
	return m_seconds;
}

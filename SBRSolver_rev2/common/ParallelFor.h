#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <functional>
#include <iostream>
#include <mutex>
#include <system_error>
#include <thread>
#include <vector>

#include "common/Profiling.h"

// Number of worker threads to use when the configuration asks for "all
// cores" (0): the hardware concurrency, or 1 if it is unknown.
inline int DefaultThreadCount()
{
	const unsigned int n = std::thread::hardware_concurrency();
	return n == 0 ? 1 : int(n);
}

// Runs body(i) for every i in [0, count) on `threads` threads (the calling
// thread included). Iterations are handed out dynamically, so uneven work per
// iteration balances itself; each iteration must only write state that belongs
// to its own index. An exception thrown by any iteration is rethrown on the
// calling thread once every worker has stopped. If the system refuses to
// start a worker thread, the loop runs on the threads that could be created.
inline void ParallelFor(int count, int threads, const std::function<void(int)>& body)
{
	if (count <= 0) return;
	threads = std::max(1, std::min(threads, count));

	if (threads == 1)
	{
		for (int i = 0; i < count; i++) body(i);
		return;
	}

	// 64-bit so the hand-out can run `threads` past `count` without wrapping
	// for any int count.
	std::atomic<long long> next{ 0 };
	std::atomic<bool> failed{ false };
	std::exception_ptr error;
	std::mutex errorMutex;

	auto worker = [&]()
	{
		try
		{
			for (long long i = next++; i < count && !failed; i = next++)
			{
				body(int(i));
			}
		}
		catch (...)
		{
			std::lock_guard<std::mutex> lock(errorMutex);
			if (!error) error = std::current_exception();
			failed = true;
		}
	};

	std::vector<std::thread> pool;
	pool.reserve(threads - 1);
	for (int t = 1; t < threads; t++)
	{
		try
		{
			pool.emplace_back([&worker]()
			{
				SBR_THREAD_NAME("sbr worker");   // the calling thread keeps its name
				worker();
			});
		}
		catch (...)
		{
			// std::thread's constructor throws (std::system_error, or
			// std::bad_alloc for its state) when the OS refuses another thread:
			// RLIMIT_NPROC, a cgroup pids limit, an address-space limit. No
			// thread was started by the failed call, so carry on with the
			// workers that exist plus the calling thread; the counter balances
			// the work over however many there are.
			std::cerr << "\tnote: only " << pool.size() + 1 << " of " << threads
				<< " threads could be created; continuing with those" << std::endl;
			break;
		}
	}
	worker();
	for (std::thread& t : pool) t.join();

	if (error) std::rethrow_exception(error);
}

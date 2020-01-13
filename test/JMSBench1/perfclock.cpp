#include "perfclock.h"



perfclock::perfclock()
{
}


perfclock::~perfclock()
{
}

std::chrono::time_point<std::chrono::steady_clock> perfclock::start()
{
	cl_start = std::chrono::steady_clock::now();
	return cl_start;
}

long long perfclock::getnsec(bool restart)
{
	std::chrono::time_point<std::chrono::steady_clock> cl_now = std::chrono::steady_clock::now();
	long long dur = std::chrono::duration_cast<std::chrono::nanoseconds>(cl_now - cl_start).count();
	if (restart) cl_start = cl_now;
	return dur;
}

long long perfclock::getusec(bool restart)
{
	std::chrono::time_point<std::chrono::steady_clock> cl_now = std::chrono::steady_clock::now();
	long long dur = std::chrono::duration_cast<std::chrono::microseconds>(cl_now - cl_start).count();
	if (restart) cl_start = cl_now;
	return dur;
}

long long perfclock::getmsec(bool restart)
{
	std::chrono::time_point<std::chrono::steady_clock> cl_now = std::chrono::steady_clock::now();
	long long dur = std::chrono::duration_cast<std::chrono::milliseconds>(cl_now - cl_start).count();
	if (restart) cl_start = cl_now;
	return dur;
}

long long perfclock::getsec(bool restart)
{
	std::chrono::time_point<std::chrono::steady_clock> cl_now = std::chrono::steady_clock::now();
	long long dur = std::chrono::duration_cast<std::chrono::seconds>(cl_now - cl_start).count();
	if (restart) cl_start = cl_now;
	return dur;
}
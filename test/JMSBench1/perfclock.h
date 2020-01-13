#pragma once
#include <time.h>
#include <chrono>
class perfclock
{
public:
	perfclock();
	~perfclock();
protected:
	std::chrono::time_point<std::chrono::steady_clock> cl_start;
public:
	std::chrono::time_point<std::chrono::steady_clock> start(void);
	long long getnsec(bool restart = false);
	long long getusec(bool restart = false);
	long long getmsec(bool restart = false);
	long long getsec(bool restart = false);
};


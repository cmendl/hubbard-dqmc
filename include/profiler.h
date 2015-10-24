#ifndef PROFILER_H
#define PROFILER_H


#ifdef PROFILE_ENABLE

#ifdef _WIN32

#include <windows.h>

#define TIME_TYPE LARGE_INTEGER
#define GET_TICKS(t) QueryPerformanceCounter((t))
#define DELTA_TICKS(a, b) ((a.QuadPart)-(b.QuadPart))

#else

#include <time.h>

#define TIME_TYPE struct timespec
#define GET_TICKS(t) clock_gettime(CLOCK_MONOTONIC, (t))
#define TICKS_PER_SEC 1000000000LL
#define DELTA_TICKS(a, b) \
	(TICKS_PER_SEC * ((a).tv_sec - (b).tv_sec) + (long long)((a).tv_nsec - (b).tv_nsec))

#endif


#define PROFILE_BEGIN(name) \
	TIME_TYPE name##_start_, name##_end_; \
	GET_TICKS(&name##_start_)

#define PROFILE_END(name) \
	GET_TICKS(&name##_end_); \
	Profile_Add(#name, DELTA_TICKS(name##_end_, name##_start_))


void Profile_Start(void);

void Profile_Add(const char *name, long long delta);

void Profile_Stop(void);

#else

#define PROFILE_BEGIN(n)
#define PROFILE_END(n)
#define Profile_Start()
#define Profile_Stop()

#endif



#endif

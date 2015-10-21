#ifndef PROFILER_H
#define PROFILER_H


#ifdef PROFILE_ENABLE

#ifdef _WIN32
	#error No profiling support for Windows yet.
#endif

#include <time.h>

#define CLOCK_TYPE CLOCK_PROCESS_CPUTIME_ID

#define DELTA_NS(a, b) \
	(1000000000LL * ((a).tv_sec - (b).tv_sec) + (long long)((a).tv_nsec - (b).tv_nsec))

#define PROFILE_BEGIN(name) \
	struct timespec name##_start_, name##_end_; \
	clock_gettime(CLOCK_TYPE, &name##_start_)

#define PROFILE_END(name) \
	clock_gettime(CLOCK_TYPE, &name##_end_); \
	Profile_Add(#name, DELTA_NS(name##_end_, name##_start_))


void Profile_Start(void);

void Profile_Add(const char *name, long long delta);

void Profile_Stop(void);

#else

#define PROFILE_BEGIN(n)
#define PROFILE_END(n)
#define Profile_Start()
#define Profile_Add(n, d)
#define Profile_Stop()

#endif



#endif

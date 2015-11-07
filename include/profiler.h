#ifndef PROFILER_H
#define PROFILER_H


#ifdef PROFILE_ENABLE

void Profile_Start(void);

void Profile_Begin(const char *name);

void Profile_End(const char *name);

void Profile_Stop(void);

#else

#define Profile_Start()
#define Profile_Begin(n)
#define Profile_End(n)
#define Profile_Stop()

#endif



#endif

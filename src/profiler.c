#ifdef PROFILE_ENABLE

#include "profiler.h"
#include "hash_table.h"
#include "dupio.h"
#include <mkl.h>
#include <omp.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

typedef struct
{
	long long start_tick;
	int calls;
	long long total;
}
profile_entry_t;

// hash table of all keys and profile entries.
static ht_t profile_table;

// for total wall time calculation
static long long main_start_tick;


//________________________________________________________________________________________________________________________
///
/// \brief Get current time tick
///
static inline long long get_ticks(void)
{
#ifdef _WIN32
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	return t.QuadPart;
#else
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return 1000000000LL * t.tv_sec + tv_nsec;
#endif
}


//________________________________________________________________________________________________________________________
///
/// \brief Call at beginning of main() to initialize profiling.
///
void Profile_Start(void)
{
	const int n_buckets = 32;
	htInit(&profile_table, n_buckets);
	main_start_tick = get_ticks();
}


//________________________________________________________________________________________________________________________
///
/// \brief Begin a profiling block
///
void Profile_Begin(const char *name)
{
	// append a -<thread number> to name if not called from main thread.
	char name2[128];
	int thread_num = omp_get_thread_num();
	if (thread_num > 0)
	{
		snprintf(name2, 128, "%s-%d", name, thread_num);
		name = name2;
	}

	#pragma omp critical(profiler)
	{
		// pointer to corresponding profile entry.
		profile_entry_t *p = (profile_entry_t *)htGet(&profile_table, name);

		// initializes everything the first time around
		if (p == NULL)
		{
			p = (profile_entry_t *)MKL_calloc(1, sizeof(profile_entry_t), MEM_DATA_ALIGN);
			htInsert(&profile_table, name, p);
		}

		// don't do anything if a profiling block with the same name is already running
		if (p->start_tick == 0)
		{
			p->calls++;
			p->start_tick = get_ticks();
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief End a profiling block
///
void Profile_End(const char *name)
{
	// append a -<thread number> to name if not called from main thread.
	char name2[128];
	int thread_num = omp_get_thread_num();
	if (thread_num > 0)
	{
		snprintf(name2, 128, "%s-%d", name, thread_num);
		name = name2;
	}

	#pragma omp critical(profiler)
	{
		// pointer to corresponding profile entry.
		profile_entry_t *p = (profile_entry_t *)htGet(&profile_table, name);

		// p should never be NULL unless Profile_End is called before Profile_Begin.
		if (p != NULL)
		{
			p->total += get_ticks() - p->start_tick;
			p->start_tick = 0; // reset to 0 to indicate completion of a profiling block
		}
	}
}


// sort by total time
/*
static int entry_compare(const void *a, const void *b)
{
	ht_entry_t *ea = *(ht_entry_t **)a;
	ht_entry_t *eb = *(ht_entry_t **)b;
	profile_entry_t *pa = (profile_entry_t *)ea->val;
	profile_entry_t *pb = (profile_entry_t *)eb->val;
	long long diff = pb->total - pa->total;
	return (diff > 0) ? 1 : (diff < 0) ? -1 : 0; // in case casting to int overflows.
}
*/

// sort by name
#include <string.h>
static int entry_compare(const void *a, const void *b)
{
	ht_entry_t *ea = *(ht_entry_t **)a;
	ht_entry_t *eb = *(ht_entry_t **)b;
	return strcmp(ea->key, eb->key);
}


//________________________________________________________________________________________________________________________
///
/// \brief Call at end of main() to stop profiling and print report.
///
void Profile_Stop(void)
{
	// get total wall time
	long long main_total = get_ticks() - main_start_tick;

	// get the tick resolution
#ifdef _WIN32
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	const double ticks_per_sec = (double)freq.QuadPart;
#else
	const double ticks_per_sec = (double)1000000000;
#endif

	// place (pointers to) every entry into a simple array, and sort with entry_compare.
	ht_entry_t **entries_table = (ht_entry_t **)MKL_malloc(profile_table.n_entries * sizeof(ht_entry_t *), MEM_DATA_ALIGN);
	int i, j = 0;
	for (i = 0; i < profile_table.n_buckets; i++)
	{
		ht_entry_t *entry;
		for (entry = profile_table.buckets[i]; entry != NULL; entry = entry->next)
			entries_table[j++] = entry;
	}
	qsort(entries_table, profile_table.n_entries, sizeof(ht_entry_t *), entry_compare);

	// print report
	duprintf("==============================Profiling report=================================\n");
	duprintf("%-31s%-8s%-10s%-12s%s\n", "name(-thread number)", "calls", "% of all", "total (s)", "time per call (us)");
	duprintf("%-49s%g\n","total wall time", main_total / ticks_per_sec);
	for (i = 0; i < profile_table.n_entries; i++)
	{
		profile_entry_t *p = (profile_entry_t *)entries_table[i]->val;
		duprintf("%-31s%-8d%-10g%-12g%-g\n",
				entries_table[i]->key,
				p->calls,
				(100. * p->total) / main_total,
				p->total / ticks_per_sec,
				p->total / (ticks_per_sec / 1000000. * p->calls));
	}
	duprintf("_______________________________________________________________________________\n");

	// free all allocated memory used for profiling.
	MKL_free(entries_table);
	htFree(&profile_table);
}



#endif

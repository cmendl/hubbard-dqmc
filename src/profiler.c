#ifdef PROFILE_ENABLE

#include "hash_table.h"
#include "dupio.h"
#include "profiler.h"
#include <mkl.h>
#include <omp.h>

// hash table of all keys and profile entries.
static ht_t profile_table;

// for total wall time calculation
static TIME_TYPE main_start_, main_end_;

typedef struct
{
	int calls;
	long long total;
}
profile_entry;


//________________________________________________________________________________________________________________________
///
/// \brief Call at beginning of main() to initialize profiling.
///
void Profile_Start(void)
{
	const int n_buckets = 32;
	htInit(&profile_table, n_buckets);
	GET_TICKS(&main_start_);
}


//________________________________________________________________________________________________________________________
///
/// \brief Adds delta ticks to the profile entry corresponding to name.
///
void Profile_Add(const char *name, long long delta)
{
	#pragma omp critical(profile_add)
	{
	// append a -<thread number> to name if Profile_Add is not called from main thread.
	char name2[128];
	int thread_num = omp_get_thread_num();
	if (thread_num > 0)
	{
		sprintf(name2, "%s-%d", name, thread_num);
		name = name2;
	}

	// pointer to corresponding profile entry.
	profile_entry *p = (profile_entry *)htGet(&profile_table, name);

	// initializes everything the first time around
	if (p == NULL)
	{
		p = (profile_entry *)MKL_malloc(sizeof(profile_entry), MEM_DATA_ALIGN);
		*p = (profile_entry) { 0 };
		htInsert(&profile_table, name, p);
	}

	// update profile entry
	p->calls++;
	p->total += delta;
	}
}


// sort by total time
/*
static int entry_compare(const void *a, const void *b)
{
	ht_entry_t *ea = *(ht_entry_t **)a;
	ht_entry_t *eb = *(ht_entry_t **)b;
	profile_entry *pa = (profile_entry *)ea->val;
	profile_entry *pb = (profile_entry *)eb->val;
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
	// stop total wall time counters
	GET_TICKS(&main_end_);
	long long main_total = DELTA_TICKS(main_end_, main_start_);

	// get the tick resolution
#ifdef _WIN32
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	const double ticks_per_sec = (double)freq.QuadPart;
#else
	const double ticks_per_sec = (double)TICKS_PER_SEC;
#endif

	// place (pointers to) every entries into a simple array, and sort with entry_compare.
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
	duprintf("===============================Profiling report=================================\n");
	duprintf("%-32s%-8s%-10s%-12s%s\n", "name(-thread number)", "calls", "% of all", "total (s)", "time per call (us)");
	duprintf("%-50s%g\n","total wall time", main_total / ticks_per_sec);
	for (i = 0; i < profile_table.n_entries; i++)
	{
		profile_entry *p = (profile_entry *)entries_table[i]->val;
		duprintf("%-32s%-8d%-10g%-12g%-g\n",
				entries_table[i]->key,
				p->calls,
				(100. * p->total) / main_total,
				p->total / ticks_per_sec,
				p->total / (ticks_per_sec / 1000000. * p->calls));
	}
	duprintf("--------------------------------------------------------------------------------\n");

	// free all allocated memory used for profiling.
	MKL_free(entries_table);
	htFree(&profile_table);
}



#endif

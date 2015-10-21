#ifdef PROFILE_ENABLE

#include "aa.h"
#include "dupio.h"
#include "profiler.h"
#include <mkl.h>
#include <stdlib.h>

#define PROFILE_MAX_ENTRIES 128

static aa_t profile_table;
static struct timespec main_start_, main_end_;

typedef struct
{
	int calls;
	long long total;
}
profile_entry;

void Profile_Start(void)
{
	aaInit(&profile_table, PROFILE_MAX_ENTRIES);
	clock_gettime(CLOCK_TYPE, &main_start_);
}

void Profile_Add(const char *name, long long delta)
{
	profile_entry **p;
	int thread_num = omp_get_thread_num();
	if (thread_num > 0)
	{
		char name2[128];
		sprintf(name2, "%s-%d", name, thread_num);
		p = (profile_entry **)aaGet(&profile_table, name2);
	}
	else
	{
		p = (profile_entry **)aaGet(&profile_table, name);
	}

	if (p == NULL)
		return;
	if (*p == NULL)
	{
		*p = (profile_entry *)MKL_malloc(sizeof(profile_entry), MEM_DATA_ALIGN);
		**p = (profile_entry){0};
	}
	(*p)->calls++;
	(*p)->total += delta;
}

/*
// sort by total time
static int entry_compare(const void *a, const void *b)
{
	int ia = *(int *)a;
	int ib = *(int *)b;
	profile_entry *pa = (profile_entry *)profile_table.nodes[ia].val;
	profile_entry *pb = (profile_entry *)profile_table.nodes[ib].val;
	long long diff = pb->total - pa->total;
	return (diff > 0) ? 1 : (diff < 0) ? -1 : 0; // in case casting to int overflows.
}
*/

// sort by name
#include <string.h>
static int entry_compare(const void *a, const void *b)
{
	int ia = *(int *)a;
	int ib = *(int *)b;
	char *na = profile_table.nodes[ia].key;
	char *nb = profile_table.nodes[ib].key;
	return strcmp(na, nb);
}


void Profile_Stop(void)
{

	clock_gettime(CLOCK_TYPE, &main_end_);
	long long main_total = DELTA_NS(main_end_, main_start_);

	int *sorted_i = (int *)MKL_malloc(profile_table.n * sizeof(int), MEM_DATA_ALIGN);
	int j;
	for (j = 0; j < profile_table.n; j++) {
		sorted_i[j] = j;
	}
	qsort(sorted_i, profile_table.n, sizeof(int), entry_compare);

	duprintf("===============================Profiling report=================================\n");
	duprintf("%-32s%-8s%-10s%-12s%s\n", "name(-thread number)", "calls", "% of all", "total (s)", "time per call (us)");
	duprintf("%-50s%g\n","total wall time", main_total/1000000000.);
	for (j = 0; j < profile_table.n; j++)
	{
		int i = sorted_i[j];
		profile_entry *p = (profile_entry *)profile_table.nodes[i].val;
		duprintf("%-32s%-8d%-10g%-12g%-g\n",
				profile_table.nodes[i].key,
				p->calls,
				(100. * p->total) / main_total,
				p->total / 1000000000.,
				p->total / (1000. * p->calls));
	}
	duprintf("--------------------------------------------------------------------------------\n");
	MKL_free(sorted_i);
	aaFree(&profile_table);
}



#endif

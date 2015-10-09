#include "aa.h"
#include <mkl.h>
#include <string.h>

void aaInit(aa_t *aa, int max)
{
	aa->nodes = (aa_node_t *)MKL_malloc(max * sizeof(aa_node_t), MEM_DATA_ALIGN);
	memset(aa->nodes, 0, max * sizeof(aa_node_t));
	aa->n = 0;
	aa->max = max;
}

void aaFree(aa_t *aa)
{
	int i;
	for (i = 0; i < aa->n; i++)
	{
		MKL_free(aa->nodes[i].key);
		MKL_free(aa->nodes[i].val);
	}
	MKL_free(aa->nodes);
}

// if key is found, returns index. if not found, returns -(1 + appropriate index for inserting new node).
static int aaFindIndex(aa_t *aa, const char *key)
{
	// binary search. O(log N).
	int mid, lo = 0, hi = aa->n;
	while (hi > lo)
	{
		mid = (lo + hi) / 2;
		int cmp = strcmp(key, aa->nodes[mid].key);
		if (cmp < 0)
			hi = mid;
		else if (cmp > 0)
			lo = mid + 1;
		else
			return mid; //found
	}
	return -(1 + lo); //not found
}

static int aaInsertIndex(aa_t *aa, int index, const char *key)
{
	if (index > aa->n || aa->n == aa->max)
		return 1;
	// shift everything above index up
	int i;
	for (i = aa->n; i > index; i--)
		aa->nodes[i] = aa->nodes[i - 1];
	aa->n++;
	// insert
	aa->nodes[index].key = (char *)MKL_malloc(strlen(key) + 1, MEM_DATA_ALIGN);
	strcpy(aa->nodes[index].key, key);
	aa->nodes[index].val = NULL;
	return 0;
}

static int aaDeleteIndex(aa_t *aa, int index)
{
	if (index >= aa->n || aa->n == 0)
		return 1;
	free(aa->nodes[index].key);
	free(aa->nodes[index].val);
	// shift everything above index down
	aa->n--;
	int i;
	for (i = index; i < aa->n; i++)
		aa->nodes[i] = aa->nodes[i + 1];
	aa->nodes[aa->n] = (aa_node_t){0};
	return 0;
}

// returns a pointer to val, which itself is a pointer to whatever
void **aaGet(aa_t *aa, const char *key)
{
	int index = aaFindIndex(aa, key);
	if (index < 0) // not found, insert new node
	{
		index = -index - 1;
		if (aaInsertIndex(aa, index, key) != 0)
			return NULL;
	}
	return &aa->nodes[index].val;
}

int aaDelete(aa_t *aa, const char *key)
{
	int index = aaFindIndex(aa, key);
	if (index < 0) // not found
		return 1;
	else
		return aaDeleteIndex(aa, index);
}

#ifndef AA_H
#define AA_H

// Associative array implemented with a sorted array and a binary search for lookup.
// All operations are O(N) except for lookup, which is O(log N).
typedef struct
{
	char *key;
	void *val;
}
aa_node_t;

typedef struct
{
	aa_node_t *nodes;
	int n; // number of nodes
	int max; // maximum possible number of nodes
}
aa_t;

void aaInit(aa_t *aa, int max);

void aaFree(aa_t *aa);

void **aaGet(aa_t *aa, const char *key);

int aaDelete(aa_t *aa, const char *key);



#endif

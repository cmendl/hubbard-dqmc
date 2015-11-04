#ifndef HASH_TABLE_H
#define HASH_TABLE_H

//associative array using a hash table and linked lists for collisions
//NOT THREAD-SAFE

typedef struct ht_entry_t
{
	struct ht_entry_t *next;
	char *key;
	void *val;
}
ht_entry_t;

typedef struct
{
	ht_entry_t **buckets; //each bucket is a linked list of entries
	int n_entries;
	int n_buckets; //number of buckets
}
ht_t;

void htInit(ht_t *ht, int n_buckets);

void htFree(ht_t *ht);

int htInsert(ht_t *ht, const char *key, void *val);

void *htGet(ht_t *ht, const char *key);

int htDelete(ht_t *ht, const char *key);



#endif

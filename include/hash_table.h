#ifndef HASH_TABLE_H
#define HASH_TABLE_H

//associative array using a hash table and linked lists for collisions
//NOT THREAD-SAFE

typedef struct ht_entry_s
{
	struct ht_entry_s *next;
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

void htInit(ht_t *ht, const int n_buckets);

void htFree(ht_t *ht);

void *htInsert(ht_t *ht, const char *key, void *val);

void *htGet(const ht_t *ht, const char *key);

int htDelete(ht_t *ht, const char *key);



#endif

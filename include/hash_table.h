#ifndef HASH_TABLE_H
#define HASH_TABLE_H


//________________________________________________________________________________________________________________________
///
/// \brief Hash table entry, forming a linked list
///
typedef struct ht_entry_s
{
	struct ht_entry_s *next;	//!< pointer to next entry
	char *key;					//!< key, as character string
	void *val;					//!< pointer to corresponding value
}
ht_entry_t;


//________________________________________________________________________________________________________________________
///
/// \brief Associative array using a hash function and linked lists for collisions;
/// not thread-safe
///
typedef struct
{
	ht_entry_t **buckets;	//!< each bucket is a linked list of entries
	int n_entries;			//!< number of entries
	int n_buckets;			//!< number of buckets
}
hash_table_t;


void AllocateHashTable(hash_table_t *ht, const int n_buckets);

typedef void FreeFuncPtr(void *ptr);
void DeleteHashTable(hash_table_t *ht, FreeFuncPtr *free_func);


void *HashTableInsert(hash_table_t *ht, const char *key, void *val);

void *HashTableGet(const hash_table_t *ht, const char *key);

void *HashTableRemove(hash_table_t *ht, const char *key);



#endif

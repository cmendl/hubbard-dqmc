#include "hash_table.h"
#include "util.h"
#include <stdint.h>
#include <string.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for a hash table
///
void AllocateHashTable(hash_table_t *ht, const int n_buckets)
{
	ht->n_buckets = n_buckets;
	ht->n_entries = 0;

	ht->buckets = (ht_entry_t **)algn_calloc(n_buckets, sizeof(ht_entry_t *));
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory for everything in a hash table
///
void DeleteHashTable(hash_table_t *ht, FreeFuncPtr *free_func)
{
	int i;
	for (i = 0; i < ht->n_buckets; i++)
	{
		ht_entry_t *entry = ht->buckets[i];
		while (entry != NULL)
		{
			algn_free(entry->key);
			free_func(entry->val);
			ht_entry_t *next = entry->next;
			algn_free(entry);
			entry = next;
		}
	}

	algn_free(ht->buckets);
}


//________________________________________________________________________________________________________________________
///
/// \brief Public domain FNV-1a (32-bit) hash algorithm from http://www.isthe.com/chongo/tech/comp/fnv/index.html
///
static uint32_t FNV1a(const char *key)
{
	uint32_t hash = 2166136261;
	while (*key)
	{
		hash = (hash ^ (*key)) * 16777619;
		key++;
	}
	return hash;
}


//________________________________________________________________________________________________________________________
///
/// \brief Insert an entry into the hash table; if key already exists,
/// replace previous value and return pointer to old value, otherwise return NULL
///
void *HashTableInsert(hash_table_t *ht, const char *key, void *val)
{
	// compute hash value of 'key'
	const int i = FNV1a(key) % ht->n_buckets;

	ht_entry_t **p_entry = &ht->buckets[i];     // p_entry is whatever was pointing at entry
	for (; (*p_entry) != NULL; p_entry = &(*p_entry)->next)
	{
		// current entry
		ht_entry_t *entry = *p_entry;

		if (strcmp(key, entry->key) == 0) // entry already exists!
		{
			void *ret = entry->val;
			entry->val = val;
			return ret;
		}
	}

	// key not found, insert
	ht_entry_t *entry = (ht_entry_t *)algn_malloc(sizeof(ht_entry_t));
	// copy key string
	entry->key = (char *)algn_malloc(strlen(key) + 1);
	strcpy(entry->key, key);
	entry->val = val;
	entry->next = NULL;
	*p_entry = entry;
	ht->n_entries++;

	return NULL;
}


//________________________________________________________________________________________________________________________
///
/// \brief Return a pointer to the value corresponding to 'key'; if the key is not found, return NULL
///
void *HashTableGet(const hash_table_t *ht, const char *key)
{
	// compute hash value of 'key'
	const int i = FNV1a(key) % ht->n_buckets;

	ht_entry_t *entry;
	for (entry = ht->buckets[i]; entry != NULL; entry = entry->next)
	{
		if (strcmp(key, entry->key) == 0) {
			return entry->val;
		}
	}

	// not found
	return NULL;
}


//________________________________________________________________________________________________________________________
///
/// \brief Remove entry with given key from hash table and return corresponding value; if the key cannot be found, return NULL
///
void *HashTableRemove(hash_table_t *ht, const char *key)
{
	// compute hash value of 'key'
	const int i = FNV1a(key) % ht->n_buckets;

	ht_entry_t **p_entry = &ht->buckets[i];     // p_entry is whatever was pointing at current entry
	for (; (*p_entry) != NULL; p_entry = &(*p_entry)->next)
	{
		// current entry
		ht_entry_t *entry = *p_entry;

		if (strcmp(key, entry->key) == 0)
		{
			// key found

			(*p_entry) = entry->next;       // redirect pointer, effectively removing current entry from linked list
			void *val = entry->val;         // keep reference to current value
			algn_free(entry->key);          // delete current entry
			algn_free(entry);
			ht->n_entries--;

			return val;
		}
	}

	// not found
	return NULL;
}

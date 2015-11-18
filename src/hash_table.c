#include "hash_table.h"
#include <mkl.h>
#include <stdint.h>
#include <string.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocates memory for a hash table
///
void htInit(ht_t *ht, const int n_buckets)
{
	ht->buckets = (ht_entry_t **)MKL_malloc(n_buckets * sizeof(ht_entry_t *), MEM_DATA_ALIGN);
	memset(ht->buckets, 0, n_buckets * sizeof(ht_entry_t *));
	ht->n_entries = 0;
	ht->n_buckets = n_buckets;
}


//________________________________________________________________________________________________________________________
///
/// \brief Frees memory for everything in a hash table
///
void htFree(ht_t *ht)
{
	int i;
	for (i = 0; i < ht->n_buckets; i++)
	{
		ht_entry_t *entry, *temp = ht->buckets[i];
		while ((entry = temp) != NULL)
		{
			MKL_free(entry->key);
			MKL_free(entry->val); //assumes this is sufficient for deallocating val
			temp = entry->next;
			MKL_free(entry);
		}
	}
	MKL_free(ht->buckets);
}


//________________________________________________________________________________________________________________________
///
/// \brief Public domain FNV-1a (32-bit) hash algorithm from http://www.isthe.com/chongo/tech/comp/fnv/index.html
///
static uint32_t fnv_1a(const char *key)
{
	uint32_t hash = 2166136261;
	while (*key)
		hash = (hash ^ *key++) * 16777619;
	return hash;
}


//________________________________________________________________________________________________________________________
///
/// \brief Inserts an entry into hash table; if key already exists,
// replace previous value and return pointer to old value, otherwise return NULL
///
void *htInsert(ht_t *ht, const char *key, void *val)
{
	const int i = fnv_1a(key) % ht->n_buckets;
	ht_entry_t **p_entry, *entry; //p_entry is whatever was pointing at entry.
	for (p_entry = ht->buckets + i; (entry = *p_entry) != NULL; p_entry = &(entry->next))
	{
		if (!strcmp(key, entry->key)) // entry already exists!
		{
			void *ret = entry->val;
			entry->val = val;
			return ret;
		}
	}
	//not found, insert.
	entry = (ht_entry_t *)MKL_malloc(sizeof(ht_entry_t), MEM_DATA_ALIGN);
	entry->key = (char *)MKL_malloc(strlen(key) + 1, MEM_DATA_ALIGN);
	strcpy(entry->key, key);
	entry->val = val;
	entry->next = NULL;
	*p_entry = entry;
	ht->n_entries++;

	return NULL;
}


//________________________________________________________________________________________________________________________
///
/// \brief Returns val, which is a pointer to the actual value. If key not found, return NULL.
///
void *htGet(const ht_t *ht, const char *key)
{
	const int i = fnv_1a(key) % ht->n_buckets;
	ht_entry_t *entry;
	for (entry = ht->buckets[i]; entry != NULL; entry = entry->next)
	{
		if (!strcmp(key, entry->key))
			return entry->val;
	}
	//not found
	return NULL;
}

//________________________________________________________________________________________________________________________
///
/// \brief Removes and deallocates entry from hash table. If key not found, return 1.
///
int htDelete(ht_t *ht, const char *key)
{
	const int i = fnv_1a(key) % ht->n_buckets;
	ht_entry_t **p_entry, *entry; //p_entry is whatever was pointing at entry.
	for (p_entry = ht->buckets + i; (entry = *p_entry) != NULL; p_entry = &(entry->next))
	{
		if (!strcmp(key, entry->key))
		{
			*p_entry = entry->next;
			MKL_free(entry->key);
			MKL_free(entry->val);
			MKL_free(entry);
			ht->n_entries--;
			return 0;
		}
	}
	//not found
	return 1;
}

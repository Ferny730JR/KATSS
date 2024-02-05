// hash_table.h

#ifndef BPP_HASH_TABLE_H
#define BPP_HASH_TABLE_H

#include <stddef.h>

typedef struct {
    double *values;
} Data;

typedef struct {
    char *key;
    Data data;
} Entry;

typedef struct {
    size_t size;
    Entry *entries;
} HashTable;

unsigned int hash(const char *key);

HashTable init_hash_table(int kmer);

void free_hash_table(HashTable *hash_table);

double *get(HashTable *hash_table, const char *key);

#endif  // BPP_HASH_TABLE_H
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
} bppHashTable;

unsigned int hash(const char *key);

bppHashTable init_hash_table(int kmer);

void free_hash_table(bppHashTable *hash_table);

double *get(bppHashTable *hash_table, const char *key);

void addValue(bppHashTable *hash_table, const char *key, double value, int value_index);

void printBPPHashTable(bppHashTable *hash_table, int kmer);

#endif  // BPP_HASH_TABLE_H
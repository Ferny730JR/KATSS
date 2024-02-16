#ifndef KMER_HASH_TABLE_H
#define KMER_HASH_TABLE_H

#include <stddef.h>

/**
 *  @brief Entries for kmerHashTable data structure. Stores the key value pair, where key is a 
 *  kmer and the value is a double array
*/
typedef struct {
    char *key;
    double *values;
} Entry;


/**
 *  @brief Hash table data structure to store kmer's and associated information.
*/
typedef struct {
    size_t  capacity;
    size_t  cols;
    Entry   **entries;
} kmerHashTable;


kmerHashTable *init_kmer_table(int kmer, int cols);


kmerHashTable *init_bpp_table(int kmer);


void free_kmer_table(kmerHashTable *hash_table);


double *kmer_get(kmerHashTable  *hash_table, 
                 const char     *key);


void kmer_add_value(kmerHashTable   *hash_table, 
                    const char      *key, 
                    double          value, 
                    int             value_index);


void print_kmer_table(kmerHashTable *hash_table);


#endif  // BPP_HASH_TABLE_H
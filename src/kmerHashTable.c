#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "kmerHashTable.h"
#include "utils.h"
#include "string_utils.h"

Entry *create_entry(const char  *key, 
                    int         col);

unsigned int hash(const char    *key);

void free_entry(Entry *entry);


kmerHashTable *init_kmer_table(int kmer, int cols) {
    size_t table_size = 1 << (2 * kmer); // 4^kmer
    kmerHashTable *hash_table = s_malloc(sizeof *hash_table);

    hash_table->capacity = table_size;
    hash_table->cols = cols;
    hash_table->entries = s_malloc(table_size * sizeof(Entry*));
    pthread_mutex_init(&hash_table->lock, NULL);

    for (size_t i = 0; i < table_size; i++) {
        hash_table->entries[i] = NULL;
    }

    return hash_table;
}


kmerHashTable *init_bpp_table(int kmer) {
    return init_kmer_table(kmer, kmer+1);
}


double *kmer_get(kmerHashTable *hash_table, const char *key) {
    unsigned int index = hash(key);
    if (strcmp(hash_table->entries[index]->key, key) == 0) {
        return hash_table->entries[index]->values;
    }
    error_message("Unable to get values from table.\n"
                  "Expected key: '%s' - Actual key: '%s'",key, hash_table->entries[index]->key);
    exit(EXIT_FAILURE);
}


void kmer_add_value(kmerHashTable *hash_table, const char *key, double value, int value_index) {
    if(value_index >= hash_table->cols) {
        error_message("value_index '%d' is greater than length of value array, which is '%d'.",
        value_index, hash_table->cols);
        exit(EXIT_FAILURE);
    }

    unsigned int index = hash(key);

    if(hash_table->entries[index] == NULL) { // make new entry if not initialized
        Entry *new_item = create_entry(key, hash_table->cols);
        hash_table->entries[index] = new_item;
    }

    if (strcmp(hash_table->entries[index]->key, key) == 0) {
        pthread_mutex_lock(&hash_table->lock);
        hash_table->entries[index]->values[value_index] += value;
        pthread_mutex_unlock(&hash_table->lock);
        return;
    }

    error_message("Unable to add value to table.\n"
                  "Expected key: '%s' - Actual key: '%s'",key, hash_table->entries[index]->key);
    exit(EXIT_FAILURE);
}


Entry *create_entry(const char *key, int col) {
    Entry *entry = s_malloc(sizeof *entry);
    entry->key = strdup(key);
    entry->values = s_calloc(col, sizeof(double));

    return entry;
}


unsigned int hash(const char *key) {
    // Assuming the key is composed of 'A', 'U/T', 'C', 'G' characters
    unsigned int hash_value = 0;
    while (*key) {
        switch(*key) {
            case 'A': hash_value = hash_value * 4;     break;
            case 'C': hash_value = hash_value * 4 + 1; break;
            case 'G': hash_value = hash_value * 4 + 2; break;
            default : hash_value = hash_value * 4 + 3; break;
        }
        key++;
    }
    return hash_value;
}


void free_entry(Entry *entry) {
    free(entry->key);
    free(entry->values);
    free(entry);
}


void free_kmer_table(kmerHashTable *hash_table) {
    for (size_t i = 0; i < hash_table->capacity; i++) {
        if(hash_table->entries[i]) {
            free_entry(hash_table->entries[i]);
        }
    }
    pthread_mutex_destroy(&hash_table->lock);
    free(hash_table->entries);
    free(hash_table);
}


void print_table_to_file(kmerHashTable *table, FILE *table_file, char sep) {
    for(size_t i = 0; i < table->capacity; i++) {
        if(table->entries[i] == NULL) {
            continue;
        }

        fprintf(table_file, "%s", table->entries[i]->key);
        for(size_t j = 0; j < table->cols; j++) {
            fprintf(table_file, "%c%9.6f", sep, table->entries[i]->values[j]);
        }
        fprintf(table_file,"\n");
    }
}


void kmerHashTable_to_file(kmerHashTable *table, char *name, char file_delimiter) {
    char *filename = concat(name, ".dsv");
    
    FILE *table_file = fopen(filename, "w");
    if (table_file == NULL) {
        error_message("Could not write to file '%s'\n",filename);
    }
    
    print_table_to_file(table, table_file, file_delimiter);

    free(filename);
    fclose(table_file);
}


void print_kmer_table(kmerHashTable *hash_table) {
    printf("--- BEGIN KMER TABLE ---\n");
    for(size_t i = 0; i < hash_table->capacity; i++) {
        if(hash_table->entries[i] == NULL) {
            continue;
        }
        
        printf("Key: %s, Values: ", hash_table->entries[i]->key);
        for(int j = 0; j < hash_table->cols; j++) {
            printf("%f ",hash_table->entries[i]->values[j]);
        }
        printf("\n");
    }
    printf("---- END KMER TABLE ----\n");
}
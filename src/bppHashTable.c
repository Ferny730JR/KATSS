#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bppHashTable.h"

#define BASES "ACGU"

unsigned int hash(const char *key) {
    // Assuming the key is composed of 'A', 'U', 'C', 'G' characters
    unsigned int hash_value = 0;
    while (*key) {
        hash_value = hash_value * 4 + (*key == 'A' ? 0 : (*key == 'C' ? 1 : (*key == 'G' ? 2 : 3)));
        key++;
    }
    return hash_value;
}

bppHashTable *init_hash_table(int kmer) {
    size_t table_size = 1 << (2 * kmer); // 4^kmer
    bppHashTable *hash_table = (bppHashTable*)malloc(sizeof(bppHashTable));
    hash_table->size = table_size;
    hash_table->entries = malloc(table_size * sizeof(Entry));

    if (!hash_table->entries) {
        perror("Failed to allocate memory for hash table entries");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < table_size; i++) {
        hash_table->entries[i] = NULL;
    }

    return hash_table;
}

Entry *create_entry(const char *key) {
    Entry *entry = (Entry*)malloc(sizeof(Entry));
    entry->key = strdup(key);
    entry->values = calloc((strlen(key)+1), sizeof(float));

    if(!entry->values) {
        perror("Failed to allocate memory for values array");
        exit(EXIT_FAILURE);
    }

    return entry;
}

void free_entry(Entry *entry) {
    free(entry->key);
    free(entry->values);
    free(entry);
}

void free_hash_table(bppHashTable *hash_table) {
    for (size_t i = 0; i < hash_table->size; i++) {
        if(hash_table->entries[i]) {
            free_entry(hash_table->entries[i]);
        }
    }
    free(hash_table->entries);
    free(hash_table);
}

float *get(bppHashTable *hash_table, const char *key) {
    unsigned int index = hash(key);
    if (strcmp(hash_table->entries[index]->key, key) == 0) {
        return hash_table->entries[index]->values;
    }
    exit(EXIT_FAILURE);
}

void addValue(bppHashTable *hash_table, const char *key, float value, int value_index) {
    unsigned int index = hash(key);

    if(hash_table->entries[index] == NULL) { // make new entry if not initialized
        Entry *new_item = create_entry(key);
        hash_table->entries[index] = new_item;
    }

    if (strcmp(hash_table->entries[index]->key, key) == 0) {
        hash_table->entries[index]->values[value_index] += value;
        return;
    }
    exit(EXIT_FAILURE);
}

void printBPPHashTable(bppHashTable *hash_table, int kmer) {
    printf("--- BEGIN HASH TABLE ---\n");
    for(size_t i = 0; i < hash_table->size; i++) {
        if(hash_table->entries[i] == NULL) {
            continue;
        }
        
        printf("Key: %s, Values: ", hash_table->entries[i]->key);
        for(int j = 0; j < kmer+1; j++) {
            printf("%f ",hash_table->entries[i]->values[j]);
        }
        printf("\n");
    }
    printf("---- END HASH TABLE ----\n");
}
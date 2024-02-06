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

bppHashTable init_hash_table(int kmer) {
    size_t table_size = 1 << (2 * kmer); // 4^kmer
    bppHashTable hash_table;
    hash_table.size = table_size;
    hash_table.entries = malloc(table_size * sizeof(Entry));

    if (!hash_table.entries) {
        perror("Failed to allocate memory for hash table");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < table_size; i++) {
        hash_table.entries[i].key = NULL;
        hash_table.entries[i].data.values = malloc(kmer * sizeof(double));

        if (!hash_table.entries[i].data.values) {
            perror("Failed to allocate memory for values array");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < kmer; j++) {
            hash_table.entries[i].data.values[j] = 0.0;
        }
    }

    //int key_length = 1 << (2 * kmer);
    for (int i = 0; i < table_size; i++) {
        hash_table.entries[i].key = malloc((kmer + 1) * sizeof(char));
        if (!hash_table.entries[i].key) {
            perror("Failed to allocate memory for key");
            exit(EXIT_FAILURE);
        }

        int index = i;
        for (int j = kmer - 1; j >= 0; j--) {
            hash_table.entries[i].key[j] = BASES[index % 4];
            index /= 4;
        }
        hash_table.entries[i].key[kmer] = '\0';
    }

    return hash_table;
}

void free_hash_table(bppHashTable *hash_table) {
    for (size_t i = 0; i < hash_table->size; i++) {
        free(hash_table->entries[i].key);
        free(hash_table->entries[i].data.values);
    }
    free(hash_table->entries);
}

double *get(bppHashTable *hash_table, const char *key) {
    unsigned int index = hash(key);
    while (hash_table->entries[index].key != NULL) {
        if (strcmp(hash_table->entries[index].key, key) == 0) {
            return hash_table->entries[index].data.values;
        }
        index = (index + 1) % hash_table->size; // Linear probing for collisions
    }
    return NULL; // Key not found
}

void printBPPHashTable(bppHashTable hash_table, int kmer) {
    for(size_t i = 0; i < hash_table.size; i++) {
        printf("Key: %s, Values: ", hash_table.entries[i].key);
        for(int j = 0; j < kmer; j++)
            printf("%f ",hash_table.entries[i].data.values[j]);
        printf("\n");
    }
}
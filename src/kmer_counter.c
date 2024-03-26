#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_counter.h"
#include "utils.h"

KmerCounter *
init_kcounter(unsigned int k_mer) 
{
	KmerCounter *kctr = s_malloc(sizeof *kctr);
	kctr->k_mer = k_mer;
	kctr->capacity = 1 << (2 * k_mer);   // 4^k_mer
	kctr->entries = s_calloc(kctr->capacity, sizeof(unsigned int));
	kctr->total_count = 0;

	return kctr;
}


void
free_kcounter(KmerCounter *kmer_counter)
{
	free(kmer_counter->entries);
	free(kmer_counter);
}


long
kctr_hash(const char *key, int start, int length) 
{
	unsigned int hash_value = 0;
	for(int i=start; i<start+length; i++) {
		switch(key[i]) {
			case 'A': hash_value = hash_value * 4;     break;
			case 'C': hash_value = hash_value * 4 + 1; break;
			case 'G': hash_value = hash_value * 4 + 2; break;
			case 'T': hash_value = hash_value * 4 + 3; break;
			case 'U': hash_value = hash_value * 4 + 3; break;
			default : return -1;
		}
	}

	return hash_value;
}


unsigned int
kctr_hash_key(const char *key)
{
	return kctr_hash(key, 0, (unsigned int)strlen(key));
}


char* 
kctr_unhash(unsigned int hash_value, int length) 
{
	char* key = s_malloc((length + 1) * sizeof(char));
	key[length] = '\0'; // Null-terminate the string

	for (int i = length - 1; i >= 0; i--) {
		switch (hash_value % 4) {
			case 0: key[i] = 'A'; break;
			case 1: key[i] = 'C'; break;
			case 2: key[i] = 'G'; break;
			case 3: key[i] = 'U'; break; // Assuming 'U' for the default case
		}	// todo: make it so case 3 can be either 'U' or 'T' based on specifications
		hash_value /= 4;
	}

	return key;
}


void 
kctr_increment(KmerCounter *kcounter, char *sequence)
{
	int num_kmers_in_seq = strlen(sequence) - kcounter->k_mer + 1;

	long index;
	for(int i=0; i<num_kmers_in_seq; i++) {
		index = kctr_hash(sequence,i,kcounter->k_mer);
		if(index == -1) {
			continue;
		}
		kcounter->entries[index]++;
	}
}


unsigned long
kctr_total(KmerCounter *kcounter)
{
	if(kcounter->total_count) {
		return kcounter->total_count;
	}

	unsigned long total = 0;
	for(unsigned int i=0; i<kcounter->capacity; i++) {
		total += kcounter->entries[i];
	}

	kcounter->total_count = total;
	return total;
}


void 
kctr_increment_substr(KmerCounter *kmer_counter, char *sequence,
                      unsigned int start, unsigned int length) 
{
	kmer_counter->entries[kctr_hash(sequence, start, length)]++;
}


unsigned int 
kctr_get(KmerCounter *kmer_counter, char *key) 
{
	unsigned int index = kctr_hash_key(key);
	return kmer_counter->entries[index];
}


char *
kctr_get_key(KmerCounter *kmer_counter, unsigned int index) 
{
	return kctr_unhash(index, kmer_counter->k_mer);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kmer_counter.h"
#include "string_utils.h"
#include "utils.h"

typedef struct {
	int shift;
	int start;
} DecrementValues;


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
			case 3: key[i] = 'T'; break; // Assuming 'T' for the default case
		}	// todo: make it so case 3 can be either 'U' or 'T' based on specifications
		hash_value /= 4;
	}

	return key;
}


void 
kctr_increment(KmerCounter *kcounter, const char *sequence)
{
	int num_kmers_in_seq = strlen(sequence) - kcounter->k_mer + 1;

	long index;
	for(int i=0; i<num_kmers_in_seq; i++) {
		index = kctr_hash(sequence,i,kcounter->k_mer);
		if(index == -1) {
			continue;
		}
		kcounter->entries[index]++;
		kcounter->total_count++;
	}
}


static DecrementValues
decrement_kmer(KmerCounter *kcounter, const char *sequence, 
               const char *pat, int min_start, int max_end)
{
	DecrementValues vals = {.shift = max_end, .start = 0};

	int pat_indx = subindx(sequence, pat);
	if(pat_indx == -1) {
		return vals;
	}

	int end = pat_indx + kcounter->k_mer;
	end = end > max_end ? max_end : end;

	int start = pat_indx - kcounter->k_mer + 1;
	start = start < min_start ? min_start : pat_indx - kcounter->k_mer + 1; 

	long index;
	while(start < end) {
		index = kctr_hash(sequence, start, kcounter->k_mer);
		if(index == -1) { start++; continue; }
		kcounter->entries[index]--;
		kcounter->total_count--;

		start++;
	}

	vals.shift = pat_indx + kcounter->k_mer;
	vals.start = end - vals.shift;

	return vals;
}


void
kctr_decrement(KmerCounter *kcounter, const char *sequence, const char *pat)
{
	int shift = 0;
	int num_kmers_in_seq = strlen(sequence) - kcounter->k_mer + 1;
	DecrementValues vals = {.shift = 0, .start = 0};

	while(shift < num_kmers_in_seq) {
		vals = decrement_kmer(kcounter, sequence+shift, pat, vals.start, num_kmers_in_seq-shift);
		shift += vals.shift;
	}
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


void
kctr_empty(KmerCounter *kmer_counter, char *key)
{
	unsigned int index = kctr_hash_key(key);
	kmer_counter->entries[index] = 0;
}


char *
kctr_get_key(KmerCounter *kmer_counter, unsigned int index) 
{
	return kctr_unhash(index, kmer_counter->k_mer);
}

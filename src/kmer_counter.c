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

/* Function declarations */
static void initialize_mapping(KmerCounter *kcounter);

KmerCounter *
init_kcounter(unsigned int k_mer) 
{
	KmerCounter *kctr = s_malloc(sizeof *kctr);
	kctr->k_mer = k_mer;
	kctr->capacity = 1 << (2 * k_mer);   // 4^k_mer
	kctr->entries = s_calloc(kctr->capacity, sizeof(unsigned int));
	kctr->total_count = 0;
	kctr->is_t = 1; // assume file is T based by default
	initialize_mapping(kctr);

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


void
kctr_unhash(char *key, unsigned int hash_value, int length, int is_t) 
{
	// char* key = s_malloc((length + 1) * sizeof(char));
	key[length] = '\0'; // Null-terminate the string

	for (int i = length - 1; i >= 0; i--) {
		switch (hash_value % 4) {
			case 0: key[i] = 'A'; break;
			case 1: key[i] = 'C'; break;
			case 2: key[i] = 'G'; break;
			case 3: key[i] = is_t ? 'T' : 'U'; break;
		}
		hash_value /= 4;
	}

	// return key;
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


void
kctr_fincrement(KmerCounter *kcounter, const char *sequence)
{
	int sequence_length = strlen(sequence);
	int k = kcounter->k_mer;
	unsigned int word_mask = (1 << 2*k) - 1;

	unsigned int hash_value = kctr_hash(sequence, 0, k);
	kcounter->entries[hash_value]++;
	kcounter->total_count++;

	for (int i = k; i < sequence_length; i++) {
		unsigned int x = kcounter->nucleotide_to_number[(unsigned char)sequence[i]];
		hash_value = ((hash_value << 2) | x) & word_mask;
        
		kcounter->entries[hash_value]++;
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
	start = start < min_start ? min_start : start; 

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


void
kctr_get_key(KmerCounter *kmer_counter, char *key_ptr, unsigned int index) 
{
	kctr_unhash(key_ptr, index, kmer_counter->k_mer, kmer_counter->is_t);
}


void
kctr_set_t_or_u(KmerCounter *kmer_counter, int is_t)
{
	kmer_counter->is_t = is_t;
}


// Function to initialize the nucleotide_to_number array
static void
initialize_mapping(KmerCounter *kcounter)
{
	// Initialize all elements to 0
	for (int i = 0; i < 256; i++) {
		kcounter->nucleotide_to_number[i] = (unsigned int)0;
	}

	// Assign values for nucleotides
	kcounter->nucleotide_to_number['A'] = (unsigned int)0;
	kcounter->nucleotide_to_number['C'] = (unsigned int)1;
	kcounter->nucleotide_to_number['G'] = (unsigned int)2;
	kcounter->nucleotide_to_number['T'] = (unsigned int)3;
	kcounter->nucleotide_to_number['U'] = (unsigned int)3;

	// Handling for lowercase nucleotides
	kcounter->nucleotide_to_number['a'] = (unsigned int)0;
	kcounter->nucleotide_to_number['c'] = (unsigned int)1;
	kcounter->nucleotide_to_number['g'] = (unsigned int)2;
	kcounter->nucleotide_to_number['t'] = (unsigned int)3;
	kcounter->nucleotide_to_number['u'] = (unsigned int)3;
}

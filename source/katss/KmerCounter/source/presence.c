#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "seqfile.h"
#include "katss_core.h"

#define BUFFER_SIZE (1024*1024) // 1 megabyte, aka 1 million NT limit

static inline bool
seen_test_and_set(uint8_t *seen, uint32_t hash)
{
    uint8_t *byte = &seen[hash >> 3];
    uint8_t mask = (uint8_t)(UINT8_C(1) << (hash & 7U));

    bool already_seen = (*byte & mask) != 0;
    *byte |= mask;

    return already_seen;
}


KatssCounter *
katss_count_presence(const char *filename, int kmer)
{
	warning_message("in count overlaps...");

	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		return NULL;
	}

	/* Open file for reading */
	SeqFile file = seqfopen(filename, "r");
	if(file == NULL) {
		katss_free_counter(counter);
		error_message("%s", seqfstrerror(seqferrno));
		return NULL;
	}

	/* Open hasher */
	KatssHasher *hasher = katss_init_hasher(kmer, 's');

	/* Allocate seen set */
	size_t seen_size = (((size_t)counter->capacity) + 8U) / 8U;
	uint8_t *seen = s_malloc(seen_size * sizeof *seen);

	/* Set variables */
	char *buffer = s_malloc(BUFFER_SIZE * sizeof *buffer);
	uint32_t hash_value;
	uint32_t seqnum = 0;
	bool dingaling;

	/* Begin counting overlaps */
	while(seqfgets_unlocked(file, buffer, BUFFER_SIZE)) {
		/* Reset seen set */
		memset(seen, 0, seen_size);
		dingaling = false;

		/* Count overlaps */
		hasher->has_previous = false; // don't use previous sequence context
		katss_set_seq(hasher, buffer);
		while(katss_get_fh(hasher, &hash_value)) {
			if(seen_test_and_set(seen, hash_value))
				continue;
			if(hash_value == 0 && dingaling) {
				printf("Line: %d -- seing k-mer \"AAAAA\"", seqnum);
			} else if(hash_value == 0) {
				printf("First: %d\n", seqnum);
				dingaling = true;
			}
			katss_increment(counter, hash_value);
		}
		seqnum++;
	}

	free(seen);
	free(buffer);
	free(hasher);
	seqfclose(file);
	return counter;
}

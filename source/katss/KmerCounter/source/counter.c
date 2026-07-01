#ifdef _WIN32
#define _CRT_RAND_S
#define rand_r rand_s
#endif
#include <stdlib.h> // for random functions 
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) && !defined(__STDC_NO_THREADS__)
#  include <threads.h>
#else
#  include <tinycthread.h>
#endif

#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"
#include "ushuffle.h"
#include "seqfile.h"
#include "thread_safe_rand.h"
#define BUFFER_SIZE 65536U

struct threadinfo {
	SeqFile seqfile;
	KatssCounter *counter;
	unsigned int kmer;
	int sample;
	unsigned int *seed;
	char filetype;
};
typedef struct threadinfo threadinfo;

/*============ Counting Function Declarations ============*/
static int
count_file_mt(void *arg);

/*============= Actual Functions Declarations =============*/
KatssCounter *
katss_count_kmers(const char *filename, unsigned int kmer)
{
	KatssCounter *counter = NULL;

	/* Open file for reading, with auto detection mode (for hasher) */
	SeqFile read_file = seqfopen(filename, "r");
	if(read_file == NULL) {
		error_message("%s", seqfstrerror(seqferrno));
		goto exit;
	}

	KatssHasher *hasher = katss_init_hasher(kmer, seqftype(read_file));
	if(hasher == NULL)
		goto cleanup_file;

	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* Prepare file reading & hash int */
	seqfsettype(read_file, 'b'); // set to binary reading
	char buffer[BUFFER_SIZE+1] = { 0 };
	size_t still_reading;
	uint32_t hash_value;

	do {
		still_reading = seqfread_unlocked(read_file, buffer, BUFFER_SIZE);
		buffer[still_reading] = '\0';

		katss_set_seq(hasher, buffer);
		while(katss_get_fh(hasher, &hash_value)) {
			katss_increment(counter, hash_value);
		}
	} while(still_reading == BUFFER_SIZE);

	/* If error was encountered while reading report and return NULL */
	if(still_reading == 0 && seqferrno) {
		katss_free_counter(counter);
		counter = NULL;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	return counter;
}


KatssCounter *
katss_count_kmers_mt(const char *filename, unsigned int kmer, int threads)
{
	/* Threads should be at least one, and at most 128 */
	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* If one thread, use single threaded computation */
	if(threads == 1)
		return katss_count_kmers(filename, kmer);

	/* Open SeqFile for reading */
	SeqFile file = seqfopen(filename, "r");
	if(file == NULL) {
		error_message("%s", seqfstrerror(seqferrno));
		return NULL;
	}

	/* Initialize counter */
	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		seqfclose(file);
		return NULL;
	}

	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	for(int i=0; i<threads; i++) {
		jobarg[i].seqfile = file;
		jobarg[i].counter = counter;
		jobarg[i].kmer = kmer;
		jobarg[i].filetype = seqftype(file);

		/* Start threads */
		thrd_create(&jobs[i], count_file_mt, &jobarg[i]);
	}

	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], NULL);
	}

	/* Free resources */
	seqfclose(file);
	free(jobs);
	free(jobarg);

	return counter;
}


KatssCounter *
katss_count_kmers_bootstrap(const char *filename, unsigned int kmer,
                            int sample, unsigned int *seed)
{
	KatssCounter *counter = NULL;

	/* Open file and hasher */
	SeqFile read_file = seqfopen(filename, "r");
	if(read_file == NULL) {
		error_message("%s", seqfstrerror(seqferrno));
		goto exit;
	}

	/* Hasher in sequences mode since we are using `seqfgets` */
	KatssHasher *hasher = katss_init_hasher(kmer, 's');
	if(hasher == NULL)
		goto cleanup_file;

	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	unsigned int local_seed;
	if(seed == NULL) {
		local_seed = time(NULL);
		seed = &local_seed;
	}

	/* Initialize buffer */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	uint32_t hash_value;
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		if(rand_r(seed) % 100000 >= sample)
			continue;
		katss_set_seq(hasher, buffer);
		while(katss_get_fh(hasher, &hash_value)) {
			katss_increment(counter, hash_value);
		}
	}

	if(seqferrno) {
		error_message("katss: sample: %s\n", seqfstrerror_r(seqferrno, buffer, BUFFER_SIZE));
		katss_free_counter(counter);
		counter = NULL;
	}

	free(buffer);
cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	return counter;
}


static int
count_file_bootstrap_mt(void *arg)
{
	threadinfo *args = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE * sizeof *buffer);
	thread_safe_rand_t *tsr = thread_safe_rand_init();

	KatssHasher *hasher = katss_init_hasher(args->kmer, 's');
	if(hasher == NULL)
		return 1;

	/* Megabyte to store counts */
	size_t num_counts = 250000;
	uint32_t *hash_values = s_malloc(num_counts * sizeof *hash_values);
	size_t cur_hash = 0;

	/* Begin counting */
	while(seqfgets(args->seqfile, buffer, BUFFER_SIZE)) {
		if(thread_safe_rand_r(tsr, args->seed) % 100000 >= args->sample)
			continue;
		katss_set_seq(hasher, buffer);
		while(katss_get_fh(hasher, &hash_values[cur_hash])) {
			if(++cur_hash == num_counts) { // begin flushing
				katss_increments(args->counter, hash_values, cur_hash);
				cur_hash = 0;
			}
		}
	}

	/* Flush values */
	katss_increments(args->counter, hash_values, cur_hash);

	/* Free resources */
	thread_safe_rand_free(tsr);
	free(hasher);
	free(buffer);
	free(hash_values);

	return 0;
}


KatssCounter *
katss_count_kmers_bootstrap_mt(const char *filename, unsigned int kmer,
                               int sample, unsigned int *seed, int threads)
{
	threads = MAX2(threads, 1);
	threads = MIN2(threads, 128);

	/* Process single-threaded computation */
	if(threads == 1)
		return katss_count_kmers_bootstrap(filename, kmer, sample, seed);

	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	/* Open SeqFile for reading */
	SeqFile file = seqfopen(filename, "r");
	if(file == NULL) {
		error_message("%s", seqfstrerror(seqferrno));
		return NULL;
	}

	/* Initialize counter */
	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		seqfclose(file);
		return NULL;
	}

	threadinfo *jobarg = s_malloc(threads * sizeof *jobarg);
	thrd_t *jobs = s_malloc(threads * sizeof *jobs);
	for(int i=0; i<threads; i++) {
		jobarg[i].seqfile = file;
		jobarg[i].counter = counter;
		jobarg[i].kmer = kmer;
		jobarg[i].filetype = seqftype(file);
		jobarg[i].sample = sample;
		jobarg[i].seed = seed;

		/* Start threads */
		thrd_create(&jobs[i], count_file_bootstrap_mt, &jobarg[i]);
	}

	for(int i=0; i<threads; i++) {
		thrd_join(jobs[i], NULL);
	}

	/* Free resources */
	seqfclose(file);
	free(jobs);
	free(jobarg);

	return counter;
}


static int
count_file_mt(void *arg)
{
	threadinfo *args = (threadinfo *)arg;
	char *buffer = s_malloc(BUFFER_SIZE * sizeof *buffer);

	KatssHasher *hasher = katss_init_hasher(args->kmer, args->filetype);
	if(hasher == NULL) {
		return 1;
	}

	/* Megabyte to store counts */
	size_t num_counts = 250000;
	uint32_t *hash_values = s_malloc(num_counts * sizeof *hash_values);
	size_t cur_hash = 0;

	/* Begin counting */
	while(seqfread(args->seqfile, buffer, BUFFER_SIZE)) {
		katss_set_seq(hasher, buffer);
		while(katss_get_fh(hasher, &hash_values[cur_hash])) {
			if(++cur_hash == num_counts) { // begin flushing
				katss_increments(args->counter, hash_values, cur_hash);
				cur_hash = 0;
			}
		}
	}

	/* Flush values */
	katss_increments(args->counter, hash_values, cur_hash);

	/* Free resources */
	free(hasher);
	free(buffer);
	free(hash_values);

	return 0;
}

/*==============================================================================
 Ushuffle counting functions
==============================================================================*/
KatssCounter *
katss_count_kmers_ushuffle(const char *filename, unsigned int kmer, int klet)
{
	if(klet < 1)
		return NULL;

	KatssCounter *counter = NULL;

	/* Open file and prepare counter & hasher */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	char *shuf   = s_calloc(BUFFER_SIZE, sizeof *shuf);
	uint32_t hash_value;

	/* Open file */
	SeqFile read_file = seqfopen(filename, "r");
	if(read_file == NULL)
		goto exit;

	KatssHasher *hasher = katss_init_hasher(kmer, 's');
	if(hasher == NULL)
		goto cleanup_file;

	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	srand(1); // reset rand seed for shuffle
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		int seqlen = strlen(buffer);
		shuffle(buffer, shuf, seqlen, klet);
		shuf[seqlen] = '\0';
		katss_set_seq(hasher, shuf);
		while(katss_get_fh(hasher, &hash_value)) {
			katss_increment(counter, hash_value);
		}
	}

	/* If error was encountered while reading report and return NULL */
	if(seqferrno) {
		katss_free_counter(counter);
		counter = NULL;
		error_message("katss: %d: %s", seqferrno, seqfstrerror(seqferrno));
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	free(buffer);
	free(shuf);
	return counter;
}

KatssCounter *
katss_count_kmers_ushuffle_bootstrap(const char *filename, unsigned int kmer,
                                     int klet, int sample, unsigned int *seed)
{
	/* sample should be between 1-100000 */
	sample = MAX2(sample, 1);
	sample = MIN2(sample, 100000);

	/* If not subsampling, just do regular ushuffle */
	if(sample == 100000)
		return katss_count_kmers_ushuffle(filename, kmer, klet);

	/* Check klet */
	if(klet < 1)
		return NULL;
	KatssCounter *counter = NULL;

	/* Initialize buffer */
	char *buffer = s_calloc(BUFFER_SIZE, sizeof *buffer);
	char *shuf   = s_calloc(BUFFER_SIZE, sizeof *shuf);
	uint32_t hash_value;

	/* Open file and hasher */
	SeqFile read_file = seqfopen(filename, "r");
	if(read_file == NULL) {
		error_message("%s", seqfstrerror(seqferrno));
		goto exit;
	}
	
	KatssHasher *hasher = katss_init_hasher(kmer, 's');
	if(hasher == NULL)
		goto cleanup_file;
	
	counter = katss_init_counter(kmer);
	if(counter == NULL)
		goto cleanup_hasher;

	/* int to subsample from rand() */
	unsigned int local_seed;
	if(seed == NULL) {
		local_seed = time(NULL);
		seed = &local_seed;
	}

	srand(1); // reset rand seed for shuffle
	while(seqfgets_unlocked(read_file, buffer, BUFFER_SIZE)) {
		/* Pick random sequences */
		if(rand_r(seed) % 100000 >= sample)
			continue;
		/* Shuffle sequences */
		int seqlen = strlen(buffer);
		shuffle(buffer, shuf, strlen(buffer), klet);
		shuf[seqlen] = '\0'; // add null terminator since shuffle uses strncpy
		katss_set_seq(hasher, shuf);
		while(katss_get_fh(hasher, &hash_value)) {
			katss_increment(counter, hash_value);
		}
	}

	if(seqferrno) {
		error_message("katss: sample: %s\n", seqfstrerror_r(seqferrno, buffer, BUFFER_SIZE));
		katss_free_counter(counter);
		counter = NULL;
	}

cleanup_hasher:
	free(hasher);
cleanup_file:
	seqfclose(read_file);
exit:
	free(buffer);
	free(shuf);
	return counter;
}

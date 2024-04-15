#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

typedef struct KmerCounter {
	unsigned int k_mer;         /**< Length of k-mers to be counted. */
	unsigned int capacity;      /**< Capacity of the k-mer counter (4^k_mer). */
	unsigned int *entries;      /**< Array to store counts of each k-mer. */
	unsigned long total_count;
	unsigned int nucleotide_to_number[256];
} KmerCounter;


/**
 * @brief Initializes a k-mer counter.
 *
 * Allocates memory for a new KmerCounter structure and initializes its fields.
 *
 * @param k_mer Length of k-mers to be counted.
 * @return A pointer to the initialized KmerCounter.
 */
KmerCounter *init_kcounter(unsigned int k_mer);


/**
 * @brief Frees the memory allocated for a KmerCounter structure.
 *
 * Deallocates all memory associated with the provided KmerCounter structure,
 * including the entries array and the structure itself.
 *
 * @param kmer_counter Pointer to the KmerCounter structure to be freed.
 */
void free_kcounter(KmerCounter *kmer_counter);


/**
 * @brief Increments k-mer counts for all substrings in a given sequence.
 *
 * Processes the input sequence and increments the counts of all k-mers in the k-mer counter.
 *
 * @param kcounter Pointer to the KmerCounter structure.
 * @param sequence Input sequence to process.
 */
void kctr_increment(KmerCounter *kcounter, const char *sequence);


/**
 * @brief Increments k-mer counts for all substrings in a given sequence.
 * 
 * Nearly the same as kctr_increment but faster by using a rolling hash. The drawback is if it
 * encounters a non-nucleotide character, the rolling hash will not work as intended. As such,
 * you should only use this on sequences you know will only contain the characters ATCG.
 * 
 * @note Use for sequences containing only `ACTGU` characters.
 * 
 * @param kcounter Pointer to the KmerCounter structure.
 * @param sequence Input sequence to process.
*/
void kctr_fincrement(KmerCounter *kcounter, const char *sequence);


/**
 * @brief Decrements the the counts of all k-mers affected by the k-mer pat
 * 
 * @param kcounter Pointer to the KmerCounter structure.
 * @param sequence Input sequence to process.
 * @param pat k-mer pattern to decrement from.
*/
void kctr_decrement(KmerCounter *kcounter, const char *sequence, const char *pat);


/**
 * @brief Increments the count for a specific substring in a given sequence.
 *
 * Increments the count of the k-mer formed by the substring starting at 'start' with the specified length.
 *
 * @param kmer_counter Pointer to the KmerCounter structure.
 * @param sequence Input sequence to process.
 * @param start Starting position of the substring.
 * @param length Length of the substring.
 */
void kctr_increment_substr(KmerCounter *kmer_counter, char *sequence,
                           unsigned int start, unsigned int length);

/**
 * @brief Retrieves the count of a specific k-mer.
 *
 * Retrieves the count of the specified k-mer from the k-mer counter.
 *
 * @param kmer_counter Pointer to the KmerCounter structure.
 * @param key The k-mer sequence to query.
 * @return The count of the specified k-mer.
 */
unsigned int kctr_get(KmerCounter *kmer_counter, char *key);


/**
 * @brief Empties the specific key.
 * 
 * If any counts were associated with the given key, it resets it and sets
 * the counts to 0.
 * 
 * @param kmer_counter Pointer to the KmerCounter structure.
 * @param key The k-mer key to empty
*/
void kctr_empty(KmerCounter *kmer_counter, char *key);


/**
 * @brief Retrieves the k-mer associated with a specific index.
 *
 * Retrieves the k-mer sequence associated with the specified index in the k-mer counter.
 *
 * @param kmer_counter Pointer to the KmerCounter structure.
 * @param index The index of the k-mer.
 * @return The k-mer sequence associated with the specified index.
 */
char *kctr_get_key(KmerCounter *kmer_counter, unsigned int index);

#endif // KMER_COUNTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

#include "rna_file_parser.h"
#include "kmer_counter.h"
#include "kmerHashTable.h"
#include "string_utils.h"
#include "utils.h"

#include "SKA_cmdl.h"

typedef struct options {
	char    *input_file;
	char    *bound_file;
	char    *out_filename;
	FILE    *out_file;
	int     out_given;
	int     kmer;
	int     iterations;
	char    file_delimiter;

	int     independent_probs;
	char	**motif;

	char    **top_kmer;
	int     cur_iter;
} options;


typedef struct {
	kmerHashTable *monomer_frq;
	kmerHashTable *dimer_frq;
	kmerHashTable *kmer_frq;
} FrqIndependentProbs;


/*##########################################################
#  Function Declarations                                   #
##########################################################*/

void 
process_iteration(options *opt);


KmerCounter *
count_kmers(char *filename, options *opt);


FrqIndependentProbs
process_independent_probs(char *filename, options *opt);


kmerHashTable *
predict_kmers(kmerHashTable *probs_1mer, kmerHashTable *probs_2mer, int kmer);


kmerHashTable *
get_frequencies(KmerCounter *counts);


kmerHashTable *
get_enrichment(kmerHashTable *input_frq, kmerHashTable *bound_frq, options *opt);


Entry *
kmer_max_entry(kmerHashTable *hash_table);


char
delimiter_to_char(char *user_delimiter);


void
entry_to_file(FILE *file, Entry *entry, char delimiter);


void
free_options(options *opt);


void
print_options(options *opt);


void 
init_default_options(options *opt)
{
	opt->input_file     = NULL;
	opt->bound_file     = NULL;
	opt->out_filename   = "motif";
	opt->out_file       = NULL;
	opt->out_given      = 0;
	opt->kmer           = 3;
	opt->iterations     = 1;
	opt->file_delimiter = ',';

	opt->independent_probs  = 0;
	opt->motif              = NULL;

	opt->top_kmer       = NULL;
	opt->cur_iter       = 0;
}


/*##########################################################
#  Main                                                    #
##########################################################*/
int
main(int argc, char **argv)
{
	struct SKA_args_info    args_info;
	options                 opt;

	init_default_options(&opt);

	/*##########################################################
	#  Parse Command Line Arguments                            #
	##########################################################*/

	if (SKA_cmdline_parser(argc, argv, &args_info) != 0) {
		exit(EXIT_FAILURE);
	}

	if(args_info.input_given) {
		opt.input_file = strdup(args_info.input_arg);
		if(access(opt.input_file, F_OK|R_OK) != 0) {
			error_message("Unable to open input file '%s' for reading.",args_info.input_arg);
			SKA_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
	}

	if(args_info.bound_given) {
		opt.bound_file = strdup(args_info.bound_arg);
		if(access(opt.bound_file, F_OK|R_OK) != 0) {
			error_message("Unable to open bound file '%s' for reading.",args_info.bound_arg);
			SKA_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
	}

	if(args_info.output_given) {
		opt.out_given    = 1;
		opt.out_filename = strdup(args_info.output_arg);
	}

	if(args_info.kmer_given) {
		if(args_info.kmer_arg <= 0) {
			error_message("option '--kmer=%d' must be a value greater than 0.",args_info.kmer_arg);
			SKA_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.kmer = args_info.kmer_arg;
	}

	if(args_info.iterations_given) {
		if(args_info.iterations_arg <= 0) {
			error_message("option '--iterations=%d' must be a value greater than 0.",
			 args_info.kmer_arg);
			SKA_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.iterations = args_info.iterations_arg;

		if(opt.iterations > 1<<(2*opt.kmer)) { /* iterations can't be greater than k-mers */
			opt.iterations = 1<<(2*opt.kmer);
		}
	}

	if(args_info.file_delimiter_given) {
		opt.file_delimiter = delimiter_to_char(args_info.file_delimiter_arg);
	}

	if(args_info.independent_probs_given) {
		opt.independent_probs = 1;
	}

	if(args_info.motif_given) {
		opt.motif = args_info.motif_arg;
	}

	char *filename = concat(opt.out_filename, ".dsv");
	opt.out_file = fopen(filename, "w");
	if (opt.out_file == NULL) {
		error_message("Could not write to file '%s'\n",filename);
	}

	opt.top_kmer = s_malloc(opt.iterations * sizeof(char *));
	for(int i=0; i<opt.iterations; i++) {
		opt.top_kmer[i] = NULL;
	}

	SKA_cmdline_parser_free(&args_info);
	free(filename);

	/*##########################################################
	#  Computations                                            #
	##########################################################*/

	while(opt.cur_iter < opt.iterations) {
		process_iteration(&opt);
	}

	/* Clean up */
	free_options(&opt);

	return 0;
}


/*##########################################################
#  Main Functions                                          #
##########################################################*/

void 
process_iteration(options *opt)
{
	kmerHashTable *input_table;
	kmerHashTable *bound_table;
	kmerHashTable *enrichments_table;

	if(opt->independent_probs) {
		FrqIndependentProbs kmer_data;
		kmer_data = process_independent_probs(opt->bound_file, opt);

		input_table = predict_kmers(kmer_data.monomer_frq, kmer_data.dimer_frq, opt->kmer);
		bound_table = kmer_data.kmer_frq;

		free_kmer_table(kmer_data.monomer_frq);
		free_kmer_table(kmer_data.dimer_frq);
	} else if(opt->motif) {

	} else {
		KmerCounter *input_counts = count_kmers(opt->input_file, opt);
		input_table = get_frequencies(input_counts);
		free_kcounter(input_counts);

		KmerCounter *bound_counts = count_kmers(opt->bound_file, opt);
		bound_table = get_frequencies(bound_counts);
		free_kcounter(bound_counts);
	}

	enrichments_table = get_enrichment(input_table, bound_table, opt);

	Entry *max_entry = kmer_max_entry(enrichments_table);
	entry_to_file(opt->out_file, max_entry, opt->file_delimiter);

	opt->top_kmer[opt->cur_iter] = strdup(max_entry->key);
	opt->cur_iter++;

	/* free tables */
	free_kmer_table(enrichments_table);
	free_kmer_table(input_table);
	free_kmer_table(bound_table);
}


KmerCounter *
count_kmers(char *filename, options *opt)
{
	RNA_FILE *read_file = rnaf_open(filename);
	KmerCounter *counter = init_kcounter(opt->kmer);

	char *sequence;
	while(1) {
		sequence = rnaf_get(read_file);

		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 1);
		for(int i=0; i<opt->cur_iter; i++) {
			cross_out(sequence, opt->top_kmer[i]);
		}

		kctr_increment(counter, sequence);
		free(sequence);
	}
	rnaf_close(read_file);

	return counter;
}


FrqIndependentProbs
process_independent_probs(char *filename, options *opt)
{
	FrqIndependentProbs     kmer_data;
	RNA_FILE *read_file = rnaf_open(filename);

    KmerCounter *monomers_cnt = init_kcounter(1);
    KmerCounter *dimers_cnt = init_kcounter(2);
    KmerCounter *counts_cnt = init_kcounter(opt->kmer);

	/* Count di/monomers & k-mers */
	char *sequence = NULL;
	while(1) {
		sequence = rnaf_get(read_file);
		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 1);
		for(int i=0; i<opt->cur_iter; i++) {
			cross_out(sequence, opt->top_kmer[i]);
		}

		kctr_increment(monomers_cnt, sequence);
		kctr_increment(dimers_cnt, sequence);
		kctr_increment(counts_cnt, sequence);

		free(sequence);
	}
	rnaf_close(read_file);

	/* Get frequencies of the counts */
	kmer_data.monomer_frq = get_frequencies(monomers_cnt);
	kmer_data.dimer_frq = get_frequencies(dimers_cnt);
	kmer_data.kmer_frq = get_frequencies(counts_cnt);

	/* Clean up data */
	free_kcounter(monomers_cnt);
	free_kcounter(dimers_cnt);
	free_kcounter(counts_cnt);

    return kmer_data;
}


kmerHashTable *
predict_kmers(kmerHashTable *probs_1mer, kmerHashTable *probs_2mer, int kmer)
{
    kmerHashTable   *predicted_kmers;
    double          dinucleotides_prob;
    double          monomers_probs;
    double          *prob;
	char			*k_str;
    char            *k_substr;

    predicted_kmers = init_kmer_table(kmer,1);
	
	k_str = s_malloc((kmer + 1) * sizeof(char));
	k_str[kmer] = '\0';

    /* Get the probability of kmer being in reads */
    for(unsigned long i = 0; i < predicted_kmers->capacity; i++) {
		int index = i;	/* Create next k-mer */
		for(int j = kmer-1; j >= 0; j--) {
			k_str[j] = BASES[index % 4];
			index /= 4;
		}

		/* Calculate probability of created k-mer string */
        dinucleotides_prob = 1;
        monomers_probs     = 1;
        for(int j=0; j < kmer-1; j++) {
            k_substr = substr(k_str, j, 2);
            prob = kmer_get(probs_2mer, k_substr);
            free(k_substr);
            dinucleotides_prob*=prob[0];
        }
        for(int j=1; j < kmer-1; j++) {
            k_substr = substr(k_str, j, 1);
            prob = kmer_get(probs_1mer, k_substr);
            free(k_substr);
            monomers_probs*=prob[0];
        }
        kmer_add_value(predicted_kmers, k_str, dinucleotides_prob/monomers_probs, 0);
    }
	free(k_str);

    return predicted_kmers;
}


void
count_motifs(char *filename, options *opt)
{
	/*
		Idea:
		**motifs is an array of char arrays, as such to access the strings you need
		to index 'motifs'. The kcounter structure will use these indeces as the mapping
		to the counter, which will allow for fast lookup in the addition for multiple patterns.

		Alongside this, we can use the boyer-moore algorithm to find the number of matches. We can
		also build a regexp interpreter to handle regular expressions.

		To exploit the boyer-moore, we need to modify RNA_FILE struct to include a large buffer
		that we can search through. We can do this by having an rnaf_open function that takes
		a buffer size parameter.
	*/
}


kmerHashTable *
get_frequencies(KmerCounter *counts)
{
	kmerHashTable *frq_table = init_kmer_table(counts->k_mer, 1);
	unsigned long total_kcounts = counts->total_count;

	unsigned int kcount;
	for(unsigned int i=0; i<counts->capacity; i++) {
		if(!(kcount=counts->entries[i])) {
			continue;
		}

		char *key = kctr_get_key(counts, i);
		double k_frequency = (double)kcount/total_kcounts; // calculate frequency
		kmer_add_value(frq_table, key, k_frequency, 0);
		free(key);
	}

	return frq_table;
}


kmerHashTable *
get_enrichment(kmerHashTable *input_frq, kmerHashTable *bound_frq, options *opt)
{
    char            *key;
    double          *input_values;
    double          *bound_values;
    double          enrichment;

    kmerHashTable *enrichments_table = init_kmer_table(opt->kmer, 1);
    
    // Get the log2 fold change for each kmer
    for(unsigned long i = 0; i < enrichments_table->capacity; i++) {
        if(bound_frq->entries[i] == NULL || input_frq->entries[i] == NULL) {
            continue;
        }

        key = bound_frq->entries[i]->key;

        bound_values = kmer_get(bound_frq, key);
        input_values = kmer_get(input_frq, key);

        if(bound_values[0] == 0.0f || input_values[0] == 0.0f) {
            enrichment = 0.0f;
        } else {
            enrichment = logf(bound_values[0]/input_values[0])/logf(2.0f);
        }
        kmer_add_value(enrichments_table, key, enrichment, 0);
    }

    return enrichments_table;
}

/*##########################################################
#  Helper Functions                                        #
##########################################################*/

Entry *
kmer_max_entry(kmerHashTable *hash_table)
{
    double max_value = -DBL_MAX;
    Entry *max_entry = NULL;
    for(unsigned long i=0; i<hash_table->capacity; i++) {
        if(hash_table->entries[i] == NULL) {
            continue;
        }
        
        if(hash_table->entries[i]->values[0] > max_value) {
            max_value = hash_table->entries[i]->values[0];
            max_entry = hash_table->entries[i];
        }
    }

    return max_entry;
}


char 
delimiter_to_char(char *user_delimiter)
{
    char delimiter;
    char provided_val = user_delimiter[0];
    if(strlen(user_delimiter)>1) {
        provided_val = 0;
    }

    switch(provided_val) {
        case ',':
            delimiter = ',';
            break;
        case 't':
            delimiter = '\t';
            break;
        case ':':
            delimiter = ':';
            break;
        case '|':
            delimiter = '|';
            break;
        case ' ':
            delimiter = ' ';
            break;
        default:
            warning_message("Provided delimiter '%s' is not valid. Defaulting to ','",
            user_delimiter);
            delimiter = ',';
    }

    return delimiter;
}


void 
entry_to_file(FILE *file, Entry *entry, char delimiter)
{
    fprintf(file, "%s%c%f\n", entry->key, delimiter, entry->values[0]);
}


void 
free_options(options *opt)
{
    if(opt->input_file && !opt->independent_probs) {
        free(opt->input_file);
    }
    if(opt->bound_file) {
        free(opt->bound_file);
    }
    if(opt->out_given) {
        free(opt->out_filename);
    }
    for(int i=0; i<opt->iterations; i++) {
        free(opt->top_kmer[i]);
    }
    free(opt->top_kmer);
    fclose(opt->out_file);
}


void
print_options(options *opt)
{
    printf("input_file: '%s'\n",opt->input_file);
    printf("bound_file: '%s'\n",opt->bound_file);
    printf("output_file: '%s'\n",opt->out_filename);
    printf("kmer: '%d'\n",opt->kmer);
    printf("iterations: '%d'\n",opt->iterations);
    printf("probs: '%d'\n",opt->independent_probs);
}

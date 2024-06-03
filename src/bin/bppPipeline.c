#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "rna_file_parser.h"
#include "kmerHashTable.h"
#include "string_utils.h"
#include "utils.h"
#include "parallel_helpers.h"

#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/utils/basic.h"

#include "bppPipeline_cmdl.h"

typedef struct {
	char    *input_file;
	char    *bound_file;
	char    *out_file;
	int     output_given;
	int     kmer;
	char    file_delimiter;
	int     jobs;

	int     keepFolds;
	FILE    *folds_file;
	int     input_fold;
	int     bound_fold;

	int     seq_windows;
	int     window_size;
	int     bin;
	int     frq;
} options;


typedef struct {
	char *sequence;
	kmerHashTable *counts_table;
	options *opt;
} record_data;


/*##########################################################
#  Function Declarations                                   #
##########################################################*/

float *
getPositionalProbabilities(char *sequence);


void
print_options(options *opt);


void
free_options(options *opt);


int
compare(const void  *a, 
        const void  *b);


void
print_table_to_file(kmerHashTable  *table,
                    FILE           *table_file,
                    char           sep);


kmerHashTable *
bpp_kmer_frequency(char *filename,
                   options *opt,
                   int folds_provided);


void
process_record(record_data *record, int folds_provided);


void
process_seq(record_data *record);


int
process_line_with_bpp(char          *line, 
                      kmerHashTable *counts_table, 
                      int           kmer);


void
process_windows(record_data *record);


void
getFrequencies(kmerHashTable *counts_table);


kmerHashTable *
getBPPEnrichment(kmerHashTable *control_frq, 
                 kmerHashTable *bound_frq, 
                 int kmer);


record_data *
copy_record_data(record_data *data);


void
stream_folds(char *sequence, float *positional_probabilities, options *opt);


char
delimiter_to_char(char *user_delimiter);


void
init_default_options(options *opt)
{
	opt->input_file     = NULL;
	opt->bound_file     = NULL;
	opt->out_file       = "rna";
	opt->output_given   = 0;
	opt->kmer           = 3;
	opt->file_delimiter = ',';
	opt->jobs           = 0;

	opt->keepFolds      = 0;
	opt->input_fold     = 0;
	opt->bound_fold     = 0;
	opt->folds_file     = NULL;

	opt->seq_windows    = 0;
	opt->window_size    = 20;
	opt->bin            = 0;
	opt->frq            = 0;
}


/*##########################################################
#  Main                                                    #
##########################################################*/
int
main(int argc, char **argv)
{

	// Declare variables
	struct bppPipeline_args_info    args_info;
	options                         opt;
	kmerHashTable                   *bounds_table;
	kmerHashTable                   *control_table;
	kmerHashTable                   *enrichments_table;

	init_default_options(&opt);

	/*##########################################################
	#  Parse Command Line Arguments                            #
	##########################################################*/

	if (bppPipeline_cmdline_parser(argc, argv, &args_info) != 0) {
		exit(EXIT_FAILURE);
	}

	if(args_info.input_given) {
		opt.input_file = strdup(args_info.input_arg);
		if(access(opt.input_file, F_OK|R_OK) != 0) {
			error_message("Unable to open input file '%s' for reading.",args_info.input_arg);
			bppPipeline_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
	}

	if(args_info.bound_given) {
		opt.bound_file = strdup(args_info.bound_arg);
		if(access(opt.bound_file, F_OK|R_OK) != 0) {
			error_message("Unable to open bound file '%s' for reading.",args_info.bound_arg);
			bppPipeline_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
	}

	if(args_info.output_given) {
		opt.out_file = strdup(args_info.output_arg);
		opt.output_given = 1;
	}

	if(args_info.kmer_given) {
		if(args_info.kmer_arg <= 0) {
			error_message("option '--kmer=%d' must be a value greater than 0.",args_info.kmer_arg);
			bppPipeline_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.kmer = args_info.kmer_arg;
	}

	if(args_info.file_delimiter_given) {
		opt.file_delimiter = delimiter_to_char(args_info.file_delimiter_arg);
	}

	if(args_info.keepFolds_given) {
		opt.keepFolds = 1;
	}

	if(args_info.input_fold_given) {
		opt.input_fold = 1;
		opt.input_file = strdup(args_info.input_fold_arg);
	}

	if(args_info.bound_fold_given) {
		opt.bound_fold = 1;
		opt.bound_file = strdup(args_info.bound_fold_arg);
	}

	if(args_info.seq_windows_given) {
		if(args_info.seq_windows_arg <= 0) {
			error_message("option '--seq-windows=%d' must be greater than 0.",
			args_info.seq_windows_arg);
			bppPipeline_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.seq_windows = 1;
		opt.window_size = args_info.seq_windows_arg;
	}

	if(args_info.jobs_given) {
		int thread_max = max_user_threads();

		if (args_info.jobs_arg == 0) {
		/* use maximum of concurrent threads */
			int proc_cores, proc_cores_conf;
			if (num_proc_cores(&proc_cores, &proc_cores_conf)) {
				opt.jobs = MIN2(thread_max, proc_cores_conf);
			} else {
				warning_message("Could not determine number of available processor cores!\n"
				                "Defaulting to serial computation");
				opt.jobs = 1;
			}
		} else {
			opt.jobs = MIN2(thread_max, args_info.jobs_arg);
		}
		opt.jobs = MAX2(1, opt.jobs);
	}

	if(args_info.frq_given) {
		opt.frq = 1;
	}

	if(args_info.bin_given) {
		opt.bin = 1;
	}

	bppPipeline_cmdline_parser_free(&args_info);

	/*##########################################################
	#  Computations                                            #
	##########################################################*/

	INIT_PARALLELIZATION(opt.jobs);
	control_table = bpp_kmer_frequency(opt.input_file, &opt, opt.input_fold);
	bounds_table  = bpp_kmer_frequency(opt.bound_file, &opt, opt.bound_fold);
	UNINIT_PARALLELIZATION

	enrichments_table = getBPPEnrichment(control_table, bounds_table, opt.kmer);

	kmerHashTable_to_file(enrichments_table, opt.out_file, opt.file_delimiter);

	/* Clean up */
	free_kmer_table(control_table);
	free_kmer_table(bounds_table);
	free_kmer_table(enrichments_table);
	free_options(&opt);

	return 0;
}


/*##########################################################
#  Main Functions                                          #
##########################################################*/

float *
getPositionalProbabilities(char *sequence)
{
	vrna_ep_t   *ptr, *pair_probabilities = NULL;
	int         seq_length                = strlen(sequence);
	char        *propensity               = (char *)vrna_alloc(sizeof(char) * (seq_length + 1));
	float       *positional_probabilities = s_calloc(seq_length, sizeof(float));
	float       probability;

	/* Get the pair probabilities */
	vrna_pf_fold(sequence, propensity, &pair_probabilities);

	/* Move pair probabilities into array */
	for(ptr = pair_probabilities; ptr->i != 0; ptr++) {
		probability = ptr->p;
		positional_probabilities[ptr->i-1]+=probability;
		positional_probabilities[ptr->j-1]+=probability;
	}

	/* Clean up memory */
	free(propensity);
	free(pair_probabilities);

	return positional_probabilities;
}


kmerHashTable *
bpp_kmer_frequency(char *filename, options *opt, int folds_provided)
{
	record_data     *record;
	kmerHashTable   *counts_table;
	RNA_FILE        *read_file;
	char            *sequence;

	counts_table = init_bpp_table(opt->kmer);
	read_file = rnaf_open(filename);

	record = s_malloc(sizeof *record);
	record->counts_table = counts_table;
	record->opt = opt;

	/* Create folds file to write to */
	if(opt->keepFolds) {
		char *filename_prefix = basename_prefix(filename);
		char *folds_filename = concat(filename_prefix, ".folds");
		opt->folds_file = fopen(folds_filename, "w");

		free(filename_prefix);
		free(folds_filename);
	}

	/* Main loop to process each line in file */
	while (1) {
		sequence = rnaf_get(read_file);

		if(sequence == NULL) {
			break;  /* Reached EOF */
		}

		/* Pre-processing on sequence */
		seq_to_RNA(sequence);
		str_to_upper(sequence);
		remove_escapes(sequence);
		
		record->sequence = sequence;

		/* Function that will determine how to process the record */
		process_record(record, folds_provided);

		free(record->sequence);
	}

	/* Clean up*/
	rnaf_close(read_file);
	free(record);
	WAIT_FOR_THPOOL

	if(opt->keepFolds) {
		fclose(opt->folds_file);
	}

	/* Calculate the %frequencies */
	getFrequencies(counts_table);

	return counts_table;
}


void
process_record(record_data *record, int folds_provided)
{
	/* If folds file from --keepFolds was provided, then process line using bpp values*/
	if(folds_provided) {
		process_line_with_bpp(record->sequence, record->counts_table, record->opt->kmer);

	/* If --seq-windows was provided, use the sliding window algorithm for calculations */
	} else if(record->opt->seq_windows) {
		WAIT_FOR_FREE_SLOT((2 * record->opt->jobs) - 1 );
		RUN_IN_PARALLEL(process_windows, copy_record_data(record));

	/* Default algorithm for getting BPP frequencies */
	} else {
		WAIT_FOR_FREE_SLOT((2 * record->opt->jobs) - 1);
		RUN_IN_PARALLEL(process_seq, copy_record_data(record));
	}
	/* ...if more algorithms to process the record are added, add them here. */
}


/*##########################################################
#  Algorithms - functions responsible for computations     #
##########################################################*/

void
process_seq(record_data *record)
{
	options *opt = record->opt;
	char    *sequence = record->sequence;

	float   *positional_probabilities;
	char    *k_substr;
	int     seq_length = strlen(sequence);
	int     num_kmers_in_seq = seq_length - opt->kmer + 1;

	positional_probabilities = getPositionalProbabilities(sequence);
	THREADSAFE_STREAM_OUTPUT(stream_folds(sequence, positional_probabilities, opt));

	/* Count kmers and their associated base-pair probability */
	for(int i=0; i<num_kmers_in_seq; i++) {
		k_substr = substr(sequence, i, opt->kmer);

		kmer_add_value(record->counts_table, k_substr, 1, opt->kmer);
		
		// Loop through bpp values in file
		for(int j=i; j<opt->kmer+i; j++) {
			kmer_add_value(record->counts_table, k_substr, positional_probabilities[j], j-i);
		}

		free(k_substr);
	}

	free(positional_probabilities);
	free(record->sequence);
	free(record);
}


int
process_line_with_bpp(char *line, kmerHashTable *counts_table, int kmer)
{
	char    *data;
	char    *sequence;
	char    *k_substr;
	double  bpp;
	int     seq_length;
	int     num_kmers_in_seq;

	data = strtok(line, " ");
	seq_length = strlen(data);
	sequence = strdup(data);

	num_kmers_in_seq = seq_length - kmer + 1;
	double token_bpp_values[num_kmers_in_seq];

	data = strtok(NULL, " ");
	int token_count = 0;
	for(int i = 0; i<num_kmers_in_seq; i++) {
		if(data == NULL) {
			error_message("File is formatted incorrectly.\nMake sure you are using the '.folds' file produced when running the pipeline with the '--keepFolds' command.");
			exit(EXIT_FAILURE);
		}

		bpp = atof(data);
		token_bpp_values[token_count++]=bpp;
		data = strtok(NULL, " ");
	}

	for(int i=0; i<num_kmers_in_seq; i++) {
		k_substr = substr(sequence, i, kmer);
		kmer_add_value(counts_table, k_substr, 1, kmer); // increase count by 1 for kmer
		
		// Loop through bpp values in file
		for(int j=i; j<kmer+i; j++) {
			kmer_add_value(counts_table, k_substr, token_bpp_values[j], j-i);
		}

		free(k_substr);
	}
	free(sequence);

	return 0; // ran without error
}


void
process_windows(record_data *record)
{
	char *sequence  = record->sequence;
	options *opt    = record->opt;
	kmerHashTable *counts_table = record->counts_table;

	char    *k_substr;
	char    *window_seq;
	float   *window_probabilities;
	float   mean_probability;
	float   sum_probability;
	int     num_windows;
	int     seq_length;
	int     count_probs;

	// Initialize probability matrix with -1
	seq_length = strlen(sequence);
	num_windows = seq_length - opt->window_size + 1;
	float probability_matrix[num_windows][seq_length];
	for(int row=0; row<num_windows; row++){
		for(int col=0; col<seq_length; col++) {
			probability_matrix[row][col]=-1;
		}
	}

	// Fill the probability matrix with probabilities
	for(int i = 0; i<num_windows; i++) {
		window_seq = substr(sequence, i, opt->window_size);
		window_probabilities = getPositionalProbabilities(window_seq);

		for(int j = 0; j < opt->window_size; j++) {
			probability_matrix[i][j+i] = window_probabilities[j];
		}

		free(window_seq);
		free(window_probabilities);
	}

	// Get the average of each column in probability matrix
	float positional_probabilities[seq_length];
	for(int col=0; col<seq_length; col++) {
		sum_probability = 0;
		count_probs     = 0;

		for(int row=0; row<num_windows; row++) {
			if(probability_matrix[row][col] == -1) {
				continue;
			}
			sum_probability+=probability_matrix[row][col];
			count_probs++;
		}

		mean_probability = sum_probability/count_probs;
		positional_probabilities[col] = mean_probability;
	}

	/* Output positional probabilities to file if necessary */
	THREADSAFE_STREAM_OUTPUT(stream_folds(sequence, positional_probabilities, opt));

	/* Fill counts_table with positional probabilities */
	int num_kmers_in_seq = seq_length - opt->kmer + 1;
	for(int i=0; i<num_kmers_in_seq; i++) {
		k_substr = substr(sequence, i, opt->kmer);

		kmer_add_value(counts_table, k_substr, 1, opt->kmer);
		
		// Loop through bpp values in file
		for(int j=i; j<opt->kmer+i; j++) {
			kmer_add_value(counts_table, k_substr, positional_probabilities[j], j-i);
		}

		free(k_substr);
	}
}


void 
getFrequencies(kmerHashTable *counts_table)
{
	int num_columns = counts_table->cols-1;
	int total_count;

	for(int i=0; i<counts_table->capacity; i++) {
		if( !counts_table->entries[i] || counts_table->entries[i]->values[num_columns]<1) {
			continue;
		}
		for(int j=0; j<num_columns; j++) {
			total_count=(int)counts_table->entries[i]->values[num_columns];
			counts_table->entries[i]->values[j]/=total_count;
		}
	}
}


kmerHashTable *
getBPPEnrichment(kmerHashTable *control_frq, kmerHashTable *bound_frq, int kmer)
{
	kmerHashTable   *enrichments_table;
	char            *key;
	double          enrichment;
	double          *bound_values;
	double          *control_values;
	double          *enrichment_values;
	double          mean_enrichment;

	enrichments_table = init_bpp_table(kmer);

	// Get the log2 fold change for each kmer
	for(size_t i = 0; i < bound_frq->capacity; i++) {
		if(bound_frq->entries[i] == NULL || control_frq->entries[i] == NULL) {
			continue;
		}

		key = bound_frq->entries[i]->key;

		bound_values   = kmer_get(bound_frq, key);
		control_values = kmer_get(control_frq, key);

		for(int j = 0; j < kmer; j++) {
			if(bound_values[j] == 0. || control_values[j] == 0.) {
				enrichment = 0.;
			} else {
				enrichment = log2(bound_values[j]/control_values[j]);
			}
			kmer_add_value(enrichments_table, key, enrichment, j);
		}
	}

	for(size_t i = 0; i < enrichments_table->capacity; i++) {
		if(enrichments_table->entries[i] == NULL) {
			continue;
		}

		key = enrichments_table->entries[i]->key;
		enrichment_values = kmer_get(enrichments_table, key);

		mean_enrichment = 0.0f;
		for(int j = 0; j < kmer; j++) {
			mean_enrichment += enrichment_values[j];
		}

		mean_enrichment/=kmer;
		kmer_add_value(enrichments_table, key, mean_enrichment, kmer);
	}

	qsort(enrichments_table->entries, enrichments_table->capacity,
			sizeof(*enrichments_table->entries), compare);

	return enrichments_table;
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

record_data *
copy_record_data(record_data *data)
{
	record_data *new_data = s_malloc(sizeof *new_data);

	new_data->sequence      = strdup(data->sequence);
	new_data->opt           = data->opt;
	new_data->counts_table  = data->counts_table;

	return new_data;
}


void
stream_folds(char *sequence, float *positional_probabilities, options *opt)
{
	if(!opt->keepFolds) { // check if folds should be streamed
		return;
	}

	int seq_length = strlen(sequence);

	fprintf(opt->folds_file, "%s", sequence);
	for(int i=0; i<seq_length; i++) {
		fprintf(opt->folds_file, " %8.6f", positional_probabilities[i]);
	}
	fprintf(opt->folds_file, "\n");
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


int
compare(const void *a, const void *b)
{
	const Entry *entryA = *(Entry **)a;
	const Entry *entryB = *(Entry **)b;

	if(!entryA && !entryB) {
		return 0;
	}
	if(!entryA) {
		return 1;
	}
	if(!entryB) {
		return -1;
	}

	int mean_index = strlen(entryA->key);
	if(entryA->values[mean_index] > entryB->values[mean_index]) {
		return -1;
	}
	if(entryA->values[mean_index] < entryB->values[mean_index]) {
		return 1;
	}

	return 0;
}


void
free_options(options *opt)
{
	if(opt->input_file) {
		free(opt->input_file);
	}
	if(opt->bound_file) {
		free(opt->bound_file);
	}
	if(opt->output_given) {
		free(opt->out_file);
	}
}


void
print_options(options *opt)
{
	printf("Input file: \"%s\"\n",opt->input_file);
	printf("Bound file: \"%s\"\n",opt->bound_file);
	printf("Output file: \"%s\"\n",opt->out_file);
	printf("Kmer: \"%d\"\n",opt->kmer);
	printf("Seq Windows: \"%d\"\n",opt->seq_windows);
	printf("Bins: \"%d\"\n",opt->bin);
	printf("Include frq: \"%d\"\n",opt->frq);
}

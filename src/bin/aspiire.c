#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "rna_file_parser.h"
#include "ire_fold.h"
#include "Regex.h"
#include "levenshtein.h"
#include "string_utils.h"
#include "utils.h"
#include "parallel_helpers.h"

#include "ViennaRNA/fold.h"

#include "aspiire_cmdl.h"

#define BUFFER_SIZE 16777216


typedef struct SearchInfo {
	char *sequence;
	size_t sequence_length;
	size_t search_index;
	int cur_motif;
} SearchInfo;


typedef struct Options {
	char    *input_file;
	FILE    *out_file;
	char    *out_filename;
	bool     out_given;

	int      kmer;
	char     out_delimiter;
	int      threshold;

	char    *format;
	bool     format_given;

	bool     no_rnafold;

	int      jobs;
} Options;


typedef struct RecordData {
	char    *sequence;
	char    *header;
	Regex   *regex;
	Options *opt;
} RecordData;

/*==================== Function Declarations ====================*/
void    process_sequence(RecordData *record);
static Matcher search_for_ire(SearchInfo *search, Regex *regex);
static void    compare_rnafold(IreStructure *ire_structure);
static float   compare_structures(const char *struct_iref, const char *struct_vrna);

static int     get_num_threads(struct aspiire_args_info args_info);
static Regex  *init_regex(void);
RecordData    *copy_record_data(RnaInfo rna_info, Regex *regex, Options *opt);
static char    delimiter_to_char(char *user_delimiter);
static void    free_options(Options *opt);
static const char *get_contraints(const char *predicted);


void
init_deafult_options(Options *opt)
{
	opt->input_file     = NULL;
	opt->out_file       = NULL;
	opt->out_filename   = "motif";
	opt->out_given      = false;

	opt->kmer           = 30;
	opt->out_delimiter  = '\t';
	opt->threshold      = 78;

	opt->format         = "seq";
	opt->format_given   = false;

	opt->no_rnafold     = false;

	opt->jobs           = 0;
}

int
main(int argc, char *argv[])
{
	struct aspiire_args_info  args_info;
	Options                   opt;

	init_deafult_options(&opt);

	/*========================================================
	|  Parse Command Line Arguments                          |
	========================================================*/
	if (aspiire_cmdline_parser(argc, argv, &args_info) != 0) {
		exit(EXIT_FAILURE);
	}

	if(!args_info.input_given) {
		printf("Usage: SKA [OPTIONS] [<input.fa>] [<bound.fa>]\n");
		printf("Try 'aspiire --help' for more information.\n");
		aspiire_cmdline_parser_free(&args_info);
		exit(EXIT_FAILURE);
	} else {
		opt.input_file = strdup(args_info.input_arg);
		if(access(opt.input_file, F_OK|R_OK) != 0) {
			error_message("Unable to open input file '%s' for reading.",args_info.input_arg);
			aspiire_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
	}

	if(args_info.output_given) {
		opt.out_given    = true;
		opt.out_filename = strdup(args_info.output_arg);
	}

	if(args_info.kmer_given) {
		if(args_info.kmer_arg <= 0) {
			error_message("option '--kmer=%d' must be a value greater than 0.",args_info.kmer_arg);
			aspiire_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.kmer = args_info.kmer_arg;
	}

	if(args_info.file_delimiter_given) {
		opt.out_delimiter = delimiter_to_char(args_info.file_delimiter_arg);
	}

	if(args_info.threshold_given) {
		if(args_info.threshold_arg < 0) {
			error_message("option '--threshold=%d' must be a value greater than 0.",args_info.threshold_arg);
			aspiire_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		if(args_info.threshold_given > 100) {
			error_message("option '--threshold=%d' must be a value less than 100.",args_info.threshold_arg);
			aspiire_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		opt.threshold = args_info.threshold_arg;
	}

	if(args_info.format_given) {
		opt.format = strdup(args_info.format_arg);
		opt.format_given = true;
	}

	if(args_info.no_rnafold_given) {
		opt.no_rnafold = true;
	}

	opt.jobs = get_num_threads(args_info);

	char *filename = concat(opt.out_filename, ".dsv");
	opt.out_file = fopen(filename, "w");
	if (opt.out_file == NULL) {
		error_message("Could not write to file '%s'\n",filename);
		aspiire_cmdline_parser_free(&args_info);
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	aspiire_cmdline_parser_free(&args_info);
	free(filename);

	/*========================================================
	|  Main Calculations                                     |
	========================================================*/

	Regex *regex = init_regex();
	RNA_FILE *read_file = rnaf_open(opt.input_file);
	if(read_file == NULL) {
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	INIT_PARALLELIZATION(opt.jobs);
	while(true) {
		RnaInfo info = rnaf_geti(read_file);

		if(info.sequence == NULL) {
			break;
		}

		RUN_IN_PARALLEL(process_sequence, copy_record_data(info, regex, &opt));
	}
	rnaf_close(read_file);
	UNINIT_PARALLELIZATION

	/* Clean up */
	free(regex);
	free_options(&opt);
	exit(EXIT_SUCCESS);
}

/*====================  Main functions  ====================*/
void
process_sequence(RecordData *record)
{
	SearchInfo *search      = s_malloc(sizeof *search);
	search->sequence        = record->sequence;
	search->sequence_length = strlen(record->sequence);
	search->search_index    = 0;
	search->cur_motif       = 0;

	Options *opt          = record->opt;
	char    *output       = NULL;
	bool     header_shown = false;

	do {
		Matcher matcher = search_for_ire(search, record->regex);
		if(!matcher.isFound) {
			search->search_index = 0;
			search->cur_motif++; /* Match not found for specified motif, search for next one */
			continue;
		}

		/* Get the IRE motif that was found */
		int ire_start = (search->search_index + matcher.foundAtIndex) - 7;
		if(ire_start < 0) { 
			search->search_index += matcher.foundAtIndex + matcher.matchLength;
			continue;
		}
		char *ire_sequence = substr(record->sequence, ire_start, 32);
		clean_seq(ire_sequence, true);

		/* Calculate the structure of the IRE motif */
		IreStructure *ire_structure = predict_ire(ire_sequence);

		if(ire_structure->quality > 0) {
			compare_rnafold(ire_structure);
			if(search->cur_motif < 3) ire_structure->quality+=2;
		}

		/* Show header if IRE is found */
		if(ire_structure->quality > 0 && !header_shown && record->header) {
			output = strdup(record->header);
			header_shown = true;
		}

		/* Store current IRE if above threshold */
		if(ire_structure->quality > 0)	{
			char out[50] = {0};
			sprintf(out,"\n%s\t[ %0.2f ]\n",ire_structure->structure,ire_structure->quality);
			append(&output, ire_structure->sequence);
			append(&output, out);
		}

		/* Free resources */
		free(ire_sequence);
		irestruct_destroy(ire_structure);

		search->search_index += matcher.foundAtIndex + matcher.matchLength;
	} while(search->cur_motif < 18); // since there are only 18 motifs

	/* Output collected*/
	if(output) {
		THREADSAFE_STREAM_OUTPUT(fprintf(opt->out_file, "%s", output))
	}

	/* free stuff */
	free(record->sequence);
	if(record->header) free(record->header);
	free(record);
	free(search);
	free(output);
}


static Matcher
search_for_ire(SearchInfo *search, Regex *regex)
{
	Matcher matcher;
	if(search->search_index >= search->sequence_length) {
		matcher.isFound = false;
		return matcher;
	}
	
	matcher = regexMatch(&regex[search->cur_motif], search->sequence+search->search_index);
	return matcher;
}


static void
compare_rnafold(IreStructure *ire_structure)
{
	const char *sequence = ire_structure->sequence;
	char *structure = s_malloc(sizeof(char) * strlen(sequence) + 1);

	/* Get the fold compound from sequence */
	vrna_fold_compound_t *fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
	float best_mfe = vrna_mfe(fc, structure);

	/* Add structure constraints to fold */
	const char *constraint = get_contraints((const char *)ire_structure->structure);
	vrna_constraints_add(fc, constraint, VRNA_CONSTRAINT_DB_DEFAULT);

	/* Predict the minimum free energy & structure */
	float mfe = vrna_mfe(fc, structure);

	/* Update quality score based on mfe */
	float mfe_percent_change = fabsf(mfe - best_mfe) / fabsf((mfe + best_mfe)/2);
	if(mfe < 0) {
		ire_structure->quality += 1.5F*(1.0F - mfe_percent_change);
	}

	/* Update quality score based on structure similarities */
	ire_structure->quality += 1.5F - (float)levenshtein(ire_structure->structure, structure)/32.0F;

	/* Compare upper stem in ire_fold and vrna structures */
	ire_structure->quality += 2.5F * compare_structures(ire_structure->structure, structure);

	/* Cleanup */
	free((char *)constraint);
	free(structure);
	vrna_fold_compound_free(fc);
}


static float
compare_structures(const char *struct_iref, const char *struct_vrna)
{
	float score = 0;
	bool  cur_test = true;
	for(int cbulge=6; cbulge<9; cbulge++) {
		if(struct_iref[cbulge] != struct_vrna[cbulge]) cur_test = false;
	}
	if(cur_test) {
		score += 0.5f;
	}

	cur_test = true;
	for(int pair=13; pair<19; pair++) {
		if(struct_iref[pair] != struct_vrna[pair]) cur_test = false;
	}
	if(cur_test) {
		score += 0.25f;
	}

	cur_test = true;
	int left=12, right=19, count=0;
	while(count < 5) {
		if(struct_iref[left] != struct_vrna[left] || struct_iref[right] != struct_vrna[right]) {
			cur_test = false;
		}
		left--;
		right++;
		count++;
	}
	if(cur_test) {
		score += 0.25f;
	}

	return score;
}

/*==================== Helper functions ====================*/
static int
get_num_threads(struct aspiire_args_info args_info)
{
    int thread_max = max_user_threads();
    int num_threads = 0;

    if(args_info.jobs_given) {
        num_threads = MIN2(thread_max, args_info.jobs_arg);
    } else {
        /* use maximum of concurrent threads */
        int proc_cores, proc_cores_conf;
        if (num_proc_cores(&proc_cores, &proc_cores_conf)) {
            num_threads = MIN2(thread_max, proc_cores_conf);
        } else {
			warning_message("Could not determine number of available processor cores!\n"
			                "Defaulting to serial computation.");
            num_threads = 1;
        }
    }
    num_threads = MAX2(1, num_threads);

    return num_threads;
}


static Regex *
init_regex(void)
{
	const char *canonical_motif_1 = "C.....CAG[TUG]G[CUTAG]";
	const char *canonical_motif_2 = "C.....CAGAG[CATU]";
	const char *selex_motif_1 =     "C.....C[TU]G[TU]G[TUC]";
	const char *selex_motif_2 =     "C.....CCG[TU]G[ATUC]";
	const char *selex_motif_3 =     "C.....CCGAGA";
	const char *selex_motif_4 =     "C.....C[TU][TU]AGC";
	const char *selex_motif_5 =     "C.....CAA[TU]GC";
	const char *selex_motif_6 =     "C.....CAGGG[ACTUG]";
	const char *selex_motif_7 =     "C.....[TU]AG[TU]A[CTU]";
	const char *selex_motif_8 =     "C.....[TU]AGGA[TU]";
	const char *selex_motif_9 =     "C.....[TU]AGAA[TUC]";
	const char *selex_motif_10 =    "C.....[TU]AGCAG";
	const char *selex_motif_11 =    "C.....GAG[TU]C[GA]";
	const char *selex_motif_12 =    "C.....GAGCC[GA]";
	const char *selex_motif_13 =    "C.....GAGAG[TUG]";
	const char *selex_motif_14 =    "C.....GGGAG[CTUAG]";
	const char *selex_motif_15 =    "C.....GAG[TU]G[TUA]";
	const char *selex_motif_16 =    "G.....CAG[TU]GA";

	Regex *regex = s_malloc(18 * sizeof *regex);
	regexCompile(&regex[0], canonical_motif_1);
	regexCompile(&regex[1], canonical_motif_2);
	regexCompile(&regex[2], selex_motif_1);
	regexCompile(&regex[3], selex_motif_2);
	regexCompile(&regex[4], selex_motif_3);
	regexCompile(&regex[5], selex_motif_4);
	regexCompile(&regex[6], selex_motif_5);
	regexCompile(&regex[7], selex_motif_6);
	regexCompile(&regex[8], selex_motif_7);
	regexCompile(&regex[9], selex_motif_8);
	regexCompile(&regex[10], selex_motif_9);
	regexCompile(&regex[11], selex_motif_10);
	regexCompile(&regex[12], selex_motif_11);
	regexCompile(&regex[13], selex_motif_12);
	regexCompile(&regex[14], selex_motif_13);
	regexCompile(&regex[15], selex_motif_14);
	regexCompile(&regex[16], selex_motif_15);
	regexCompile(&regex[17], selex_motif_16);

	return regex;
}


RecordData *
copy_record_data(RnaInfo rna_info, Regex *regex, Options *opt) {
    RecordData *new_data = s_malloc(sizeof *new_data);
	
	new_data->sequence      = rna_info.sequence;
	new_data->header        = rna_info.header;
    new_data->regex         = regex;
	new_data->opt           = opt;

    return new_data;
}


static const char *
get_contraints(const char *predicted)
{
	char *constraint = strdup(predicted);
	constraint[7] = 'x';
	for(int i=13; i<19; i++) {
		constraint[i] = 'x';
	}

	return (const char *)constraint;
}


static char
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

static void
free_options(Options *opt)
{
	if(opt->input_file) {
		free(opt->input_file);
	}
	if(opt->out_given) {
		free(opt->out_filename);
	}
	if(opt->out_file) {
		fclose(opt->out_file);
	}
	if(opt->format_given) {
		free(opt->format);
	}
}

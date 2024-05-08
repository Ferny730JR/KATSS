#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

#include "ViennaRNA/fold.h"
#include "rna_file_parser.h"
#include "Regex.h"
#include "string_utils.h"
#include "utils.h"

#include "aspiire_cmdl.h"

#define BUFFER_SIZE 16777216


typedef struct SearchInfo {
	char *sequence;
	size_t sequence_length;
	size_t search_index;
	int cur_motif;
} SearchInfo;


typedef struct IreInfo {
	char *sequence;
	char *structure;
	char *gene_id;
	char *transcript_id;
	char *chromosome;

	uint loop_type;
	uint mismatches;
	uint bulges;
	char n25;
	uint wobble;
} IreInfo;


typedef struct IreStructure {
	char *sequence;
	char *structure;
	bool is_valid;

	bool paired;        // Determine if IRE is paired (idk tbh)
	uint no_pair;       // Number of times a nucleotide in the upper stem is unpaired
	uint lower_pairs;   // The number of nucleotide pairs in the lower stem
	int  first_pair;    // Position of the first pair
	uint second_pair;   // Position of the complementing nucleotide of the first pair
		
	uint mismatch_pair;

} IreStructure;


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
} Options;

/*==================== Function Declarations ====================*/
static void process_sequence(RnaInfo *rna_info, Regex *regex);
static Matcher search_for_ire(SearchInfo *search, Regex *regex);
static IreStructure *calculate_ire_structure(char *sequence);
static void find_ire(RnaInfo *rna_info, Regex *regex);
static void determine_UTRpair(IreStructure *ire);
static void lowerstem_UTRpair(IreStructure *ire, const uint unpaired_count);
static bool follow_constraints(uint no_pair, uint mismatch_pair);
static Regex *init_regex(void);
static char delimiter_to_char(char *user_delimiter);
static void free_options(Options *opt);


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

	while(true) {
		RnaInfo info = rnaf_geti(read_file);

		if(info.sequence == NULL) {
			break;
		}

		process_sequence(&info, regex);

		// find_ire(&info, regex);
		rnaf_destroy(&info);
	}
	rnaf_close(read_file);

	/* Clean up */
	free(regex);
	free_options(&opt);
	exit(EXIT_SUCCESS);
}

/*====================  Main functions  ====================*/
static void
process_sequence(RnaInfo *rna_info, Regex *regex)
{
	SearchInfo *search = s_malloc(sizeof *search);
	search->sequence = rna_info->sequence;
	search->sequence_length = strlen(rna_info->sequence);
	search->search_index = 0;
	search->cur_motif = 0;

	bool header_shown = false;

	do {
		Matcher matcher = search_for_ire(search, regex);
		if(!matcher.isFound) {
			search->search_index = 0;
			search->cur_motif++; /* Match not found for specified motif, search for next one */
			continue;
		}

		/* Get the IRE motif that was found */
		int ire_start = (search->search_index + matcher.foundAtIndex) - 6;
		if(ire_start < 0) { 
			search->search_index += matcher.foundAtIndex + matcher.matchLength;
			continue;
		}
		char *ire_sequence = substr(rna_info->sequence, ire_start, 31);

		/* Calculate the structure of the IRE motif */
		IreStructure *ire_structure = calculate_ire_structure(ire_sequence);

		/* TODO: Compare structure with RNAfold */

		/* Show header if IRE is found */
		if(ire_structure->is_valid && !header_shown && rna_info->header) {
			printf("%s",rna_info->header);
			header_shown = true;
		}

		/* Dump IRE information */
		if(ire_structure->is_valid) {
			printf("%s\n",ire_structure->sequence);
			printf("%s\n",ire_structure->structure);
		}

		/* Free resources */
		free(ire_structure->structure);
		free(ire_structure->sequence);
		free(ire_structure);

		search->search_index += matcher.foundAtIndex + matcher.matchLength;
	} while(search->cur_motif < 18); // since there are only 18 motifs

	/* free stuff */
	free(search);
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
find_ire(RnaInfo *rna_info, Regex *regex)
{
	size_t seq_len = strlen(rna_info->sequence);
	Matcher matcher;
	size_t search_index = 0;

	/* Search entire sequence for IREs using all motifs */
	for(int i=0; i<18; i++) {
		search_index = 0;
		while(search_index < seq_len) {
			matcher = regexMatch(&regex[i], rna_info->sequence+search_index);
			if(!matcher.isFound) { break; }

			/* Process match. Shift buffer by search_index due to search starting from offset*/
			int ire_start = (search_index + matcher.foundAtIndex) - 6;
			if(ire_start > 0) {
				char *ire_seq = substr(rna_info->sequence, ire_start, 31);
				IreStructure *ire = calculate_ire_structure(ire_seq);
				if(ire->is_valid) {
					printf("%s",rna_info->header);
					printf("%s\n",ire->sequence);
					printf("%s\n",ire->structure);
				}
				free(ire->structure);
				free(ire->sequence);
				free(ire);
			}
			search_index += matcher.foundAtIndex+matcher.matchLength;
		}
	}
}


static IreStructure *
calculate_ire_structure(char *sequence)
{
	IreStructure *ire = s_malloc(sizeof *ire);
	ire->paired = false;
	ire->no_pair = 0;
	ire->lower_pairs = 0;
	ire->first_pair = 11;
	ire->second_pair = 18;
	ire->mismatch_pair = 0;
	ire->is_valid = true;
	ire->sequence = sequence;
	ire->structure = strdup("......");

	/* Sequence is not a valid IRE, too small :( */
	if(strlen(sequence) != 31) {
		ire->is_valid = false;
		return ire;
	}

	/* Loops through all first_pair nucleotides going downstream */
	for(; 4<ire->first_pair && ire->second_pair<25; ire->first_pair-=1) {
			
		ire->lower_pairs=0;
		ire->paired=false;
		if(ire->first_pair == 6) { // Skips the C bulge in upper stem
			prepend(&ire->structure, ".");
			continue;
		}
			
		/* Used to compare first_pair with second_pair */
		while( ire->no_pair < 2 && !ire->paired ) {
			determine_UTRpair(ire);
			ire->second_pair+=1;
		}

		/* If the pairs do not match IRE pair constraints, then discard */
		if(!follow_constraints(ire->no_pair, ire->mismatch_pair)) {
			ire->is_valid = false;
			return ire;
		}
	}

	/* Check for pairs in the lower stem */
	while(ire->first_pair >= 0 && ire->second_pair < 31) {
		ire->paired = true;
		uint unp = 0;

		/* Used to compare first_pair with second_pair */
		while(ire->paired) {
			lowerstem_UTRpair(ire, unp);
			ire->second_pair++;
			unp++;
		}
		ire->first_pair--;
	}

	for(; ire->second_pair<31; ire->second_pair+=1) {
		append(&ire->structure, ".");
	}
	for(; ire->first_pair>=0; ire->first_pair--) {
		prepend(&ire->structure, ".");
	}

	if(ire->lower_pairs < 5) {
		ire->is_valid = false;
	}

	return ire;
}


static uint
check_pair(const char first_nucleotide, const char second_nucleotide)
{
	/* Define a mapping between nucleotide pairs and return values */
    const char         *pairs = "GC CG AT TA AU UA TG GT UG GU";
    const int return_values[] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
    
    /* Iterate over pairs and compare with the input nucleotides */
    for (int i = 0; i < 10; ++i) {
        if (pairs[i * 3] == first_nucleotide && pairs[i * 3 + 1] == second_nucleotide) {
            return return_values[i];
        }
    }
    
    return 0; /* If the pair is not found, return 0 */
}


static void // Determines if the given pair (i.e., GC) is a valid pair
determine_UTRpair(IreStructure *ire)
{
	const char first_nucleotide = ire->sequence[ire->first_pair];
	const char second_nucleotide = ire->sequence[ire->second_pair];

	const uint pair_type = check_pair(first_nucleotide, second_nucleotide);
	if(pair_type == 1) {
		ire->paired = true;
		append(&ire->structure, ")");
		prepend(&ire->structure, "(");

		if(ire->first_pair == 5) {
			ire->lower_pairs += 1;
		}

	} else if(pair_type == 2) {
		ire->paired = true;
		ire->mismatch_pair += 1;
		append(&ire->structure, "}");
		prepend(&ire->structure, "{");

		if(ire->first_pair == 5) {
			ire->lower_pairs+=1;
		}

	} else {
		if( 19 <= ire->second_pair && ire->second_pair <= 23) {
			ire->no_pair += 1;
			append(&ire->structure, ".");
		} else {
			ire->no_pair = 2;
		}
	}
}


/* Determines if the pair in the lower stem is a valid pair */
static void
lowerstem_UTRpair(IreStructure *ire, const uint unpaired_count)
{
	char first_nucleotide;
	char second_nucleotide;

	if(ire->first_pair >= 0) { first_nucleotide = ire->sequence[ire->first_pair]; }
	else { first_nucleotide = 'X'; }

	if(ire->second_pair < 31) { second_nucleotide = ire->sequence[ire->second_pair]; }
	else { second_nucleotide = 'X'; }


	const uint pair_type = check_pair(first_nucleotide, second_nucleotide);

	if(pair_type == 1) {
		if(unpaired_count==1) {
			prepend(&ire->structure, "(");
			append(&ire->structure, ".)");
			ire->paired = false;
			ire->lower_pairs += 1;
		} else if(unpaired_count==2) {
			prepend(&ire->structure, "(");
			append(&ire->structure, "..)");
			ire->paired = false;
			ire->lower_pairs += 1;
		} else {
			prepend(&ire->structure, "(");
			append(&ire->structure, ")");
			ire->paired = false;
			ire->lower_pairs += 1;
		}
	} else if(pair_type == 2) {
		if(unpaired_count==1) {
			prepend(&ire->structure, "{");
			append(&ire->structure, ".}");
			ire->paired = false;
			ire->lower_pairs += 1;
		} else if(unpaired_count==2) {
			prepend(&ire->structure, "{");
			append(&ire->structure, "..}");
			ire->paired = false;
			ire->lower_pairs += 1;
		} else {
			prepend(&ire->structure, "{");
			append(&ire->structure, "}");
			ire->paired = false;
			ire->lower_pairs += 1;
		}
	} else {
		if(unpaired_count > 1) {
			prepend(&ire->structure, ".");
			ire->paired = false;
			ire->second_pair -= 3;
		}
	}
}

/* Determines if the pairs that were determined follow IRE contraints */
static bool
follow_constraints(uint no_pair, uint mismatch_pair)
{
	if(no_pair==2) {
		return false;
	}
	if(mismatch_pair>2) {
		return false;
	}
	if(no_pair > 0 && mismatch_pair > 0) {
		return false;
	}
	return true;
}


/*==================== Helper functions ====================*/
static Regex *
init_regex(void)
{
	const char *canonical_motif_1 = "C.....CAGTG[CTAG]";
	const char *canonical_motif_2 = "C.....CAGAG[CAT]";
	const char *selex_motif_1 =     "C.....CTGTG[TC]";
	const char *selex_motif_2 =     "C.....CCGTG[ATC]";
	const char *selex_motif_3 =     "C.....CCGAGA";
	const char *selex_motif_4 =     "C.....CTTAGC";
	const char *selex_motif_5 =     "C.....CAATGC";
	const char *selex_motif_6 =     "C.....CAGGG[ACTG]";
	const char *selex_motif_7 =     "C.....TAGTA[CT]";
	const char *selex_motif_8 =     "C.....TAGGAT";
	const char *selex_motif_9 =     "C.....TAGAA[TC]";
	const char *selex_motif_10 =    "C.....TAGCAG";
	const char *selex_motif_11 =    "C.....GAGTC[GA]";
	const char *selex_motif_12 =    "C.....GAGCC[GA]";
	const char *selex_motif_13 =    "C.....GAGAG[TG]";
	const char *selex_motif_14 =    "C.....GGGAG[CTAG]";
	const char *selex_motif_15 =    "C.....GAGUG[TA]";
	const char *selex_motif_16 =    "G.....CAGTGA";

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

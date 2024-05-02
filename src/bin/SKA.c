#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <stdbool.h>
#include <errno.h>

#include "rna_file_parser.h"
#include "kmer_counter.h"
#include "Regex.h"
#include "string_utils.h"
#include "utils.h"

#include "SKA_cmdl.h"


typedef struct {
	KmerCounter *input_counter;
	KmerCounter *bound_counter;
} kcounts;


typedef struct IndependentProbCounts {
	KmerCounter *monomers;
	KmerCounter *dimers;
	KmerCounter *kmers;
} IndependentProbCounts;


typedef struct options {
	char    *input_file;
	char    *bound_file;
	char    *out_filename;
	FILE    *out_file;
	bool     out_given;
	int      kmer;
	int      iterations;
	char     file_delimiter;
	bool     no_log;

	bool     independent_probs;
	char   **fmotif;
	char    *motif;
	int      num_motifs;

	char  **top_kmer;
	int     cur_iter;

	unsigned long input_total;
	unsigned long bound_total;

	IndependentProbCounts *counts;
	kcounts                kmer_counts;

	bool encountered_error;
} options;


typedef struct TopKmer {
	char k_mer[16]; /* 16 because max 15mer + string terminator */
	double enrichment;
} TopKmer;


/*##########################################################
#  Function Declarations                                   #
##########################################################*/

/* SKA functions */
void 
process_iteration(options *opt);

KmerCounter *
count_kmers(char *filename, options *opt);

void
uncount_kmers(KmerCounter *counter, char *filename, options *opt);

TopKmer
get_top_kmer(KmerCounter *input_counter, KmerCounter *bound_counter);

IndependentProbCounts *
count_dimonomers(char *filename, options *opt);

void
uncount_dimonomers(char *filename, options *opt);

TopKmer
get_top_prediction(IndependentProbCounts *counts);

/* fmotif functions */
void
process_fmotifs(options *opt);

uint64_t
count_fmotifs(char *filename, char *fmotif);

/* motif pattern search functions */
void
process_motifs(options *opt);

RegexCluster *
process_motif(char *motif, options *opt);

RegexCluster *
count_motifs(char *filename, char *pattern);

void
process_motif_match(RegexCluster *cluster, Matcher matcher, char *search);

void
get_cluster_frequencies(RegexCluster *cluster);

void
get_bin_frequencies(RegexBins *bin);

void
get_pos_frequencies(BinsLenInfo *len, uint16_t num_pos);

double
get_pos_total(LenPosInfo pos);

RegexCluster *
get_cluster_enrichment(RegexCluster *input_cluster, RegexCluster *bound_cluster, char *pattern);

void
get_bin_enrichment(RegexBins *enr_bin, RegexBins in_bin, RegexBins bo_bin);

void
get_pos_enrichments(BinsLenInfo *enr_len, BinsLenInfo in_len, 
                    BinsLenInfo bo_len, uint16_t num_pos);

/* Helper functions */
char
delimiter_to_char(char *user_delimiter);

void
topkmer_to_file(TopKmer top_kmer, options *opt);

void
cluster_to_file(FILE *file, RegexCluster *cluster, char delimiter);

void
clusterlen_to_file(BinsLenInfo len, uint8_t bin, uint16_t lenlen, uint16_t 
                   maxlen, FILE *file, char delimiter);

void
free_options(options *opt);

void
print_options(options *opt);

void
print_cluster(RegexCluster *cluster);


void 
init_default_options(options *opt)
{
	opt->input_file     = NULL;
	opt->bound_file     = NULL;
	opt->out_filename   = "motif";
	opt->out_file       = NULL;
	opt->out_given      = false;
	opt->kmer           = 3;
	opt->iterations     = 1;
	opt->file_delimiter = ',';
	opt->no_log         = false;

	opt->independent_probs  = false;
	opt->fmotif             = NULL;
	opt->motif              = NULL;

	opt->top_kmer       = NULL;
	opt->cur_iter       = 0;

	opt->input_total    = 0;
	opt->bound_total    = 0;

	opt->kmer_counts.bound_counter = NULL;
	opt->kmer_counts.input_counter = NULL;
	opt->counts = NULL;

	opt->encountered_error = false;
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
	
	if(!args_info.input_given && !args_info.bound_given) {
		printf("Usage: SKA [OPTIONS] [<input.fa>] [<bound.fa>]\n");
		printf("Try 'SKA --help' for more information.\n");
		SKA_cmdline_parser_free(&args_info);
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
	} else {
		error_message("You need to provide a 'bound' file");
		SKA_cmdline_parser_free(&args_info);
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	if(args_info.output_given) {
		opt.out_given    = true;
		opt.out_filename = strdup(args_info.output_arg);
	}

	if(args_info.kmer_given) {
		if(args_info.kmer_arg <= 0) {
			error_message("option '--kmer=%d' must be a value greater than 0.",args_info.kmer_arg);
			SKA_cmdline_parser_free(&args_info);
			free_options(&opt);
			exit(EXIT_FAILURE);
		}
		if(args_info.kmer_arg > 15) {
			error_message("option '--kmer=%d' must be a value less than or equal to 15.",args_info.kmer_arg);
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

	if(args_info.no_log_given) {
		opt.no_log = true;
	}

	if(args_info.independent_probs_given) {
		opt.independent_probs = true;
	}

	if(args_info.fmotif_given) {
		opt.fmotif = s_malloc(args_info.fmotif_given * sizeof(char*));
		for(uint32_t i = 0; i < args_info.fmotif_given; i++) {
			opt.fmotif[i] = strdup(args_info.fmotif_arg[i]);
		}
		opt.num_motifs = args_info.fmotif_given;
	}

	if(args_info.motif_given) {
		opt.motif = strdup(args_info.motif_arg);
	}

	if(args_info.motif_given && args_info.fmotif_given) {
		error_message("--motif and --fmotif given! Please use only one.");
		SKA_cmdline_parser_free(&args_info);
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	char *filename = concat(opt.out_filename, ".dsv");
	opt.out_file = fopen(filename, "w");
	if (opt.out_file == NULL) {
		error_message("Could not write to file '%s'\n",filename);
		SKA_cmdline_parser_free(&args_info);
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	opt.top_kmer = s_malloc(opt.iterations * sizeof(char *));
	for(int i=0; i<opt.iterations; i++) {
		opt.top_kmer[i] = NULL;
	}

	if(!opt.independent_probs && !opt.input_file) {
		error_message("You need to provide an 'input' file");
		SKA_cmdline_parser_free(&args_info);
		free_options(&opt);
		exit(EXIT_FAILURE);
	}

	SKA_cmdline_parser_free(&args_info);
	free(filename);

	/*##########################################################
	#  Computations                                            #
	##########################################################*/

	if(opt.fmotif) {
		process_fmotifs(&opt);
	} else if (opt.motif) {
		process_motifs(&opt);
	} else {
		while(opt.cur_iter < opt.iterations) {
			process_iteration(&opt);
			if(opt.encountered_error) { break; }
		}
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
	TopKmer top_kmer;

	/* Compute independent SKA analysis (no control) */
	if(opt->independent_probs) {
		if(opt->cur_iter == 0) {
			IndependentProbCounts *counts = count_dimonomers(opt->bound_file, opt);
			opt->counts = counts;
		} else {
			uncount_dimonomers(opt->bound_file, opt);
		}
		top_kmer = get_top_prediction(opt->counts);

	/* Compute Iterative K-mer Knockout Analysis IKKA */
	} else {
		if(opt->cur_iter == 0) { // first iteration, count kmers
			opt->kmer_counts.input_counter = count_kmers(opt->input_file, opt);
			if(opt->kmer_counts.input_counter == NULL) { opt->encountered_error = true; return; }
			opt->kmer_counts.bound_counter = count_kmers(opt->bound_file, opt);
			if(opt->kmer_counts.bound_counter == NULL) { opt->encountered_error = true; return; }
		} else { // second iteration onwards, remove the previously counter kmers
			uncount_kmers(opt->kmer_counts.input_counter, opt->input_file, opt);
			uncount_kmers(opt->kmer_counts.bound_counter, opt->bound_file, opt);
		}

		/* After getting counts, we can get the top k-mer */
		top_kmer = get_top_kmer(opt->kmer_counts.input_counter,opt->kmer_counts.bound_counter);
	}

	/* Dump information into output file */
	topkmer_to_file(top_kmer, opt);

	/* Store information for next iteration */
	opt->top_kmer[opt->cur_iter] = strdup(top_kmer.k_mer);
	opt->cur_iter++;
}


KmerCounter *
count_kmers(char *filename, options *opt)
{
	RNA_FILE *read_file = rnaf_open(filename);
	if(read_file == NULL) {
		return NULL;
	}
	KmerCounter *counter = init_kcounter(opt->kmer);
	kctr_set_t_or_u(counter, read_file->is_t);

	char *sequence;
	while(true) {
		sequence = rnaf_get(read_file);

		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 0);

		kctr_fincrement(counter, sequence);
		free(sequence);
	}
	rnaf_close(read_file);

	return counter;
}


void
uncount_kmers(KmerCounter *counter, char *filename, options *opt)
{
	RNA_FILE *read_file = rnaf_open(filename);
	char     *top_kmer  = opt->top_kmer[opt->cur_iter-1];
	rnaf_rebuff(read_file, 65536); /* increase buffer size to make getm faster */

	char *sequence;
	while (true) {
		sequence = rnaf_getm(read_file, top_kmer);

		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 0);
		for(int i=0; i<opt->cur_iter-1; i++) {
			cross_out(sequence, opt->top_kmer[i]);
		}

		kctr_decrement(counter, sequence, top_kmer);
	}
	rnaf_close(read_file);
}


TopKmer
get_top_kmer(KmerCounter *input_counter, KmerCounter *bound_counter)
{
	TopKmer top_kmer = {.enrichment = -DBL_MAX};
	unsigned int top_enrichment_hash;
	double top_enrichment = -DBL_MAX;

	/* Sanity check, make sure total_count is greater than 0 */
	if(!input_counter->total_count || !bound_counter->total_count) {
		return top_kmer;
	}
	
	for(unsigned int i=0; i<input_counter->capacity; i++) {
		/* Get frequencies of input and bound */
		double input_frq = (double)input_counter->entries[i]/input_counter->total_count;
		double bound_frq = (double)bound_counter->entries[i]/bound_counter->total_count;

		/* if input_frq is 0, then div by 0 error would occur so skip */
		if(input_frq == 0) {
			continue;
		}

		/* get enrichment of current kmer */
		double cur_enrichment = log2(bound_frq/input_frq);

		if(cur_enrichment > top_enrichment) {
			top_enrichment = cur_enrichment;
			top_enrichment_hash = i;
		}
	}

	/* No top kmer was found (probs because something terrible happened) return empty struct */
	if(top_enrichment == top_kmer.enrichment) {
		return top_kmer;
	}

	/* Get k-mer string from hash */
	kctr_get_key(input_counter, top_kmer.k_mer, top_enrichment_hash);

	/* Fill struct with info */
	top_kmer.enrichment = top_enrichment;

	/* Everything (probably) worked! Hurray, now return. */
	return top_kmer;
}


IndependentProbCounts *
count_dimonomers(char *filename, options *opt)
{
	RNA_FILE *read_file = rnaf_open(filename);
	if(read_file == NULL) {
		return NULL;
	}

	IndependentProbCounts *counts = s_malloc(sizeof(IndependentProbCounts));
	counts->monomers = init_kcounter(1U);     kctr_set_t_or_u(counts->monomers, read_file->is_t);
	counts->dimers   = init_kcounter(2U);     kctr_set_t_or_u(counts->dimers, read_file->is_t);
	counts->kmers = init_kcounter(opt->kmer); kctr_set_t_or_u(counts->kmers, read_file->is_t);

	char *sequence;
	while(true) {
		sequence = rnaf_get(read_file);

		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 0);
		kctr_fincrement(counts->monomers, sequence);
		kctr_fincrement(counts->dimers, sequence);
		kctr_fincrement(counts->kmers, sequence);

		free(sequence);
	}
	rnaf_close(read_file);

	return counts;
}


void
uncount_dimonomers(char *filename, options *opt)
{
	RNA_FILE *read_file = rnaf_open(filename);
	if(read_file == NULL) {
		opt->encountered_error = true;
		return;
	}

	char *top_kmer = opt->top_kmer[opt->cur_iter-1];
	rnaf_rebuff(read_file, 65536U);

	char *sequence;
	while(true) {
		sequence = rnaf_getm(read_file, top_kmer);

		if(sequence == NULL) {
			break;
		}

		clean_seq(sequence, 0);
		for(int i=0; i<opt->cur_iter-1; i++) {
			cross_out(sequence, opt->top_kmer[i]);
		}

		kctr_decrement(opt->counts->kmers, sequence, top_kmer);
	}
	rnaf_close(read_file);
}


double
predict_kmer(char *kseq, KmerCounter *monomer_counts, KmerCounter *dimer_counts)
{
	char monoseq[2];  monoseq[1] = '\0';
	char diseq[3];    diseq[2]   = '\0';

	double monoprob = 1;
	double diprob = 1;

	int kmer = strlen(kseq);

	/* Get the cumulative probabilities for overlapping monomers */
	for(int i=1; i<kmer-1; i++) {
		monoseq[0] = kseq[i];
		double count = (double)kctr_get(monomer_counts, monoseq);
		monoprob *= count/monomer_counts->total_count;
	}

	/* Get probabilities for all di-mers in k-mer */
	for(int i=0; i<kmer-1; i++) {
		diseq[0] = kseq[i]; diseq[1] = kseq[i+1];
		double count = (double)kctr_get(dimer_counts, diseq);
		diprob *= count/dimer_counts->total_count;
	}

	/* Predicted k-mer probability is dinucleotides / overlapping monomers */
	return diprob/monoprob;
}


TopKmer
get_top_prediction(IndependentProbCounts *counts)
{
	TopKmer top_kmer = {.enrichment = -DBL_MAX};
	double top_enrichment = -DBL_MAX;
	char kseq[16];

	for(unsigned int i=0; i<counts->kmers->capacity; i++) {
		kctr_get_key(counts->kmers, kseq, i);

		/* Get actual and predicted frequencies */
		double kmer_frq = (double)counts->kmers->entries[i]/counts->kmers->total_count;
		double predicted_frq = predict_kmer(top_kmer.k_mer, counts->monomers, counts->dimers);

		/* if input_frq is 0, then div by 0 error would occur so skip */
		if(predicted_frq == 0) {
			continue;
		}

		/* get enrichment of current kmer */
		double cur_enrichment = log2(kmer_frq/predicted_frq);
		// printf("%s: log2(%f/%f) = %f\n",top_kmer.k_mer, kmer_frq,predicted_frq,cur_enrichment);

		if(cur_enrichment > top_enrichment) {
			top_enrichment = cur_enrichment;
			kctr_get_key(counts->kmers, top_kmer.k_mer, i);
		}
	}

	/* No top kmer was found (probs because something terrible happened) return empty struct */
	if(top_enrichment == top_kmer.enrichment) {
		return top_kmer;
	}

	/* Fill struct with info */
	top_kmer.enrichment = top_enrichment;

	/* Everything (probably) worked! Hurray, now return. */
	return top_kmer;
}


void
process_fmotifs(options *opt)
{
	RNA_FILE *input_file = rnaf_open(opt->input_file);
	RNA_FILE *bound_file = rnaf_open(opt->bound_file);

	uint64_t input_numlines = rnaf_numlines(input_file);
	uint64_t bound_numlines = rnaf_numlines(input_file);

	uint64_t input_total = rnaf_numchars(input_file);
	uint64_t bound_total = rnaf_numchars(bound_file);

	rnaf_close(input_file);
	rnaf_close(bound_file);

	for(int i = 0; i < opt->num_motifs; i++) {
		uint64_t input_counts = count_fmotifs(opt->input_file, opt->fmotif[i]);
		input_total -= (input_numlines * strlen(opt->fmotif[i]));
		double input_frq = (double)input_counts / input_total;

		uint64_t bound_counts = count_fmotifs(opt->bound_file, opt->fmotif[i]);
		bound_total -= (bound_numlines * strlen(opt->fmotif[i]));
		double bound_frq = (double)bound_counts / bound_total;

		fprintf(opt->out_file, "%s,%f", opt->fmotif[i], log2(bound_frq/input_frq));
	}
}


uint64_t
count_fmotifs(char *filename, char *fmotif)
{	/* Open file to read */
	RNA_FILE *read_file = rnaf_open(filename);
	rnaf_rebuff(read_file, 65536);

	uint64_t total_matches = 0;
	size_t fmotif_len = strlen(fmotif) - 1;

	/* Loop through entire file and search for fmotif */
	while(rnaf_oread(read_file, fmotif_len)) {
		total_matches += rnaf_search(read_file, fmotif);
	}
	rnaf_close(read_file);

	return total_matches;
}


void
process_motifs(options *opt)
{
	RegexCluster *enrichment;

	/* Check if motif is a file */
	if(access(opt->motif, F_OK) != 0) {
		enrichment = process_motif(opt->motif, opt);
		if(enrichment == NULL) {
			return;
		}

		fprintf(opt->out_file, "%s%c%f\n", opt->motif, opt->file_delimiter, enrichment->total);
		cluster_to_file(opt->out_file, enrichment, opt->file_delimiter);
		return;
	}

	/* Argument passed is a file */
	FILE *fp = fopen(opt->motif, "r");
	if(fp == NULL) {
		error_message("Unable to open file '%s' for reading: %s",opt->motif, strerror(errno));
		return;
	}

	char buf[1024];
	while(fgets(buf, 1024, fp)) {
		remove_escapes(buf);

		/* Calculate the cluster enrichments of the motif */
		RegexCluster *enrichment = process_motif(buf, opt);
		if(enrichment == NULL) {
			return;
		}

		/* Output the cluster enrichments to file */
		fprintf(opt->out_file, "%s%c%f\n", buf, opt->file_delimiter, enrichment->total);
		cluster_to_file(opt->out_file, enrichment, opt->file_delimiter);
		fprintf(opt->out_file, "\n");
	}
}


RegexCluster *
process_motif(char *motif, options *opt)
{
	/* Get the counts and frequencies of the motif matches in the input file */
	RegexCluster *input = count_motifs(opt->input_file, motif);
	if(input) {
		get_cluster_frequencies(input);
	} else {
		return NULL; // failed to open file, so break
	}

	/* Get the counts and frequencies of the motif matches in the bound file */
	RegexCluster *bound = count_motifs(opt->bound_file, motif);
	if(bound) {
		get_cluster_frequencies(bound);
	} else {
		freeRegexCluster(input);
		return NULL; // failed to open file, so break
	}

	RegexCluster *enrichments = get_cluster_enrichment(input, bound, motif);
	
	freeRegexCluster(input);
	freeRegexCluster(bound);
	
	return enrichments;
}


RegexCluster *
count_motifs(char *filename, char *pattern)
{	/* Open file to read */
	RNA_FILE *read_file = rnaf_open(filename);
	if(read_file == NULL) {
		return NULL;
	}

	/* Increase the buffer size */
	rnaf_rebuff(read_file, 65536);

	/* Make regex pattern */
	Regex regex;
	regexCompile(&regex, pattern);
	if(!regex.isPatternValid) {
		error_message("Could not compile '%s': %s", pattern, regex.errorMessage);
		rnaf_close(read_file);
		return NULL;
	}

	/* Make global cluster */
	RegexCluster *cluster = regexClusterInit(&regex);

	/* Search for matches */
	Matcher matcher;
	while(rnaf_oread(read_file, strlen(pattern))) {
		uint32_t search_index = 0;
		while(search_index < read_file->buffer_size) {
			matcher = regexMatch(&regex, read_file->buffer+search_index);
			if(!matcher.isFound) { break; }

			/* Process match. Shift buffer by search_index due to search starting from offset*/
			process_motif_match(cluster, matcher, read_file->buffer+search_index);
			search_index += matcher.foundAtIndex+matcher.matchLength;
		}
	}
	rnaf_close(read_file);

	/* Return cluster */
	return cluster;
}


void
process_motif_match(RegexCluster *cluster, Matcher matcher, char *search)
{
	cluster->total++;

	uint8_t clusterIndex = 0;
	while(matcher.cluster[clusterIndex].binType != BIN_END_OF_PATTERN) {
		int32_t binClusterLength = matcher.cluster[clusterIndex].clusterLength;

		for(int32_t index = 0; index < binClusterLength; index++) {
			int32_t binStartIndex = matcher.cluster[clusterIndex].startIndex + index;
			reclustAddNT(cluster, search[binStartIndex], clusterIndex, binClusterLength, index);
		}
		clusterIndex++;
	}
}


void
get_cluster_frequencies(RegexCluster *cluster)
{
	for(uint8_t num_bin = 0; num_bin < cluster->num_bins; num_bin++) {
		get_bin_frequencies(&cluster->bin[num_bin]);
	}
}


void
get_bin_frequencies(RegexBins *bin)
{
	uint32_t num_lens = (bin->maxLen - bin->minLen) + 1;
	for(uint32_t cur_len = 0; cur_len < num_lens; cur_len++) {
		/* Get the frequency of the current length */
		bin->len[cur_len].total /= bin->total;

		/* Get the frequencies of the NT in the current length*/
		uint16_t lenlen = bin->minLen + cur_len;
		get_pos_frequencies(&bin->len[cur_len], lenlen);
	}
}


void
get_pos_frequencies(BinsLenInfo *len, uint16_t num_pos)
{
	for(uint16_t cur_pos = 0; cur_pos < num_pos; cur_pos++) {
		double nt_total = get_pos_total(len->pos[cur_pos]);
		len->pos[cur_pos].A /= nt_total;
		len->pos[cur_pos].C /= nt_total;
		len->pos[cur_pos].G /= nt_total;
		len->pos[cur_pos].T /= nt_total;
	}
}


double
get_pos_total(LenPosInfo pos)
{
	return pos.A + pos.C + pos.G + pos.T;
}


RegexCluster *
get_cluster_enrichment(RegexCluster *input_cluster, RegexCluster *bound_cluster, char *pattern)
{	/* Create the enrichment cluster */
	Regex regex;
	regexCompile(&regex, pattern);
	RegexCluster *enrichment_cluster = regexClusterInit(&regex);

	/* Get overall enrichment */
	enrichment_cluster->total = bound_cluster->total / input_cluster->total;

	/* Get the enrichments for the bins */
	for(uint8_t cur_bin = 0; cur_bin < enrichment_cluster->num_bins; cur_bin++) {
		get_bin_enrichment(&enrichment_cluster->bin[cur_bin],
		                   input_cluster->bin[cur_bin],
		                   bound_cluster->bin[cur_bin]);
	}

	return enrichment_cluster;
}


void
get_bin_enrichment(RegexBins *enr_bin, RegexBins in_bin, RegexBins bo_bin)
{
	uint32_t num_lens = (enr_bin->maxLen - enr_bin->minLen) + 1;
	for(uint32_t cur_len = 0; cur_len < num_lens; cur_len++) {
		enr_bin->len[cur_len].total = bo_bin.len[cur_len].total / in_bin.len[cur_len].total;

		uint16_t lenlen = enr_bin->minLen + cur_len;
		get_pos_enrichments(&enr_bin->len[cur_len], in_bin.len[cur_len], 
		                     bo_bin.len[cur_len], lenlen);
	}
}


void
get_pos_enrichments(BinsLenInfo *enr_len, BinsLenInfo in_len, BinsLenInfo bo_len, uint16_t num_pos)
{
	for(uint16_t cur_pos = 0; cur_pos < num_pos; cur_pos++) {
		enr_len->pos[cur_pos].A = bo_len.pos[cur_pos].A / in_len.pos[cur_pos].A;
		enr_len->pos[cur_pos].C = bo_len.pos[cur_pos].C / in_len.pos[cur_pos].C;
		enr_len->pos[cur_pos].G = bo_len.pos[cur_pos].G / in_len.pos[cur_pos].G;
		enr_len->pos[cur_pos].T = bo_len.pos[cur_pos].T / in_len.pos[cur_pos].T;
	}
}
/*##########################################################
#  Helper Functions                                        #
##########################################################*/

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
topkmer_to_file(TopKmer top_kmer, options *opt)
{
	char *kseq = strdup(top_kmer.k_mer);
	clean_seq(kseq, true);
	fprintf(opt->out_file, "%s%c%f\n", kseq, opt->file_delimiter, top_kmer.enrichment);
	free(kseq);
}


void
cluster_to_file(FILE *file, RegexCluster *cluster, char delimiter)
{
	/* Get the max len */
	uint16_t cur, max_len = 0;
	for(uint8_t cur_bin = 0; cur_bin < cluster->num_bins; cur_bin++)
	{
		cur = (cluster->bin[cur_bin].maxLen - cluster->bin[cur_bin].minLen) + 1;
		if(cur>max_len) {
			max_len = cur;
		}
	}

	/* Print the header of CSV file */
	fprintf(file, "bin%clength%clength_pref", delimiter, delimiter);
	for(uint16_t i = 0; i < max_len; i++) {
		fprintf(file, "%cNT%d_pref", delimiter, i+1);
	}
	fprintf(file, "\n");

	/* Output cluster to file */
	for(uint8_t cur_bin = 0; cur_bin < cluster->num_bins; cur_bin++) {
		uint16_t num_lens = (cluster->bin[cur_bin].maxLen - cluster->bin[cur_bin].minLen) + 1;
		for(uint16_t cur_len = 0; cur_len < num_lens; cur_len++) {
			uint16_t lenlen = cur_len + cluster->bin[cur_bin].minLen;
			clusterlen_to_file(cluster->bin[cur_bin].len[cur_len], cur_bin, lenlen, max_len, file, delimiter);
		}
	}
}


void
clusterlen_to_file(BinsLenInfo len, uint8_t bin, uint16_t lenlen, uint16_t 
                   maxlen, FILE *file, char delimiter)
{
	fprintf(file, "%d,%d,%f",bin+1,lenlen,len.total);
	for(uint16_t cur = 0; cur < maxlen; cur++) {
		if(cur < lenlen) {
			fprintf(file, "%c\"A:%f", delimiter, len.pos[cur].A);
			fprintf(file, "%cC:%f", delimiter, len.pos[cur].C);
			fprintf(file, "%cG:%f", delimiter, len.pos[cur].G);
			fprintf(file, "%cT:%f\"", delimiter, len.pos[cur].T);
		} else {
			fprintf(file, "%c-", delimiter);
		}
	}
	fprintf(file, "\n");
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
		if(opt->top_kmer) {
			free(opt->top_kmer[i]);
		}
	}
	if(opt->motif) {
		free(opt->motif);
	}
	if(opt->fmotif) {
		for(int i=0; i < opt->num_motifs; i++) {
			free(opt->fmotif[i]);
		}
		free(opt->fmotif);
	}
	if(opt->top_kmer) {
		free(opt->top_kmer);
	}
	if(opt->out_file) {
		fclose(opt->out_file);
	}
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


void
print_cluster(RegexCluster *cluster)
{
	for(int i=0; i<cluster->num_bins; i++) {
		printf("bins[%d] = { max: %d, min: %d, total: %f\n",i+1,cluster->bin[i].maxLen, cluster->bin[i].minLen,cluster->bin[i].total);
		int num_lens = (cluster->bin[i].maxLen - cluster->bin[i].minLen) + 1;
		for(int len=0; len < num_lens; len++) {
			int lenlen = len + cluster->bin[i].minLen;
			printf("    len[%d] = {\n        .total=%f,\n",lenlen,cluster->bin[i].len[len].total);
			for(int pos=0; pos < lenlen; pos++) {
				printf("        .pos[%d] = {",pos+1);
				printf(".A=%f ",cluster->bin[i].len[len].pos[pos].A);
				printf(".C=%f ",cluster->bin[i].len[len].pos[pos].C);
				printf(".G=%f ",cluster->bin[i].len[len].pos[pos].G);
				printf(".T=%f ",cluster->bin[i].len[len].pos[pos].T);
				printf("}\n");
			}
			printf("    }\n");
		}
		printf("}\n");
	}
}
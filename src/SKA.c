#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "kmerHashTable.h"
#include "string_utils.h"
#include "utils.h"

#include "SKA_cmdl.h"

typedef struct {
    char    *input_file;
    char    *bound_file;
    char    *out_file;
    int     out_given;
    int     kmer;
    char    file_delimiter;

    int     independent_probs;
    int     frq;
} options;


typedef struct {
    kmerHashTable *monomer_frq;
    kmerHashTable *dimer_frq;
    kmerHashTable *kmer_frq;
} frqIndependentProbs;


kmerHashTable *count_kmers(char *filename, options *opt);


void process_counts(char *line, kmerHashTable *counts_table, int kmer);


frqIndependentProbs process_independent_probs(char *filename, options *opt);


void count_di_mono_nt(char *sequence, frqIndependentProbs kmer_data, options *opt);


kmerHashTable *predict_kmers(kmerHashTable *probs_1mer, kmerHashTable *probs_2mer, int kmer);


void getFrequencies(kmerHashTable *counts_table);


kmerHashTable *getEnrichment(kmerHashTable *input_frq, kmerHashTable *bound_frq, options *opt);


char delimiter_to_char(char *user_delimiter);


int compare(const void *a, const void *b);


void free_options(options *opt);


void print_options(options *opt);


void init_default_options(options *opt) {
    opt->input_file     = NULL;
    opt->bound_file     = NULL;
    opt->out_file       = "motif";
    opt->out_given      = 0;
    opt->kmer           = 3;
    opt->file_delimiter = ',';

    opt->independent_probs  = 0;
    opt->frq                = 0;
}


/*##########################################################
#  Main                                                    #
##########################################################*/
int main(int argc, char **argv) {
    struct SKA_args_info    args_info;
    options                 opt;

    kmerHashTable           *input_table;
    kmerHashTable           *bound_table;
    kmerHashTable           *enrichments_table;

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
        opt.out_given = 1;
        opt.out_file  = strdup(args_info.output_arg);
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

    if(args_info.file_delimiter_given) {
        opt.file_delimiter = delimiter_to_char(args_info.file_delimiter_arg);
    }

    if(args_info.independent_probs_given) {
        opt.independent_probs = 1;
    }

    SKA_cmdline_parser_free(&args_info);

    /*##########################################################
    #  Computations                                            #
    ##########################################################*/

    if(opt.independent_probs) {
        frqIndependentProbs kmer_data;
        kmer_data = process_independent_probs(opt.bound_file, &opt);

        input_table = predict_kmers(kmer_data.monomer_frq, kmer_data.dimer_frq, opt.kmer);
        bound_table = kmer_data.kmer_frq;
    } else {
        input_table = count_kmers(opt.input_file, &opt);
        bound_table = count_kmers(opt.bound_file, &opt);
        getFrequencies(input_table);
        getFrequencies(bound_table);
    }

    enrichments_table = getEnrichment(input_table, bound_table, &opt);

    qsort(enrichments_table->entries, enrichments_table->capacity,
          sizeof(*enrichments_table->entries), compare);

    kmerHashTable_to_file(enrichments_table, opt.out_file, opt.file_delimiter);

    /* Clean up */
    free_options(&opt);
    free_kmer_table(input_table);
    free_kmer_table(bound_table);
    free_kmer_table(enrichments_table);
}


kmerHashTable *count_kmers(char *filename, options *opt) {
    FILE *read_file;
    kmerHashTable *counts_table;

    counts_table = init_kmer_table(opt->kmer, 1);
    read_file    = fopen(filename, "r");

    char buffer[10000];
    while(fgets(buffer, sizeof(buffer), read_file)) {
        seq_to_RNA(buffer);
        str_to_upper(buffer);
        remove_escapes(buffer);

        process_counts(buffer, counts_table, opt->kmer);
    }
    fclose(read_file);

    return counts_table;
}


void process_counts(char *sequence, kmerHashTable *counts_table, int kmer) {
    char *k_substr;
    int seq_length = strlen(sequence);
    int num_kmers_in_seq = seq_length - kmer + 1;

    for(int i=0; i<num_kmers_in_seq; i++) {
        // Get kmer substring
        k_substr = substr(sequence, i, kmer);

        kmer_add_value(counts_table, k_substr, 1, 0);

        free(k_substr);
    }
}


frqIndependentProbs process_independent_probs(char *filename, options *opt) {
    FILE                    *read_file;
    kmerHashTable           *monomers_table;
    kmerHashTable           *dimers_table;
    kmerHashTable           *counts_table;
    frqIndependentProbs     kmer_data;

    read_file = fopen(filename, "r");

    monomers_table  = init_kmer_table(1,1);
    dimers_table    = init_kmer_table(2,1);
    counts_table    = init_kmer_table(opt->kmer, 1);

    kmer_data.monomer_frq = monomers_table;
    kmer_data.dimer_frq   = dimers_table;
    kmer_data.kmer_frq    = counts_table;

    char buffer[10000];
    while(fgets(buffer, sizeof(buffer), read_file)) {
        str_to_upper(buffer);
        seq_to_RNA(buffer);
        remove_escapes(buffer);

        count_di_mono_nt(buffer, kmer_data, opt);
    }
    fclose(read_file);

    getFrequencies(kmer_data.monomer_frq);
    getFrequencies(kmer_data.dimer_frq);
    getFrequencies(kmer_data.kmer_frq);

    return kmer_data;
}


void count_di_mono_nt(char *sequence, frqIndependentProbs kmer_data, options *opt) {
    char    *k_substr;
    int     seq_length = strlen(sequence);
    int     num_di_nt_in_seq = seq_length - 1;
    int     num_kmers_in_seq = seq_length - opt->kmer + 1;

    for(int i=0; i<seq_length; i++) {
        // count monomers
        k_substr = substr(sequence, i, 1);
        kmer_add_value(kmer_data.monomer_frq, k_substr, 1, 0);
        free(k_substr);

        // count dimers
        if(i<seq_length-1) {
            k_substr = substr(sequence, i, 2);
            kmer_add_value(kmer_data.dimer_frq, k_substr, 1, 0);
            free(k_substr);
        }

        // counts kmers
        if(i<num_kmers_in_seq) {
            k_substr = substr(sequence, i, opt->kmer);
            kmer_add_value(kmer_data.kmer_frq, k_substr, 1, 0);
            free(k_substr);
        }
    }
}


kmerHashTable *predict_kmers(kmerHashTable *probs_1mer, kmerHashTable *probs_2mer, int kmer) {
    kmerHashTable   *predicted_kmers;
    double          dinucleotides_prob;
    double          monomers_probs;
    double          *prob;
    char            **kmers;
    char            *k_substr;

    predicted_kmers = init_kmer_table(kmer,1);
    kmers = s_malloc(predicted_kmers->capacity * sizeof(char *));

    /* Produce array of strings containing all k-mers */
    for (int i = 0; i < predicted_kmers->capacity; i++) {
        kmers[i] = malloc((kmer + 1) * sizeof(char));

        int index = i;
        for (int j = kmer - 1; j >= 0; j--) {
            kmers[i][j] = BASES[index % 4];
            index /= 4;
        }
        kmers[i][kmer] = '\0';
    }

    /* Get the probability of kmer being in reads */
    for(int i = 0; i < predicted_kmers->capacity; i++) {
        dinucleotides_prob = 1;
        monomers_probs     = 1;
        for(int j=0; j < kmer-1; j++) {
            k_substr = substr(kmers[i], j, 2);
            prob = kmer_get(probs_2mer, k_substr);
            free(k_substr);
            dinucleotides_prob*=prob[0];
        }
        for(int j=1; j < kmer-1; j++) {
            k_substr = substr(kmers[i], j, 1);
            prob = kmer_get(probs_1mer, k_substr);
            free(k_substr);
            monomers_probs*=prob[0];
        }
        kmer_add_value(predicted_kmers, kmers[i], dinucleotides_prob/monomers_probs, 0);
        free(kmers[i]);
    }
    free(kmers);

    return predicted_kmers;
}


void getFrequencies(kmerHashTable *counts_table) {
    int num_columns = counts_table->cols-1;
    double total_count = 0;
    double k_count;

    for(int i=0; i<counts_table->capacity; i++) {
        total_count+=counts_table->entries[i]->values[0];
    }

    for(int i=0; i<counts_table->capacity; i++) {
        k_count = counts_table->entries[i]->values[0];
        counts_table->entries[i]->values[num_columns] = k_count/total_count;
    }
}


kmerHashTable *getEnrichment(kmerHashTable *input_frq, kmerHashTable *bound_frq, options *opt) {
    kmerHashTable   *enrichments_table;
    char            *key;
    double          *input_values;
    double          *bound_values;
    double          enrichment;

    enrichments_table = init_kmer_table(opt->kmer, 1);
    
    // Get the log2 fold change for each kmer
    for(size_t i = 0; i < enrichments_table->capacity; i++) {
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

char delimiter_to_char(char *user_delimiter) {
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


int compare(const void *a, const void *b) {
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

    if(entryA->values[0] > entryB->values[0]) {
        return -1;
    }
    if(entryA->values[0] < entryB->values[0]) {
        return 1;
    }

    return 0;
}


void free_options(options *opt) {
    if(opt->input_file && !opt->independent_probs) {
        free(opt->input_file);
    }
    if(opt->bound_file) {
        free(opt->bound_file);
    }
    if(opt->out_given) {
        free(opt->out_file);
    }
}


void print_options(options *opt) {
    printf("input_file: '%s'\n",opt->input_file);
    printf("bound_file: '%s'\n",opt->bound_file);
    printf("output_file: '%s'\n",opt->out_file);
    printf("kmer: '%d'\n",opt->kmer);
    printf("frq: '%d'\n",opt->frq);
    printf("probs: '%d'\n",opt->independent_probs);
}
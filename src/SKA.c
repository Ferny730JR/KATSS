#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

#include "rna_file_parser.h"
#include "kmerHashTable.h"
#include "string_utils.h"
#include "utils.h"
#include "parallel_helpers.h"

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
    int     jobs;

    char    **top_kmer;
    int     cur_iter;
} options;


typedef struct {
    kmerHashTable *monomer_frq;
    kmerHashTable *dimer_frq;
    kmerHashTable *kmer_frq;
} frqIndependentProbs;


typedef struct {
    kmerHashTable   *counts_table;
    char            *sequence;
    int             kmer;
} record_data;


kmerHashTable *count_kmers(char *filename, options *opt);


void process_counts(record_data *record);


frqIndependentProbs process_independent_probs(char *filename, options *opt);


void count_di_mono_nt(char *sequence, frqIndependentProbs kmer_data, options *opt);


kmerHashTable *predict_kmers(kmerHashTable *probs_1mer, kmerHashTable *probs_2mer, int kmer);


void getFrequencies(kmerHashTable *counts_table);


kmerHashTable *getEnrichment(kmerHashTable *input_frq, kmerHashTable *bound_frq, options *opt);


Entry *kmer_max_entry(kmerHashTable *hash_table);


char delimiter_to_char(char *user_delimiter);


int compare(const void *a, const void *b);


void entry_to_file(FILE *file, Entry *entry, char delimiter);


void free_options(options *opt);


void print_options(options *opt);


void init_default_options(options *opt) {
    opt->input_file     = NULL;
    opt->bound_file     = NULL;
    opt->out_filename   = "motif";
    opt->out_file       = NULL;
    opt->out_given      = 0;
    opt->kmer           = 3;
    opt->iterations     = 1;
    opt->file_delimiter = ',';

    opt->independent_probs  = 0;
    opt->jobs               = 0;

    opt->top_kmer       = NULL;
    opt->cur_iter       = 0;
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
    }

    if(args_info.file_delimiter_given) {
        opt.file_delimiter = delimiter_to_char(args_info.file_delimiter_arg);
    }

    if(args_info.independent_probs_given) {
        opt.independent_probs = 1;
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

    char *filename = concat(opt.out_filename, ".dsv");
    opt.out_file = fopen(filename, "w");
    if (opt.out_file == NULL) {
        error_message("Could not write to file '%s'\n",filename);
    }

    opt.top_kmer = s_malloc(opt.iterations * sizeof(char *));
    for(size_t i=0; i<opt.iterations; i++) {
        opt.top_kmer[i] = NULL;
    }

    SKA_cmdline_parser_free(&args_info);
    free(filename);

    /*##########################################################
    #  Computations                                            #
    ##########################################################*/

    while(opt.cur_iter < opt.iterations) {
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

        Entry *max_entry = kmer_max_entry(enrichments_table);
        entry_to_file(opt.out_file, max_entry, opt.file_delimiter);

        opt.top_kmer[opt.cur_iter] = strdup(max_entry->key);
        opt.cur_iter++;

        /* free tables */
        free_kmer_table(enrichments_table);
        free_kmer_table(input_table);
        free_kmer_table(bound_table);
    }

    /* Clean up */
    free_options(&opt);

    return 0;
}


kmerHashTable *count_kmers(char *filename, options *opt) {
    RNA_FILE *read_file;
    kmerHashTable *counts_table;
    record_data *record = s_malloc(sizeof *record);

    counts_table = init_kmer_table(opt->kmer, 1);
    read_file    = rnaf_open(filename);

    record->counts_table = counts_table;
    record->kmer = opt->kmer;

    char *sequence;
    while(1) {
        sequence = rnaf_get(read_file);

        if(sequence == NULL) {
            break;
        }

        seq_to_RNA(sequence);
        str_to_upper(sequence);
        remove_escapes(sequence);
        for(size_t i=0; i<opt->cur_iter; i++) {
            cross_out(sequence, opt->top_kmer[i]);
        }

        record->sequence = sequence;
        
        process_counts(record);
    }
    rnaf_close(read_file);
    free(record);

    return counts_table;
}


void process_counts(record_data *record) {
    char *k_substr;
    int seq_length = strlen(record->sequence);
    int num_kmers_in_seq = seq_length - record->kmer + 1;

    for(int i=0; i<num_kmers_in_seq; i++) {
        // Get kmer substring
        k_substr = substr(record->sequence, i, record->kmer);

        if(k_substr[record->kmer-1] != 'X' && k_substr[0] != 'X') {
            kmer_add_value(record->counts_table, k_substr, 1, 0);
        }

        free(k_substr);
    }

    free(record->sequence);
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
        if(counts_table->entries[i] == NULL) {
            continue;
        }
        total_count+=counts_table->entries[i]->values[0];
    }

    for(int i=0; i<counts_table->capacity; i++) {
        if(counts_table->entries[i] == NULL) {
            continue;
        }
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

Entry *kmer_max_entry(kmerHashTable *hash_table) {
    double max_value = -DBL_MAX;
    Entry *max_entry = NULL;
    for(int i=0; i<hash_table->capacity; i++) {
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


void entry_to_file(FILE *file, Entry *entry, char delimiter) {
    fprintf(file, "%s%c%f\n", entry->key, delimiter, entry->values[0]);
}


void free_options(options *opt) {
    if(opt->input_file && !opt->independent_probs) {
        free(opt->input_file);
    }
    if(opt->bound_file) {
        free(opt->bound_file);
    }
    if(opt->out_given) {
        free(opt->out_filename);
    }
    for(size_t i=0; i<opt->iterations; i++) {
        free(opt->top_kmer[i]);
    }
    free(opt->top_kmer);
    fclose(opt->out_file);
}


void print_options(options *opt) {
    printf("input_file: '%s'\n",opt->input_file);
    printf("bound_file: '%s'\n",opt->bound_file);
    printf("output_file: '%s'\n",opt->out_filename);
    printf("kmer: '%d'\n",opt->kmer);
    printf("iterations: '%d'\n",opt->iterations);
    printf("jobs: '%d'\n",opt->jobs);
    printf("probs: '%d'\n",opt->independent_probs);
}
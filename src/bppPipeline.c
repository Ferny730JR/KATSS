#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "bppHashTable.h"
#include "utils.h"

#include "bppPipeline_cmdl.h"

struct options {
    char    *input_file;
    char    *bound_file;
    char    *out_file;
    int     kmer;
    int     noBin;
    int     frq;
};

/*##########################################################
#  Function Declarations                                   #
##########################################################*/

void print_options(struct options *opt);


void free_options(struct options *opt);


bppHashTable bppCountKmers(char *filename, 
                           int kmer);


void process_line(char          *line, 
                  bppHashTable  counts_table, 
                  int           kmer);


void getFrequencies(bppHashTable counts_table);


bppHashTable getBPPEnrichment(bppHashTable  *control_frq, 
                              bppHashTable  *bound_frq, 
                              int           kmer);


void init_default_options(struct options *opt) {
    opt->input_file = NULL;
    opt->bound_file = NULL;
    opt->out_file   = NULL;
    opt->kmer       = 3;
    opt->noBin      = 0;
    opt->frq        = 0;
}


int main(int argc, char **argv) {
    
    // Declare variables
    struct bppPipeline_args_info    args_info;
    struct options                  opt;
    bppHashTable                    bounds_table;
    bppHashTable                    control_table;
    bppHashTable                    enrichments_table;

    init_default_options(&opt);

    /*##########################################################
    #  Parse Command Line Arguments                            #
    ##########################################################*/

    if (bppPipeline_cmdline_parser(argc, argv, &args_info) != 0)
        exit(EXIT_FAILURE);
    
    if(args_info.input_given) {
        opt.input_file = strdup(args_info.input_arg);
        if(access(opt.input_file, F_OK|R_OK) != 0) {
            error_message("Unable to open input file '%s' for reading.",args_info.bound_arg);
            exit(EXIT_FAILURE);
        }
    }

    if(args_info.bound_given) {
        opt.bound_file = strdup(args_info.bound_arg);
        if(access(opt.bound_file, F_OK|R_OK) != 0) {
            error_message("Unable to open bound file '%s' for reading.",args_info.bound_arg);
            exit(EXIT_FAILURE);
        }
    }

    if(args_info.output_given)
        opt.out_file = strdup(args_info.output_arg);
    
    if(args_info.kmer_given) {
        if(args_info.kmer_arg <= 0) {
            return 0;
        }
        opt.kmer = args_info.kmer_arg;
    }

    if(args_info.frq_given)
        opt.frq = 1;
    
    if(args_info.noBin_given)
        opt.noBin = 1;
    
    bppPipeline_cmdline_parser_free(&args_info);
    
    print_options(&opt);

    /*##########################################################
    #  Computations                                            #
    ##########################################################*/

    control_table = bppCountKmers(opt.input_file, opt.kmer);
    bounds_table  = bppCountKmers(opt.bound_file, opt.kmer);

    enrichments_table = getBPPEnrichment(&control_table, &bounds_table, opt.kmer);

    
    printBPPHashTable(&enrichments_table, opt.kmer);

    /* Clean up */
    free_hash_table(&control_table);
    free_hash_table(&bounds_table);
    free_hash_table(&enrichments_table);
    free_options(&opt);

    return 0;
}


bppHashTable bppCountKmers(char *filename, int kmer) {
    bppHashTable    counts_table;
    FILE            *read_file;

    counts_table = init_hash_table(kmer);
    read_file = fopen(filename, "r");

    char buffer[10000];
    while (fgets(buffer, sizeof(buffer), read_file)) {

        process_line(buffer, counts_table, kmer);
    }
    fclose(read_file);
    
    getFrequencies(counts_table);

    return counts_table;
}


void process_line(char *line, bppHashTable counts_table, int kmer) {
    char    *data;
    char    *sequence;
    char    k_substr[kmer+1];
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
        bpp = atof(data);
        token_bpp_values[token_count++]=bpp;
        data = strtok(NULL, " ");
    }

    for(int i=0; i<num_kmers_in_seq; i++) {
        // Get kmer substring
        memcpy( k_substr, &sequence[i], kmer );
  
        k_substr[kmer] = '\0';
        addValue(&counts_table, k_substr, 1, kmer);
        
        // Loop through bpp values in file
        for(int j=i; j<kmer+i; j++) {
            addValue(&counts_table, k_substr, token_bpp_values[j], j-i);
        }
    }
    free(sequence);
}


void getFrequencies(bppHashTable counts_table) {

    int num_columns = strlen(counts_table.entries[0].key);
    int total_count;

    for(int i=0; i<counts_table.size; i++) {
        if(counts_table.entries[i].data.values[num_columns]<1)
            continue;
        for(int j=0; j<num_columns; j++) {
            total_count=counts_table.entries[i].data.values[num_columns];
            counts_table.entries[i].data.values[j]/=total_count;
        }
    }
}


bppHashTable getBPPEnrichment(bppHashTable *control_frq, 
                              bppHashTable *bound_frq, 
                              int kmer) {

    bppHashTable    enrichments_table;
    char            *key; 
    double          enrichment;
    double          *bound_values;
    double          *control_values;
    double          *enrichment_values;
    double          mean_enrichment;

    enrichments_table = init_hash_table(kmer);

    // Get the log2 fold change for each kmer
    for(size_t i = 0; i < bound_frq->size; i++) {
        key = bound_frq->keys[i];

        bound_values   =    get(bound_frq, key);
        control_values =    get(control_frq, key);

        for(int j = 0; j < kmer; j++) {
            if(bound_values[j] == 0 || control_values[j] == 0) {
                enrichment = 0;
            } else {
                enrichment = log(bound_values[j]/control_values[j])/log(2);
            }
            addValue(&enrichments_table, key, enrichment, j);
        }
    }

    for(size_t i = 0; i < enrichments_table.size; i++) {
        key = enrichments_table.keys[i];
        enrichment_values = get(&enrichments_table, key);

        mean_enrichment = 0;
        for(int j = 0; j < kmer; j++) {
            mean_enrichment += enrichment_values[j];
        }

        mean_enrichment/=kmer;
        addValue(&enrichments_table, key, mean_enrichment, kmer);
    }
    return enrichments_table;
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

void print_options(struct options *opt) {
    printf("Input file: \"%s\"\n",opt->input_file);
    printf("Bound file: \"%s\"\n",opt->bound_file);
    printf("Output file: \"%s\"\n",opt->out_file);
    printf("Kmer: \"%d\"\n",opt->kmer);
    printf("No Bins: \"%d\"\n",opt->noBin);
    printf("Include frq: \"%d\"\n",opt->frq);
}

void free_options(struct options *opt) {
    free(opt->input_file);
    free(opt->bound_file);
    free(opt->out_file);
}
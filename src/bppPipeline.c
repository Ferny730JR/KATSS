#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

bppHashTable bppCountKmers(char *filename, 
                           int kmer);

void getFrequencies(bppHashTable counts_table);

void init_default_options(struct options *opt) {
    opt->input_file = NULL;
    opt->bound_file = NULL;
    opt->out_file   = NULL;
    opt->kmer       = 3;
    opt->noBin      = 0;
    opt->frq        = 0;
}

void print_options(struct options *opt) {
    printf("Input file: \"%s\"\n",opt->input_file);
    printf("Bound file: \"%s\"\n",opt->bound_file);
    printf("Output file: \"%s\"\n",opt->out_file);
    printf("Kmer: \"%d\"\n",opt->kmer);
    printf("No Bins: \"%d\"\n",opt->noBin);
    printf("Include frq: \"%d\"\n",opt->frq);
}


int main(int argc, char **argv) {
    

    struct bppPipeline_args_info    args_info;
    struct options                  opt;

    init_default_options(&opt);

    /*
    ############################################################
    # Options                                                  #
    ############################################################
    */

    if (bppPipeline_cmdline_parser(argc, argv, &args_info) != 0)
        exit(EXIT_FAILURE);
    
    if(args_info.input_given) {
        opt.input_file = strdup(args_info.input_arg);
        if(access(opt.input_file, F_OK|R_OK) != 0) {
            error_warning("Unable to open input file '%s' for reading.",args_info.bound_arg);
            exit(EXIT_FAILURE);
        }
    }

    if(args_info.bound_given) {
        opt.bound_file = strdup(args_info.bound_arg);
        if(access(opt.bound_file, F_OK|R_OK) != 0) {
            error_warning("Unable to open bound file '%s' for reading.",args_info.bound_arg);
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
    
    print_options(&opt);


    bppHashTable my_table = bppCountKmers(opt.input_file, opt.kmer);
    
    printBPPHashTable(&my_table, opt.kmer);

    // cleanup
    free_hash_table(&my_table);

    return 0;
}

bppHashTable bppCountKmers(char *filename, int kmer) {
    bppHashTable    counts_table;
    FILE            *read_file;
    char            *data;
    char            *sequence;
    char            k_substr[kmer+1];
    double          bpp;
    int             seq_length;
    int             num_kmers_in_seq;

    counts_table = init_hash_table(kmer);
    read_file = fopen(filename, "r");

    char buffer[10000];
    while (fgets(buffer, sizeof(buffer), read_file)) {

        char *data = strtok(buffer, " ");
        seq_length = strlen(data);
        sequence = strdup(data);

        num_kmers_in_seq = seq_length - kmer + 1;
        double token_bpp_values[num_kmers_in_seq];

        data = strtok(NULL, " ");
        int token_count = 0;
        while(data != NULL) {
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
        free(data);
    }

    fclose(read_file);
    
    getFrequencies(counts_table);

    return counts_table;
}

void getFrequencies(bppHashTable counts_table) {
    int num_columns = strlen(counts_table.entries[0].key);
    int total_count;
    for(int i=0; i<counts_table.size; i++) {
        if(counts_table.entries[i].data.values[num_columns]<1)
            continue;
        for(int j=0; j<num_columns; j++) {
            total_count=counts_table.entries[i].data.values[num_columns];
            counts_table.entries[i].data.values[j]/=counts_table.entries[i].data.values[num_columns];
        }
    }
}
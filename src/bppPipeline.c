#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bppHashTable.h"
#include "bppPipeline_cmdl.h"

struct options {
    char    *input_file;
    char    *bound_file;
    char    *out_file;
    int     kmer;
    int     noBin;
    int     frq;
};

void 
init_default_options(struct options *opt) {
    opt->input_file = NULL;
    opt->bound_file = NULL;
    opt->out_file   = NULL;
    opt->kmer       = 3;
    opt->noBin      = 0;
    opt->frq        = 0;
}

void 
print_options(struct options *opt) {
    printf("Input file: \"%s\"\n",opt->input_file);
    printf("Bound file: \"%s\"\n",opt->bound_file);
    printf("Output file: \"%s\"\n",opt->out_file);
    printf("Kmer: \"%d\"\n",opt->kmer);
    printf("No Bins: \"%d\"\n",opt->noBin);
    printf("Include frq: \"%d\"\n",opt->frq);
}

bppHashTable bppCountKmers(char *filename, int kmer);

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
        exit(1);
    
    if(args_info.input_given) {
        opt.input_file = strdup(args_info.input_arg);
    }

    if(args_info.bound_given)
        opt.bound_file = strdup(args_info.bound_arg);
    
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
    
    // do some stuff
    free_hash_table(&my_table);

    return 0;
}

bppHashTable bppCountKmers(char *filename, int kmer) {
    bppHashTable    frequency_table;
    FILE            *read_file;
    double          bpp;

    frequency_table = init_hash_table(kmer);
    read_file = fopen(filename, "r");

    char buffer[1000];
    while (fgets(buffer, sizeof(buffer), read_file)) {
        printf("==NEW LINE==\n");

        char *data = strtok(buffer, " ");
        printf("Seq: %s\n",data);

        data = strtok(NULL, " ");
        while (data !=NULL)
        {
            bpp = atof(data);   // base pair probabilities
            printf("\tFIELD %f\n", atof(data));

            /* todo - Further processing of data */

            data = strtok(NULL, " ");   // process next field
        }
        printf("\n");
    }

    fclose(read_file);

    return frequency_table;
}
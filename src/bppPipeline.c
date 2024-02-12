#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "bppHashTable.h"
#include "string_utils.h"
#include "utils.h"

#include "ViennaRNA/fold.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/utils/basic.h"

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

float* getPositionalProbabilities(char *sequence);


void print_options(struct options *opt);


void free_options(struct options *opt);


int compare(const void  *a, 
            const void  *b);


void bppHashTable_to_file(bppHashTable *table, 
                           char         *name, 
                           char         *file_extension);


void print_table_to_file(bppHashTable   *table, 
                         FILE           *table_file,
                         char           sep);


bppHashTable bppCountKmers(char *filename, 
                           int kmer);


void process_line(char *sequence, 
                  bppHashTable counts_table, 
                  int kmer);


void process_line_with_bpp(char          *line, 
                           bppHashTable  counts_table, 
                           int           kmer);


void getFrequencies(bppHashTable counts_table);


bppHashTable getBPPEnrichment(bppHashTable  *control_frq, 
                              bppHashTable  *bound_frq, 
                              int           kmer);


void init_default_options(struct options *opt) {
    opt->input_file = NULL;
    opt->bound_file = NULL;
    opt->out_file   = "rna";
    opt->kmer       = 3;
    opt->noBin      = 0;
    opt->frq        = 0;
}

/*##########################################################
#  Main                                                    #
##########################################################*/
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

    /*##########################################################
    #  Computations                                            #
    ##########################################################*/

    control_table = bppCountKmers(opt.input_file, opt.kmer);
    bounds_table  = bppCountKmers(opt.bound_file, opt.kmer);

    enrichments_table = getBPPEnrichment(&control_table, &bounds_table, opt.kmer);

    bppHashTable_to_file(&enrichments_table, opt.out_file, ".csv"); // TODO: implement proper method for changing file extension

    /* Clean up */
    free_hash_table(&control_table);
    free_hash_table(&bounds_table);
    free_hash_table(&enrichments_table);
    free_options(&opt);

    return 0;
}


float* getPositionalProbabilities(char *sequence) {
    vrna_ep_t   *ptr, *pair_probabilities = NULL;
    int         seq_length  = strlen(sequence);
    char        *propensity = (char *)vrna_alloc(sizeof(char) * (seq_length + 1));
    float       *positional_probabilities = malloc((seq_length + 1) * sizeof(float));
    float       probability;

    /* Initialize positional_probabilities array to 0's */
    memset(positional_probabilities, 0, (seq_length+1) * sizeof(float));

    /* Get the pair probabilities */
    vrna_pf_fold(sequence, propensity, &pair_probabilities);

    /* Move pair probabilities into array */
    for(ptr = pair_probabilities; ptr->i != 0; ptr++) {
        probability = ptr->p;
        positional_probabilities[ptr->i]+=probability;
        positional_probabilities[ptr->j]+=probability;
    }

    /* Clean up memory */
    free(propensity);
    free(pair_probabilities);

    return positional_probabilities;
}


bppHashTable bppCountKmers(char *filename, int kmer) {
    bppHashTable    counts_table;
    FILE            *read_file;

    counts_table = init_hash_table(kmer);
    read_file = fopen(filename, "r");

    char buffer[10000];
    while (fgets(buffer, sizeof(buffer), read_file)) {

        // Pre processing in the line
        seq_to_RNA(buffer);
        seq_to_upper(buffer);
        remove_escapes(buffer);

        process_line(buffer, counts_table, kmer);
    }
    fclose(read_file);
    
    getFrequencies(counts_table);

    return counts_table;
}


void process_line(char *sequence, bppHashTable counts_table, int kmer) {
    float  *positional_probabilities;
    //char    k_substr[kmer+1];
    char    *k_substr;
    int     seq_length = strlen(sequence);
    int     num_kmers_in_seq = seq_length - kmer + 1;

    positional_probabilities = getPositionalProbabilities(sequence);

    for(int i=0; i<num_kmers_in_seq; i++) {
        // Get kmer substring
        k_substr = substr(sequence, i, kmer);
  
        addValue(&counts_table, k_substr, 1, kmer);
        
        // Loop through bpp values in file
        for(int j=i+1; j<kmer+i+1; j++) {
            addValue(&counts_table, k_substr, positional_probabilities[j], j-i-1);
        }

        free(k_substr);
    }
    free(positional_probabilities);
}


void process_line_with_bpp(char *line, bppHashTable counts_table, int kmer) {
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
    qsort(enrichments_table.entries, enrichments_table.size, sizeof(Entry), compare);
    return enrichments_table;
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

int compare(const void *a, const void *b) {

    Entry *entryA = (Entry *)a;
    Entry *entryB = (Entry *)b;

    int mean_index = strlen(entryA->key);
    if(entryA->data.values[mean_index] > entryB->data.values[mean_index])
        return -1;
    if(entryA->data.values[mean_index] < entryB->data.values[mean_index])
        return 1;
    return 0;
}


void bppHashTable_to_file(bppHashTable *table, char *name, char *file_extension) {
    size_t name_len = strlen(name);
    size_t extension_len = strlen(file_extension);
    char *filename = malloc(name_len + extension_len + 1);

    memcpy(filename, name, name_len);
    memcpy(filename + name_len, file_extension, extension_len + 1);
    
    FILE *table_file = fopen(filename, "w");
    if (table_file == NULL) {
        error_message("Could not write to file '%s'\n",filename);
    }
    
    if(strcmp(file_extension, ".csv") == 0) {
        print_table_to_file(table, table_file, ',');
    } else if(strcmp(file_extension, ".tsv") == 0) {
        print_table_to_file(table, table_file, '\t');
    }

    free(filename);
    fclose(table_file);
}

void print_table_to_file(bppHashTable *table, FILE *table_file, char sep) {

    int num_columns = strlen(table->keys[0]) + 1;
    for(size_t i = 0; i < table->size; i++) {

        fprintf(table_file, "%s", table->entries[i].key);
        for(size_t j = 0; j < num_columns; j++) {
            fprintf(table_file, "%c%9.6f", sep, table->entries[i].data.values[j]);
        }
        fprintf(table_file,"\n");
    }
}


void free_options(struct options *opt) {
    free(opt->input_file);
    free(opt->bound_file);
}


void print_options(struct options *opt) {
    printf("Input file: \"%s\"\n",opt->input_file);
    printf("Bound file: \"%s\"\n",opt->bound_file);
    printf("Output file: \"%s\"\n",opt->out_file);
    printf("Kmer: \"%d\"\n",opt->kmer);
    printf("No Bins: \"%d\"\n",opt->noBin);
    printf("Include frq: \"%d\"\n",opt->frq);
}

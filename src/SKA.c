#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "kmerHashTable.h"
#include "string_utils.h"
#include "utils.h"

#include "SKA_cmdl.h"

struct options {
    char    *input_file;
    char    *bound_file;
    char    *out_file;
    int     out_given;
    int     kmer;
    char    file_delimiter;

    int     independent_probs;
    int     frq;
};


char delimiter_to_char(char *user_delimiter);


void free_options(struct options *opt);


void init_default_options(struct options *opt) {
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
    struct options          opt;

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

    SKA_cmdline_parser_free(&args_info);
    free_options(&opt);
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

void free_options(struct options *opt) {
    if(opt->input_file) {
        free(opt->input_file);
    }
    if(opt->bound_file) {
        free(opt->bound_file);
    }
    if(opt->out_given) {
        free(opt->out_file);
    }
}
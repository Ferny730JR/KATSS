#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rna_file_parser.h"
#include "string_utils.h"
#include "utils.h"

/* Function declarations */

void parse_fasta(RNA_FILE *rna_file, char **ret_seq);

void parse_fastq(RNA_FILE *rna_file, char **ret_seq);

void parse_reads(RNA_FILE *rna_file, char **ret_seq);

int is_nucleotide(char character);

int is_full(char *buffer, int buffer_size);


/*##########################################################
#  Main Functions (Used in header)                         #
##########################################################*/

RNA_FILE *rnaf_open(char* filename) {
    RNA_FILE *rna_file = s_malloc(sizeof *rna_file);
    rna_file->file = fopen(filename, "r");
    memset(rna_file->buffer, 0, sizeof(char)*MAX_SEQ_LENGTH); // init buffer to '\0'

    /* Check if we can open file for reading */
    if( (rna_file->file) == NULL ) {
        free(rna_file);
        error_message("Failed to open file '%s'",filename);
        return NULL;
    }

    /* Check if file contains a valid line */
    rna_file->end_of_file = fgets(rna_file->buffer, sizeof(rna_file->buffer), rna_file->file);
    if(rna_file->end_of_file == NULL) { 
        fclose(rna_file->file);
        free(rna_file);
        warning_message("File '%s' contains no sequences.",filename);
        return NULL;
    }

    return rna_file;
}


char *rnaf_get(RNA_FILE *rna_file) {
    char *seq = NULL;
    
    /* Check the type of file format */
    if(rna_file->buffer[0] == '>') {    // fasta file
        parse_fasta(rna_file, &seq);

    } else if(rna_file->buffer[0] == '@') { // fastq file
        parse_fastq(rna_file, &seq);

    } else if( is_nucleotide(rna_file->buffer[0]) ) {   // pure reads file
        parse_reads(rna_file, &seq);

    } else {
        warning_message("Unable to read sequence from file.\nCurrent supported file types are:"
        "FASTA, FASTQ, and files containing sequences per line.");
    }

    return seq;
}


void rnaf_close(RNA_FILE *rna_file) {
    fclose(rna_file->file);
    free(rna_file);
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

void parse_fasta(RNA_FILE *rna_file, char **ret_seq) {
    char    *seq;

    /* using 'seq = NULL' results in seg fault from append because of strlen. */
    seq = s_malloc(0); /* Init seq, will be realloc'd based on input stream */
    while (fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file)) { 
        if (rna_file->buffer[0] == '>') {
            break;      /* New sequence, so break */
        }

        if( !is_nucleotide(rna_file->buffer[0]) ) {
            continue;   /* Still reading from metadata, so skip iteration */
        }

        /* Append sequence found in buffer into seq variable */
        append(&seq, rna_file->buffer);
    }

    /* Have ret_seq point to new seq */
    *ret_seq  = seq;
}


void parse_fastq(RNA_FILE *rna_file, char **ret_seq) {
    char            *seq;
    unsigned int    reading_seq = 1;

    /* using 'seq = NULL' results in seg fault from append because of strlen. */
    seq = s_malloc(0); /* Init seq, will be realloc'd based on input stream */
    while(fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file)) {
        if(rna_file->buffer[0] == '@' ) {
            break;  /* New sequence, so break */
        }

        if ( !is_nucleotide(rna_file->buffer[0]) ) {
            reading_seq = 0;
            continue;   /* Reading from metadata, so skip iteration */
        }

        /* Append sequence found in buffer into seq variable */
        if(reading_seq) {
            append(&seq, rna_file->buffer);
        }
    }

    /* Have ret_seq point to new seq */
    *ret_seq = seq;
}


void parse_reads(RNA_FILE *rna_file, char **ret_seq) {
    char    *seq;

    /* using 'seq = NULL' results in seg fault from append because of strlen. */
    seq = NULL;
    append(&seq, rna_file->buffer);

    /* If buffer was not large enough to store sequence, keep reading */
    while(is_full(rna_file->buffer, MAX_SEQ_LENGTH)) {
        memset(&rna_file->buffer[0], 0, MAX_SEQ_LENGTH);
        fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file);
        append(&seq, rna_file->buffer);
    }
    
    /* If EOF has not been reached, return seq and read next line */
    if(rna_file->end_of_file) {
        *ret_seq = seq;
        rna_file->end_of_file = fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file);
    } else {
        free(seq);
    }
}


int is_nucleotide(char character) {
    int ret;

    switch(character) {
        case 'A':   ret = 1;    break;
        case 'a':   ret = 1;    break;
        case 'C':   ret = 1;    break;
        case 'c':   ret = 1;    break;
        case 'G':   ret = 1;    break;
        case 'g':   ret = 1;    break;
        case 'T':   ret = 1;    break;
        case 't':   ret = 1;    break;
        case 'U':   ret = 1;    break;
        case 'u':   ret = 1;    break;
        default:    ret = 0;    break;
    }

    return ret;
}


int is_full(char *buffer, int buffer_size) {
    int ret=0;

    if( buffer[buffer_size-1] == '\0' && is_nucleotide(buffer[buffer_size-2]) ) { // full
        ret=1;
    }

    return ret;
}
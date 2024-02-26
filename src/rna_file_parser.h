#ifndef FILE_PARSER
#define FILE_PARSER

#include <stdio.h>

#define MAX_SEQ_LENGTH 1000

typedef struct {
    FILE *file;
    char buffer[MAX_SEQ_LENGTH];
    char *end_of_file;
} RNA_FILE;


RNA_FILE *rnaf_open(char *file);


char *rnaf_get(RNA_FILE *rna_file);


void rnaf_close(RNA_FILE *rna_file);


#endif // FILE_PARSER
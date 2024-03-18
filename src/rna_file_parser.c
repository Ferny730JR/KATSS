#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rna_file_parser.h"
#include "string_utils.h"
#include "utils.h"

#define NO_OF_CHARS 256

/* Function declarations */
static void
parse_fasta(RNA_FILE *rna_file, char **ret_seq);

static void
parse_fastq(RNA_FILE *rna_file, char **ret_seq);

static void
parse_reads(RNA_FILE *rna_file, char **ret_seq);

static char
determine_filetype(char peek);

static int
is_nucleotide(char character);

static int
is_full(char *buffer, int buffer_size);

static void
badCharHeuristic(const char* str, int size, int badchar[NO_OF_CHARS]);


/*##########################################################
#  Main Functions (Used in header)                         #
##########################################################*/

RNA_FILE *
rnaf_open(char* filename) 
{
	RNA_FILE *rna_file = s_malloc(sizeof *rna_file);
	rna_file->file = fopen(filename, "r");
	rna_file->buffer = s_malloc(MAX_SEQ_LENGTH * sizeof(char));
	rna_file->buffer_size = MAX_SEQ_LENGTH;
	memset(rna_file->buffer, 0, MAX_SEQ_LENGTH * sizeof(char)); // init buffer to '\0'

	/* Check if we can open file for reading */
	if( (rna_file->file) == NULL ) {
		free(rna_file->buffer);
		free(rna_file);
		error_message("Failed to open file '%s'",filename);
		return NULL;
	}

	/* Check if file contains a valid line */
	if(!fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file)) { 
		fclose(rna_file->file);
		free(rna_file->buffer);
		free(rna_file);
		warning_message("File '%s' contains no sequences.",filename);
		return NULL;
	}

	/* Determine what type of file was passed */
	rna_file->filetype = determine_filetype(rna_file->buffer[0]);

	/* Reset the file to read from beginning */
	if(rna_file->filetype == 'r') {
		rewind(rna_file->file);
	}

	return rna_file;
}


char *
rnaf_get(RNA_FILE *rna_file) 
{
	char *seq = NULL;

	/* Check the type of file format */
	switch(rna_file->filetype) {
		case 'a':   parse_fasta(rna_file, &seq);    break;
		case 'q':   parse_fastq(rna_file, &seq);    break;
		case 'r':   parse_reads(rna_file, &seq);    break;
		default:
			error_message("Unable to read sequence from file.\nCurrent supported file types are:"
			" FASTA, FASTQ, and files containing sequences per line.");
			break;
	}

	return seq;
}


size_t 
rnaf_oread(RNA_FILE *rna_file, unsigned int offset)
{
	for(unsigned int i=0; i<offset; i++) {
		unsigned int indx = rna_file->buffer_size-(offset-i);
		rna_file->buffer[i] = rna_file->buffer[indx];
	}

	size_t len, ret;
	len = rna_file->buffer_size - offset;
	ret = fread(rna_file->buffer+offset, 1, len, rna_file->file);

	if(len == ret) {
		return ret;
	}

	/* EOF has been reached, fill unfilled buffer with NULL */
	memset(rna_file->buffer+offset+ret, 0, (rna_file->buffer_size-(offset+ret))*sizeof(char));
	return ret;
}


void 
rnaf_close(RNA_FILE *rna_file) 
{
	fclose(rna_file->file);
	free(rna_file->buffer);
	free(rna_file);
}


void
rnaf_rebuff(RNA_FILE *rna_file, unsigned int size)
{
	rna_file->buffer = s_realloc(rna_file->buffer, size+1);
	rna_file->buffer_size = size;
	rewind(rna_file->file);
	memset(rna_file->buffer, 0, (size+1) * sizeof(char));
}


unsigned int
rnaf_search(RNA_FILE *rna_file, const char *sequence)
{
	unsigned int counts = 0;
	unsigned int m = strlen(sequence);
	unsigned int n = rna_file->buffer_size;

	if(n==0) {
		return 0;
	}

	int badchar[NO_OF_CHARS];
	badCharHeuristic(sequence, m, badchar);

	unsigned int s = 0;
	while (s <= (n - m)) {
		int j = m - 1;

		while (j >= 0 && sequence[j] == rna_file->buffer[s + j]) {
			j--;
		}

		if (j < 0) {
			// printf("pattern occurs at shift = %d\n", s);
			counts++;
			s += (s + m < n) ? m - badchar[(uint8_t)rna_file->buffer[s + m]] : 1;
		}

		else {
			s += MAX2(1, j - badchar[(uint8_t)rna_file->buffer[s + j]]);
		}
	}

	return counts;
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

static void
parse_fasta(RNA_FILE *rna_file, char **ret_seq) 
{   /* Init seq, will be realloc'd based on input stream */
	char    *seq = NULL;

	while (fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file)) { 
		if (rna_file->buffer[0] == '>') {
			break;      /* New sequence, so break */
		}

		/* Append sequence found in buffer into seq variable */
		append(&seq, rna_file->buffer);
	}

	/* Have ret_seq point to new seq */
	*ret_seq = seq;
}


static void
parse_fastq(RNA_FILE *rna_file, char **ret_seq) 
{
	char            *seq = NULL; /* Init seq, will be realloc'd based on input stream */
	unsigned int    reading_seq = 1;

	while(fgets(rna_file->buffer, MAX_SEQ_LENGTH, rna_file->file)) {
		if(rna_file->buffer[0] == '@' ) {
			break;  /* New sequence, so break */
		}

		if ( rna_file->buffer[0] == '+' ) {
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


static void
parse_reads(RNA_FILE *rna_file, char **ret_seq) 
{
	char *seq = NULL;

	/* Get next line, and if NULL, return */
	if(!fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file)) {
		return;
	}

	/* Append the read file to seq */
	append(&seq, rna_file->buffer);

	/* If buffer was not large enough to store sequence, keep reading */
	while(is_full(rna_file->buffer, rna_file->buffer_size)) {
		memset(rna_file->buffer, 0, rna_file->buffer_size*sizeof(char));
		fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file);
		append(&seq, rna_file->buffer);
	}

	/* Have ret_seq point to new seq */
	*ret_seq = seq;
}


static char
determine_filetype(char peek) 
{
	char ret;

	if( is_nucleotide(peek) ) {
		ret = 'r';  /* reads file */
	} else if(peek == '@') {
		ret = 'q';  /* fastq file */
	} else if(peek == '>') {
		ret = 'a';  /* fasta file */
	} else {
		ret = '\0'; /* unsupported file type */
	}

	return ret;
}


static int
is_nucleotide(char character) 
{
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


int 
is_full(char *buffer, int buffer_size) 
{
	int ret=0;

	if( buffer[buffer_size-1] == '\0' && is_nucleotide(buffer[buffer_size-2]) ) { // full
		ret=1;
	}

	return ret;
}


// The preprocessing function for Boyer Moore's bad character heuristic
void 
badCharHeuristic(const char* str, int size, int badchar[NO_OF_CHARS])
{
	int i;

	for (i = 0; i < NO_OF_CHARS; i++) {
		badchar[i] = -1;
	}

	for (i = 0; i < size; i++) {
		badchar[(int)str[i]] = i;
	}
}

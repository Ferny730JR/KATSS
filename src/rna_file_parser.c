#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

#include "rna_file_parser.h"
#include "string_utils.h"
#include "utils.h"

#define NO_OF_CHARS 16

/* Struct to contain getm information */
typedef struct getm_line_info {
	int shift;
	char *line;
} getm_line_info;

/* Function declarations */
static void
parse_fasta(RNA_FILE *rna_file, char **ret_seq);

static void
parse_fastq(RNA_FILE *rna_file, char **ret_seq);

static void
parse_reads(RNA_FILE *rna_file, char **ret_seq);

static getm_line_info
getm_line(char *search, unsigned int search_len, unsigned int found_at);

static char
determine_filetype(char peek);

static bool
is_nucleotide(char character);

static bool
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
	rna_file->filename = filename;
	rna_file->buffer = s_malloc(MAX_SEQ_LENGTH * sizeof(char));
	rna_file->buffer_size = MAX_SEQ_LENGTH;
	rna_file->num_chars = 0;
	rna_file->num_lines = 0;
	rna_file->shift = 0;
	memset(rna_file->buffer, 0, MAX_SEQ_LENGTH * sizeof(char)); // init buffer to '\0'

	/* Check if we can open file for reading */
	if( (rna_file->file) == NULL ) {
		free(rna_file->buffer);
		free(rna_file);
		error_message("Failed to open file '%s': %s",filename, strerror(errno));
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


char *
rnaf_getm(RNA_FILE *rna_file, char *match)
{
	getm_line_info ret;
	size_t         still_reading = 0;
	char          *line;
	int            m = strlen(match);
	int            n = rna_file->buffer_size;

	if(n==0) {
		return NULL;
	}

	int badchar[NO_OF_CHARS];
	badCharHeuristic(match, m, badchar);
	do {
	
	while (rna_file->shift <= (n - m)) {
		int j = m - 1;

		while (j >= 0 && match[j] == rna_file->buffer[rna_file->shift + j]) {
			j--;
		}

		if (j < 0) {
			ret = getm_line(rna_file->buffer, rna_file->buffer_size, rna_file->shift);
			if(ret.line) {
				rna_file->shift += (rna_file->shift + m < n) ? ret.shift : 1;
				return ret.line;
			} else {
				rnaf_oread(rna_file, n - rna_file->shift + 1);
				rna_file->shift = 0;
			}
		} else {
			rna_file->shift+=MAX2(1,j-badchar[(unsigned int)rna_file->buffer[rna_file->shift+j]]);
		}
	}

	int offset = 0;
	while(rna_file->buffer[rna_file->buffer_size-offset-1] != '\n') {
		offset++;
	}

	rna_file->shift = 0;
	still_reading = rnaf_oread(rna_file, offset);

	} while(still_reading);

	printf("NOT FOUND ABORT\n");

	return NULL;
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
	if(ret > 0) {
		memset(rna_file->buffer+offset+ret, 0, (rna_file->buffer_size-(offset+ret))*sizeof(char));
	} else { /* Nothing left to read, empty buffer */
		memset(rna_file->buffer, 1, rna_file->buffer_size*sizeof(char));
	}
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
			counts++;
			s += (s + m < n) ? m - badchar[(unsigned int)rna_file->buffer[s + m]] : 1;
		}

		else {
			s += MAX2(1, j - badchar[(unsigned int)rna_file->buffer[s + j]]);
		}
	}

	return counts;
}


unsigned long
rnaf_numchars(RNA_FILE *rna_file)
{	/* If numchars has already been calculated, return */
	if(rna_file->num_chars) {
		return rna_file->num_chars;
	}

	/* Open a new FILE pointer to not mess up RNA_FILE's pointer */
	FILE *fp = fopen(rna_file->filename, "r");
	if(fp == NULL) {
		error_message("rnaf Failed get numchars: %s", strerror(errno));
		return 0;
	}

	/* Count the number of characters read in file */
	char buf[65536];
	unsigned long total_chars = 0;
	while(1) {
		size_t res = fread(buf, 1, 65536, fp);
		if(ferror(fp)) {
			error_message("Failed to read file: %s", strerror(errno));
			return 0;
		}

		total_chars += res;

		if(feof(fp)) {
			break;
		}
	}

	/* Update num_chars and return */
	rna_file->num_chars = total_chars;
	return total_chars;
}


unsigned long
rnaf_numlines(RNA_FILE *rna_file)
{
	if(rna_file->num_lines) {
		return rna_file->num_lines;
	}

	FILE *fp = fopen(rna_file->filename, "r");
	if(fp == NULL) {
		error_message("rnaf Failed get numlines: %s", strerror(errno));
		return 0;
	}

	char buf[65536];
	unsigned long counter = 0;
	while(1) {
		size_t res = fread(buf, 1, 65536, fp);
		if(ferror(fp)) {
			error_message("Failed to read file: %s", strerror(errno));
			return 0;
		}

		for(size_t i = 0; i < res; i++) {
			if(buf[i] == '\n') {
				counter++;
			}
		}

		if(feof(fp)) {
			break;
		}
	}

	rna_file->num_lines = counter;
	return counter;
}


/*##########################################################
#  Helper Functions                                        #
##########################################################*/

static void
parse_fasta(RNA_FILE *rna_file, char **ret_seq) 
{   /* Init seq, will be realloc'd based on input stream */
	char    *seq = NULL;

	while (fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file)) { 
		if (rna_file->buffer[0] == '>') {
			/* If buffer wasn't large enough to store header, keep reading */
			while(is_full(rna_file->buffer, rna_file->buffer_size)) {
				fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file);
			}
			break;  /* New sequence, so break */
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
	unsigned int	current_iter = 0;

	while(fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file)) {
		current_iter+=1;

		if(rna_file->buffer[0] == '@' && current_iter == 4) {
			/* If buffer wasn't large enough to store header, keep reading */
			while(is_full(rna_file->buffer, rna_file->buffer_size)) {
				fgets(rna_file->buffer, rna_file->buffer_size, rna_file->file);
			}
			break;  /* New sequence, so break */
		}

		if(rna_file->buffer[0] == '+') {
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


static getm_line_info
getm_line(char *search, unsigned int search_len, unsigned int found_at) {
	unsigned int line_start = found_at;
	unsigned int line_length = 1;

	getm_line_info info = {.shift = 0, .line=NULL};

	while(search[line_start-1] != '\n') {
		line_start--;
		if(line_start == 0) {
			break;
		}
	}

	while(search[line_start+line_length] != '\n') {
		if(search[line_start+line_length] == '\0' || line_start+line_length>65536) {
			return info;
		}
		line_length++;
	}

	info.shift = line_length - (found_at - line_start);
	info.line = substr(search, line_start, line_length);

	return info;
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


static bool
is_nucleotide(char character)
{
	switch(character) {
		case 'A':   return true;
		case 'a':   return true;
		case 'C':   return true;
		case 'c':   return true;
		case 'G':   return true;
		case 'g':   return true;
		case 'T':   return true;
		case 't':   return true;
		case 'U':   return true;
		case 'u':   return true;
		default:    return false;
	}
}


static bool
is_full(char *buffer, int buffer_size)
{
	return buffer[buffer_size-1] == '\0' && is_nucleotide(buffer[buffer_size-2]);
}


// The preprocessing function for Boyer Moore's bad character heuristic
static void 
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

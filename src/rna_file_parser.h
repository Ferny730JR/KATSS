#ifndef RNA_FILE_PARSER
#define RNA_FILE_PARSER

#include <stdio.h>

#define MAX_SEQ_LENGTH 1000

/**
 *  Represents an RNA file for reading, including file information and a character buffer.
 *  The struct is used in conjunction with RNA file parsing functions.
 */
typedef struct RNA_FILE {
    FILE *file;                     /** File pointer for the RNA file. */
    char *buffer;                   /** Character buffer to store sequence data. */
	unsigned int buffer_size;       /** Int to store the size of the buffer */
    char *end_of_file;              /** Pointer to the end of the file indicator. */
    char filetype;                  /** Character to store which file type was passed. */
} RNA_FILE;


/**
 *  @brief Opens an RNA file for reading.
 *
 *  This function opens an RNA file specified by the provided filename and returns a pointer
 *  to the RNA_FILE struct representing the opened file. This library currently only supports
 *  file reading for FASTA, FASTQ, text files with sequences.
 *
 *  @param  filename The name of the RNA file to be opened.
 *  @return A pointer to the RNA_FILE struct representing the opened file, or NULL if there was an
 *  error.
 */
RNA_FILE *rnaf_open(char *filename);


/**
 *  @brief Retrieves the next sequence from the RNA file.
 *
 *  This function reads the next sequence from the RNA file represented by the provided RNA_FILE
 *  struct. It automatically detects the file format (FASTA, FASTQ, or sequences per line) and
 *  parses accordingly. The retrieved sequence is dynamically allocated and returned.
 * 
 *  @note You have to free the string. Since memory is allocated to store the string, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 *
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file.
 *  @return A dynamically allocated string containing the sequence, or NULL if there are no more
 *  sequences or an error occurs.
 */
char *rnaf_get(RNA_FILE *rna_file);


/**
 *  @brief Closes the RNA file.
 *
 *  This function closes the RNA file represented by the provided RNA_FILE struct and frees
 *  associated resources.
 *
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file.
 */
void rnaf_close(RNA_FILE *rna_file);


/**
 *  @brief Change the size of the buffer
 * 
 *  This function reallocates the amount of memory the buffer uses.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file
 *  @param size The new size of the buffer in bytes (number of chars)
*/
void
rnaf_rebuff(RNA_FILE *rna_file, unsigned int size);


/**
 *  @brief Search for sequence in RNA_FILE
*/
unsigned int
rnaf_search(RNA_FILE *rna_file, const char *sequence);


void rnaf_read(RNA_FILE *rna_file, unsigned int offset);

#endif // FILE_PARSER

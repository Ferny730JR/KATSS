#ifndef RNA_FILE_PARSER
#define RNA_FILE_PARSER

#include <zlib.h>
#include <stdbool.h>

#define MAX_SEQ_LENGTH 1000

/**
 *  Represents an RNA file for reading, including file information and a character buffer.
 *  The struct is used in conjunction with RNA file parsing functions.
 */
typedef struct RNA_FILE {
	gzFile file;                    /** File pointer for the RNA file. */
	char *filename;                 /** String storing the name of the opened file. */
    char *buffer;                   /** Character buffer to store sequence data. */
	char *getm_ptr;                 /** Pointer that keeps track of the last position in buffer */
	unsigned int buffer_size;       /** Int to store the size of the buffer.  */
	unsigned long num_chars;        /** Number to store the total number of chars in file. */
	unsigned long num_lines;        /** Number to store the total number of lines in file. */
    char filetype;                  /** Character to store which file type was passed. */
	bool is_t;                      /** Determine if file is T or U based */
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
 *  @brief Retrieves the next sequence that contains the string `match` within it.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file.
 *  @param match The string to search for.
 * 
 *  @return A string of the matched sequence, or NULL if there are no more
 *  sequences or an error occurs.
 * 
 *  @note Recommend to set the rna_file buffer size to 65536
*/
char *rnaf_getm(RNA_FILE *rna_file, char *match);


/**
 *  @brief Read from the file until buffer is full.
 *
 *  Fills the first offset number of characters of the buffer with the last offset
 *  number of characters of the previous buffer. After which, the rest of the buffer is filled
 *  with the new reads.
 *
 *  For example, Let's assume you had a file with the line: `"Hello, world! How are you?"`, and the
 *  current buffer of size 13 was filled with `"Hello, world!"`. Calling `rnaf_oread(rna_file, 6)`
 *  would set the new buffer as: `"world! How ar"`.
 *
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file
 *  @param offset The number of characters to read from the previous buffer
 *
 *  @return Total number of elements successfully read
*/
size_t rnaf_oread(RNA_FILE *rna_file, unsigned int offset);


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
 *  @brief Change the size of the buffer.
 * 
 *  This function reallocates the amount of memory the buffer uses.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file
 *  @param size The new size of the buffer in bytes (number of chars)
*/
void
rnaf_rebuff(RNA_FILE *rna_file, unsigned int size);


/**
 *  @brief Search for sequence in RNA_FILE.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file
 *  @param sequence The sequence to search for in RNA_FILE
 * 
 *  @return The total number of times sequence was found in rna_file
*/
unsigned int
rnaf_search(RNA_FILE *rna_file, const char *sequence);


/** 
 *  @brief Count the number of characters in RNA_FILE.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file
 * 
 *  @return The total number of characters in rna_file
 *  
*/
unsigned long
rnaf_numchars(RNA_FILE *rna_file);


/** 
 *  @brief Count the number of lines in RNA_FILE.
 *
 *  This function counts the total number of newline characters present in RNA_FILE.
 * 
 *  @param rna_file A pointer to the RNA_FILE struct representing the opened file.
 * 
 *  @return The number of lines in RNA_FILE
*/
unsigned long
rnaf_numlines(RNA_FILE *rna_file);

#endif // RNA_FILE_PARSER

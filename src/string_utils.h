#ifndef STRING_UTILS_H
#define STRING_UTILS_H


char* substr(char *sequence, int start, int length);


void seq_to_upper(char *sequence);


void seq_to_RNA(char *sequence);


void remove_escapes(char *sequence);


#endif
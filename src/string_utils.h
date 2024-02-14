#ifndef STRING_UTILS_H
#define STRING_UTILS_H


char* substr(char *sequence, int start, int length);


char* prefix_of_str(char* str);


char* concat(const char *s1, const char *s2);


void seq_to_upper(char *sequence);


void seq_to_RNA(char *sequence);


void remove_escapes(char *sequence);


#endif
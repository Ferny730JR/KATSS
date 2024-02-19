#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "utils.h"
#include "string_utils.h"

char* substr(char *sequence, int start, int length) {
    char*   substring = (char*) s_malloc((length+1)*sizeof(char));
    int     seq_length = strlen(sequence);
    int     substr_max_length = seq_length - start;

    if(start < 0 || start > seq_length) {
        warning_message("'start' value of %d not valid for sequence '%s' not valid.",start,sequence);
        free(substring);
        return NULL;
    }

    if(!sequence) {
        error_message("Unable to read string %s",sequence);
    }
    
    if(length > substr_max_length) {
        memcpy(substring, &sequence[start], substr_max_length);
        substring[substr_max_length] = '\0';
    } else {
        memcpy(substring, &sequence[start], length);
        substring[length] = '\0';
    }

    return substring;
}


char* prefix_of_str(char* str) {
    char *prefix;

    prefix = strtok(str, ".");

    return strdup(prefix);
}


char* concat(const char *s1, const char *s2) {
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = malloc(len1 + len2 + 1); // +1 for the null-terminator
    if(!result) {
        error_message("Failed to allocate memory for concatenation of '%s' and '%s'.",s1,s2);
    }

    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1); // +1 to copy the null-terminator
    return result;
}


void str_to_upper(char *str) {
    if(!str) {
        error_message("Unable to read string %s",str);
    }

    for(int i=0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}


void seq_to_RNA(char *sequence) {
    unsigned int i;

    if(!sequence) {
        error_message("Unable to read string %s",sequence);
    }

    for (i = 0; sequence[i]; i++) {
        if (sequence[i] == 'T') {
            sequence[i] = 'U';
        } else if (sequence[i] == 't') {
            sequence[i] = 'u';
        }
    }
}


void remove_escapes(char *sequence) {

    int ln = strlen(sequence)-1;

    if (sequence[ln] == '\n') {  // remove trailing new line character
        sequence[ln] = '\0';
    }
    
    while(isspace(*sequence)) {  // move pointer past white space
        ++sequence;
    }
}
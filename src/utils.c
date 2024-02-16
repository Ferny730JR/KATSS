#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "utils.h"

void _error_message(const char *format, va_list args);


void _warning_message(const char *format, va_list args);


void error_message(const char *format, ...) {
    va_list args;

    va_start(args, format);
    _error_message(format, args);
    va_end(args);
}


void _error_message(const char *format, va_list args) {
    fprintf(stderr, ANSI_COLOR_RED "ERROR: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
}


void warning_message(const char *format, ...) {
    va_list args;

    va_start(args, format);
    _warning_message(format, args);
    va_end(args);
}


void _warning_message(const char *format, va_list args) {
    fprintf(stderr, ANSI_COLOR_RED "WARNING: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
}


void *s_malloc(size_t mem_size) {
    void *pointer = malloc(mem_size);

    if(!pointer && mem_size) {
        error_message("Could not allocate memory.");
        exit(EXIT_FAILURE);
    }

    return pointer;
}


void *s_calloc(size_t count, size_t mem_size) {
    void *pointer = calloc(count, mem_size);

    if(!pointer && mem_size && count) {
        error_message("Could not allocate memory.");
        exit(EXIT_FAILURE);
    }

    return pointer;
}
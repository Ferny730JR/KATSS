#include <stdio.h>
#include <stdarg.h>

#include "utils.h"

void _error_warning(const char *format, va_list args);

void error_warning(const char *format, ...) {
    va_list args;

    va_start(args, format);
    _error_warning(format, args);
    va_end(args);
}

void _error_warning(const char *format, va_list args) {
    fprintf(stderr, ANSI_COLOR_RED "ERROR: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
}
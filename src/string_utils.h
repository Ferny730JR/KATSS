#ifndef STRING_UTILS_H
#define STRING_UTILS_H


/**
 *  @brief Returns a substring of sequence that is length-characters long starting at 
 *  character start.
 *  
 *  The start parameter defines the starting index of the sequence string
 * 
 *  @attention You have to free the substring. Since memory is allocated to store the substring, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 * 
 *  @param sequence The sequence to get the substring from
 *  @param start    The starting index of sequence for substring
 *  @param length   The number of characters in the substring
 * 
 *  @return Pointer to the substring
*/
char *substr(const char *sequence, const int start, const int length);


/**
 *  @brief Get the basename prefix of a file path.
 * 
 *  This function removes all characters before the last occurrence of the character `/`, and all
 *  characters following the first `'.'`. Assuming `full_path` were to be: `/usr/bin/file.txt`, 
 *  this function would return the string `file`. If the string were to contain no `/` or `.`
 *  characters, then it returns a duplicate of the string.
 * 
 *  @note You have to free the returned string. Since memory is allocated to store the string, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param full_path    The character pointer containing the file path.
 * 
 *  @return char pointer to the basename prefix
*/
char *basename_prefix(const char *full_path);


/**
 *  @brief Concatenate two strings.
 * 
 *  The original strings passed in the argument will remain unaffected, since a new string
 *  containing the combined contents will be returned. 
 * 
 *  @note You have to free the returned string. Since memory is allocated to store the new string,
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param s1   String to concatenate to
 *  @param s2   Appends string to `s1`
 * 
 *  @return char pointer to concatenated string.
*/
char *concat(const char *s1, const char *s2);


/**
 *  @brief Appends the contents of the second string to the end of the first string.
 * 
 *  The function allocates memory for the combined result and updates the first string accordingly.
 *  If the first string is null or empty, it essentially duplicates the contents of `s2` into `s1`.
 *  Similarly, if `s2` is null or empty, then the contents of `s1` will remain the same.
 *
 *  @param s1 The pointer to the first string (modifiable).
 *  @param s2 The second string to append.
 */
void append(char **s1, const char *s2);


/**
 *  @brief Finds the starting index of the first occurrence of a substring in a given string.
 *
 *  This function searches for the first occurrence of the substring specified by `s2` within
 *  the string `s1` and returns the starting index of that occurrence. If the substring is not
 *  found, the function returns -1.
 *
 *  @param s1 The null-terminated string in which the substring is searched.
 *  @param s2 The null-terminated substring to be located within the string `s1`.
 *  @return The starting index of the first occurrence of the substring, or -1 if not found.
 */
int subindx(const char *s1, const char *s2);


/**
 *  @brief Capitalizes every lower case letter in string.
 * 
 *  @param str  String to capitalize
*/
void str_to_upper(char *str);


/**
 *  @brief Converts from DNA alphabet to RNA.
 * 
 *  The function will substitute all 'T' and 't' characters with 'U' and 'u' characters,
 *  respectively.
 * 
 *  @param sequence The sequence to substitute characters
*/
void seq_to_RNA(char *sequence);


void remove_escapes(char *sequence);


#endif
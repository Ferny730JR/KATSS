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
char *substr(char *sequence, int start, int length);


/**
 *  @brief Return basename of string containing filepath.
*/
char *prefix_of_str(char* str);


/**
 *  @brief Get the basename prefix of a file path.
 * 
 *  This function removes all characters before the last occurrence of the character '/', and all
 *  characters following the first '.'. Assuming the full_path were to be: '/usr/bin/file.txt', 
 *  this function would return the string 'file'. If the string were to contain no '/' or '.'
 *  characters, then it returns a copy of the string as is.
 * 
 *  @attention You have to free the string. Since memory is allocated to store the string, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param full_path    The character pointer containing the file path.
 * 
 *  @return Pointer to the basename prefix
*/
char *basename_prefix(char *full_path);


/**
 *  @brief Concatenate two strings.
 *  @attention You have to free the new string. Since memory is allocated to store the new string, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param s1   String to concatenate to
 *  @param s2   Appends string to s1
 * 
 *  @return Pointer to concatenated string.
*/
char *concat(const char *s1, const char *s2);


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
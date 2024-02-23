/** @file SKA_cmdl.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef SKA_CMDL_H
#define SKA_CMDL_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef SKA_CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define SKA_CMDLINE_PARSER_PACKAGE "Streaming K-mer Analysis"
#endif

#ifndef SKA_CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define SKA_CMDLINE_PARSER_PACKAGE_NAME "Streaming K-mer Analysis"
#endif

#ifndef SKA_CMDLINE_PARSER_VERSION
/** @brief the program version */
#define SKA_CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct SKA_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *detailed_help_help; /**< @brief Print help, including all details and hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * input_arg;	/**< @brief Read the controls sequences file.
.  */
  char * input_orig;	/**< @brief Read the controls sequences file.
 original value given at command line.  */
  const char *input_help; /**< @brief Read the controls sequences file.
 help description.  */
  char * bound_arg;	/**< @brief Read the protein bound RNA sequences file.
.  */
  char * bound_orig;	/**< @brief Read the protein bound RNA sequences file.
 original value given at command line.  */
  const char *bound_help; /**< @brief Read the protein bound RNA sequences file.
 help description.  */
  char * output_arg;	/**< @brief Set the name of the output files.
.  */
  char * output_orig;	/**< @brief Set the name of the output files.
 original value given at command line.  */
  const char *output_help; /**< @brief Set the name of the output files.
 help description.  */
  int kmer_arg;	/**< @brief Set the length of k-mers.
 (default='5').  */
  char * kmer_orig;	/**< @brief Set the length of k-mers.
 original value given at command line.  */
  const char *kmer_help; /**< @brief Set the length of k-mers.
 help description.  */
  char * file_delimiter_arg;	/**< @brief Set the delimiter used to separate the values in the output file.
 (default=',').  */
  char * file_delimiter_orig;	/**< @brief Set the delimiter used to separate the values in the output file.
 original value given at command line.  */
  const char *file_delimiter_help; /**< @brief Set the delimiter used to separate the values in the output file.
 help description.  */
  int independent_probs_flag;	/**< @brief Calculate the enrichments without the input reads.
 (default=off).  */
  const char *independent_probs_help; /**< @brief Calculate the enrichments without the input reads.
 help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int detailed_help_given ;	/**< @brief Whether detailed-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int bound_given ;	/**< @brief Whether bound was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int kmer_given ;	/**< @brief Whether kmer was given.  */
  unsigned int file_delimiter_given ;	/**< @brief Whether file-delimiter was given.  */
  unsigned int independent_probs_given ;	/**< @brief Whether independent-probs was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct SKA_cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure SKA_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure SKA_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *SKA_args_info_purpose;
/** @brief the usage string of the program */
extern const char *SKA_args_info_usage;
/** @brief the description string of the program */
extern const char *SKA_args_info_description;
/** @brief all the lines making the help output */
extern const char *SKA_args_info_help[];
/** @brief all the lines making the detailed help output (including hidden options and details) */
extern const char *SKA_args_info_detailed_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int SKA_cmdline_parser (int argc, char **argv,
  struct SKA_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use SKA_cmdline_parser_ext() instead
 */
int SKA_cmdline_parser2 (int argc, char **argv,
  struct SKA_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int SKA_cmdline_parser_ext (int argc, char **argv,
  struct SKA_args_info *args_info,
  struct SKA_cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int SKA_cmdline_parser_dump(FILE *outfile,
  struct SKA_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int SKA_cmdline_parser_file_save(const char *filename,
  struct SKA_args_info *args_info);

/**
 * Print the help
 */
void SKA_cmdline_parser_print_help(void);
/**
 * Print the detailed help (including hidden options and details)
 */
void SKA_cmdline_parser_print_detailed_help(void);
/**
 * Print the version
 */
void SKA_cmdline_parser_print_version(void);

/**
 * Initializes all the fields a SKA_cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void SKA_cmdline_parser_params_init(struct SKA_cmdline_parser_params *params);

/**
 * Allocates dynamically a SKA_cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized SKA_cmdline_parser_params structure
 */
struct SKA_cmdline_parser_params *SKA_cmdline_parser_params_create(void);

/**
 * Initializes the passed SKA_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void SKA_cmdline_parser_init (struct SKA_args_info *args_info);
/**
 * Deallocates the string fields of the SKA_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void SKA_cmdline_parser_free (struct SKA_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int SKA_cmdline_parser_required (struct SKA_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SKA_CMDL_H */
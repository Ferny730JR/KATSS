/* seqf_core.h - Header to access seqf internal definition and public seqf*
 * functions
 * 
 * Copyright (c) 2024-2025 Francisco F. Cavazos
 * Subject to the MIT License
 * 
 * This file should not be used in applications. It is used to implement the
 * seqf library and is subject to change.
 */

#ifndef SEQFCORE_H
#define SEQFCORE_H

#if (_HAS_ISA_L_ == 1)
#	include <igzip_lib.h>
#else
#	include <zlib.h>
#endif

#if SEQF_USE_C11_THREADS
#	include <threads.h>
#else
#	include "tinycthread.h"
#endif

#include "seqfile.h"

extern _Thread_local int seqferrno_;

typedef enum SEQF_COMPRESSION { GZIP, ZLIB, PLAIN } SEQF_COMPRESSION;

struct seqf_state {
	int fd;                       /** File descriptor */
	SEQF_COMPRESSION compression; /** Type of compression, if any */
	unsigned char mode;           /** File mode, e.g read, write, append, etc. */
	unsigned char type;           /** Type of file, e.g FASTA, FASTQ, or reads */
#if defined _IGZIP_H              /** Use isa-l if available, otherwise zlib */
	struct inflate_state stream;  /** ISA-L Decompressor */
#else
	z_stream stream;     /** ZLIB Decompressor */
	bool stream_is_init; /** Check is stream is initialized */
#endif

	unsigned char *in_buf;  /** Input buffer*/
	size_t in_bufsiz;       /** Size of the input buffer */
	unsigned char *out_buf; /** Output buffer */
	size_t out_bufsiz;      /** Size of the output buffer */
	unsigned char *next;    /** Next available byte in output buffer */
	size_t have;            /** Numberof bytes available in next */

	mtx_t mutex;        /** Mutex for thread safe functions */
	bool mutex_is_init; /** Check if mutex is initialized (for rnafclose) */

	bool eof; /** Flag to test if at end of rnafile */
};

typedef struct seqf_state *seqf_statep;

/* Define file mode macros */
#define SEQF_READ   0x00000001
#define SEQF_WRITE  0x00000002

#endif // SEQFCORE_H

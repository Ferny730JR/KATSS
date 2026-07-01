/* seqflib.c - seqf functions for opening a SeqFile stream
 *
 * Copyright (c) 2024-2025 Francisco F. Cavazos
 * Subject to the MIT License
 */

#include <stdbool.h>
#include <stdlib.h>

#include <fcntl.h>
#include <string.h>
#ifdef _WIN32
#	include <windows.h>
#	include <io.h>
#	define open     _open
#	define close    _close
#	define read     _read
#	define lseek    _lseek
#	define O_RDONLY _O_RDONLY
#	define O_WRONLY _O_WRONLY
#	define O_RDWR   _O_RDWR
#	define O_CREAT  _O_CREAT
#else
#	include <unistd.h>
#endif

#include "seqfile.h"
#include "seqf_core.h"
#include "seqf_read.h"

#define EXIT_AND_SETERR(state, _seqferrno)                                                         \
	do {                                                                                           \
		seqferrno_ = _seqferrno;                                                                   \
		seqfclose((SeqFile)state);                                                                 \
		return NULL;                                                                               \
	} while (0)

#define EXIT_ERR(state)                                                                            \
	do {                                                                                           \
		seqfclose((SeqFile)state);                                                                 \
		return NULL;                                                                               \
	} while (0)

static void
init_seqfstatep(seqf_statep state)
{
	state->fd          = -1;
	state->compression = PLAIN;
	state->mode        = 0;
	state->type        = 'd';
#ifndef _IGZIP_H
	state->stream_is_init = false;
#endif
	state->in_buf        = NULL;
	state->out_buf       = NULL;
	state->next          = NULL;
	state->have          = 0;
	state->mutex_is_init = false;
	state->eof           = false;
}

static int
getfdflags(const char *mode, int *flags)
{
	int open_flags = 0; // If file should be read/write
	int file_flags = 0; // How do file should be modified

	/* Parse the first char: can only be one of rwa */
	switch (*mode) {
	case 'r':
		open_flags = O_RDONLY;
		file_flags = 0; // no extra flags
		break;
	case 'w':
		open_flags = O_WRONLY;
		file_flags = O_CREAT | O_TRUNC;
		break;
	case 'a':
		open_flags = O_WRONLY;
		file_flags = O_CREAT | O_APPEND;
		break;
	default:
		seqferrno = SEQF_E_INVALM;
		return 1;
	}
	/* Check for [rwa]+, which opens file for reading and writing */
	while (*++mode) {
		switch (*mode) {
		case '+':
			open_flags = O_RDWR;
		default:
			break;
		}
	}

	*flags = open_flags | file_flags;
	return 0;
}

static int
getsfflags(seqf_statep state, const char *mode)
{
	/* Parse the first char: can only be one of rwa */
	switch (*mode++) {
	case 'r':
		state->mode = SEQF_READ;
		break;
	case 'w':
		state->mode = SEQF_WRITE;
		break;
	case 'a':
		state->mode = SEQF_WRITE;
		break;
	default:
		seqferrno = SEQF_E_INVALM;
		return 1;
	}

	while (*mode) {
		switch (*mode++) {
		case '+':
			state->mode |= SEQF_READ;
			state->mode |= SEQF_WRITE;
			break;
		case 'b':
			state->type = 'b';
			break;
		case 's':
			state->type = 's';
			break;
		case 'f':
			if (*mode == 'a')
				state->type = 'a';
			else if (*mode == 'q')
				state->type = 'q';
			else
				seqferrno = SEQF_E_INVALT;
			mode++;
			break;
		default:
			seqferrno = SEQF_E_INVALT;
			return 1;
		}
	}
	return 0;
}

static bool
hasflag(seqf_statep state, int flag)
{
	return (state->mode & flag) == flag;
}

static int
determine_compression(seqf_statep state)
{
	/* Get magic bytes */
	size_t nread = 0;
	do {
		size_t n = read(state->fd, state->in_buf, 2);
		if (n == -1) {
			seqferrno = SEQF_E_ERRNO;
			return 1;
		}
		if (n == 0)
			break; // reached EOF before reading magic bytes
		nread += n;
	} while (nread != 2);

	/* 0 or 1 bytes read, PLAIN text file */
	if (nread < 2) {
		state->compression = PLAIN;
	}

	/* Check magic bytes to determine type of compression used */
	if (state->in_buf[0] == 0x1F && state->in_buf[1] == 0x8B) {
		state->compression = GZIP;
	} else if (state->in_buf[0] == 0x78
	           && (state->in_buf[1] == 0x01 || state->in_buf[1] == 0x5E || state->in_buf[1] == 0x9C
	               || state->in_buf[1] == 0xDA)) {
		state->compression = ZLIB;
	} else {
		state->compression = PLAIN;
	}

	/* Rewind to beginning of file */
	lseek(state->fd, 0, SEEK_SET);
	return 0;
}

/* IUPAC nucleotide alphabet (upper and lower case) used to validate sequence
 * lines in determine_type(). Includes all standard ambiguity codes. */
static bool
is_iupac_nt(unsigned char c)
{
	switch (c | 0x20) { /* fold to lower case */
	case 'a':
	case 'c':
	case 'g':
	case 't':
	case 'u':
	case 'r':
	case 'y':
	case 's':
	case 'w':
	case 'k':
	case 'm':
	case 'b':
	case 'd':
	case 'h':
	case 'v':
	case 'n':
	case '-':
	case '.':
		return true;
	default:
		return false;
	}
}

/**
 * @brief Detect the sequence file type from the file's content and set
 * state->type accordingly.
 *
 * Detection rules:
 *   '>'  as first non-empty byte -> FASTA  (type = 'a')
 *   '@'  as first non-empty byte -> FASTQ  (type = 'q')
 *   IUPAC nucleotide character   -> sequences (type = 's')
 *   Empty file                   -> SEQF_E_EMPTY, return 1
 *   Anything else                -> SEQF_E_UNKNT, return 1
 *
 * For FASTA:  reads lines until it finds a header ('>'), then verifies that
 *   the very next non-empty line consists entirely of IUPAC characters.
 * For FASTQ:  verifies that a '@'-prefixed header is followed by an IUPAC
 *   sequence line, then '+', then a quality line of the same length.
 * For sequences: scans every non-empty line and confirms it contains only
 *   IUPAC nucleotide characters.
 *
 * File position is rewound to 0 before returning in all cases.
 *
 * @param state  SeqFile state
 * @return 0 on success (state->type set), 1 on error (seqferrno set).
 */
static int
determine_type(seqf_statep state)
{
	if (seqf_fetch(state) != 0)
		return 1; // error

	if (state->have == 0) {
		seqferrno = SEQF_E_EMPTY;
		return 1;
	}

	const unsigned char *buf = state->out_buf;
	const unsigned char *end = buf + state->have;

	/* Skip leading blank lines */
	while (buf < end && (*buf == '\n' || *buf == '\r' || *buf == ' ' || *buf == '\t'))
		buf++;
	if (buf == end) {
		seqferrno = SEQF_E_EMPTY;
		return 1;
	}

	/* Skip FASTA comment lines (uncommon, but check just in case) */
	while (*buf == ';') {
		while (buf < end && *buf != '\n')
			buf++;
		if (buf == end) { // very long comment line, probably not fasta
			seqferrno = SEQF_E_UNKNT;
			return 1;
		}
		buf++; // move past '\n'
	}

	unsigned char first = *buf;

	/* --- FASTA --- */
	if (first == '>') {
		while (buf < end && *buf != '\n')
			buf++;
		if (buf == end) {
			state->type = 'a';
			return 0;
		}
		buf++;

		/* Skip blank lines after header */
		while (buf < end && (*buf == '\n' || *buf == '\r'))
			buf++;
		if (buf == end) {
			state->type = 'a';
			return 0;
		}

		/* Allow leading spaces before the sequence */
		while (buf < end && (*buf == ' ' || *buf == '\t'))
			buf++;

		if (buf < end && !is_iupac_nt(*buf)) {
			seqferrno = SEQF_E_UNKNT;
			return 1;
		}
		state->type = 'a';
		return 0;
	}

	/* --- FASTQ --- */
	if (first == '@') {
		/* Skip header line */
		while (buf < end && *buf != '\n')
			buf++;
		if (buf == end) {
			state->type = 'q';
			return 0;
		}
		buf++;

		/* Sequence line */
		const unsigned char *seq_start = buf;
		while (buf < end && *buf != '\n' && *buf != '\r') {
			if (!is_iupac_nt(*buf)) {
				seqferrno = SEQF_E_UNKNT;
				return 1;
			}
			buf++;
		}
		size_t seq_len = (size_t)(buf - seq_start);

		if (buf == end || seq_len == 0) {
			state->type = 'q'; // very long sequence, assume fastq
			return 0;
		}

		/* Skip newlines */
		while (buf < end && (*buf == '\n' || *buf == '\r'))
			buf++;

		/* Check for '+' */
		if (buf >= end) {
			state->type = 'q';
			return 0;
		}
		if (*buf != '+') {
			seqferrno = SEQF_E_UNKNT;
			return 1;
		}

		/* Skip '+' line */
		while (buf < end && *buf != '\n')
			buf++;
		if (buf == end) {
			state->type = 'q';
			return 0;
		}
		buf++;

		/* Quality line */
		size_t qual_len = 0;
		while (buf < end && *buf != '\n' && *buf != '\r') {
			if (*buf < 0x21 || *buf > 0x7E) { // not alphanumeric
				seqferrno = SEQF_E_UNKNT;
				return 1;
			}
			qual_len++;
			buf++;
		}

		/* If we hit the end of the buffer, we can't expect qual_len == seq_len */
		if (buf < end && qual_len != seq_len) {
			seqferrno = SEQF_E_UNKNT;
			return 1;
		}

		state->type = 'q';
		return 0;
	}

	/* --- SEQUENCES (Raw) --- */
	if (is_iupac_nt(first)) {
		const unsigned char *p = buf;
		while (p < end) {
			if (*p == '\n' || *p == '\r' || *p == ' ' || *p == '\t') {
				p++;
				continue;
			}
			if (!is_iupac_nt(*p)) {
				seqferrno = SEQF_E_UNKNT;
				return 1;
			}
			p++;
		}
		state->type = 's';
		return 0;
	}

	seqferrno = SEQF_E_UNKNT;
	return 1;
}

SeqFile
seqfdopen(int fd, const char *mode)
{
	seqferrno_ = 0; // no error encountered. yet.

	/* Initialize SeqFile */
	seqf_statep seq_file = malloc(sizeof *seq_file);
	if (seq_file == NULL)
		EXIT_AND_SETERR(seq_file, SEQF_E_OOM);

	/* Init deafult values */
	init_seqfstatep(seq_file);

	/* Open file and check for errors */
	seq_file->fd = fd;
	if (fd < 0)
		EXIT_AND_SETERR(seq_file, SEQF_E_ERRNO);

	/* Set the file flags */
	if (getsfflags(seq_file, mode) != 0)
		EXIT_ERR(seq_file);

	/* Initialize mutex */
	if (mtx_init(&seq_file->mutex, mtx_plain) != thrd_success)
		EXIT_AND_SETERR(seq_file, SEQF_E_MUTEX);
	seq_file->mutex_is_init = true;

	/* Set input and output buffers */
	seq_file->in_buf = malloc(SEQFBUFSIZ * sizeof *seq_file->in_buf);
	if (seq_file->in_buf == NULL)
		EXIT_AND_SETERR(seq_file, SEQF_E_OOM);
	seq_file->in_bufsiz = SEQFBUFSIZ;
	seq_file->out_buf   = malloc(2 * SEQFBUFSIZ * sizeof *seq_file->out_buf);
	if (seq_file->out_buf == NULL)
		EXIT_AND_SETERR(seq_file, SEQF_E_OOM);
	seq_file->out_bufsiz = 2 * SEQFBUFSIZ;
	seq_file->next       = seq_file->out_buf;

	/* Determine type of compression, if any */
	if (hasflag(seq_file, SEQF_READ) && determine_compression(seq_file) != 0)
		EXIT_ERR(seq_file);

	/* Initialize decompressor */
	if (hasflag(seq_file, SEQF_READ) && seq_file->compression != PLAIN) {
#if defined _IGZIP_H
		isal_inflate_init(&seq_file->stream);
		seq_file->stream.crc_flag = seq_file->compression == GZIP ? ISAL_GZIP : ISAL_ZLIB;
		seq_file->stream.next_in  = seq_file->in_buf;
#else
		/* allocate inflate state */
		int ret                   = Z_ERRNO;
		seq_file->stream.zalloc   = Z_NULL;
		seq_file->stream.zfree    = Z_NULL;
		seq_file->stream.opaque   = Z_NULL;
		seq_file->stream.avail_in = 0;
		seq_file->stream.next_in  = Z_NULL;
		if (seq_file->compression == GZIP) {
			ret = inflateInit2(&seq_file->stream, 16 + MAX_WBITS);
		} else if (seq_file->compression == ZLIB) {
			ret = inflateInit(&seq_file->stream);
		} else {
			seqferrno_ = 1;
			seqfclose((SeqFile)seq_file);
		}
		if (ret != Z_OK)
			EXIT_AND_SETERR(seq_file, 1);
		seq_file->stream.next_in = seq_file->in_buf;
		seq_file->stream_is_init = true;
#endif
	}

	/* Auto-detect file type if the caller did not specify one */
	if (hasflag(seq_file, SEQF_READ) && seq_file->type == 'd') {
		if (determine_type(seq_file) != 0)
			EXIT_ERR(seq_file);
	}

	return (SeqFile)seq_file;
}

SeqFile
seqfopen(const char *path, const char *mode)
{
	/* Validate input */
	if (path == NULL || mode == NULL)
		return NULL;

	int flags;
	if (getfdflags(mode, &flags) != 0)
		return NULL;

	/* Currently, only read-only is supported */
	if (flags != O_RDONLY) {
		seqferrno = SEQF_E_INVALM;
		return NULL;
	}
#ifdef _WIN32
	flags |= O_BINARY;
#endif
	int fd = open(path, flags); // currently only support reading
	if (fd == -1) {
		seqferrno_ = SEQF_E_ERRNO;
		return NULL;
	}
	return seqfdopen(fd, mode);
}

int
seqfclose(SeqFile file)
{
	if (file == NULL)
		return 1;
	int return_code   = 0;
	seqf_statep state = (seqf_statep)file;
	if (state->fd > 2 && close(state->fd) == -1)
		return_code = seqferrno_ = 1;
	if (state->mutex_is_init)
		mtx_destroy(&state->mutex);
	if (state->in_buf)
		free(state->in_buf);
	if (state->out_buf)
		free(state->out_buf);
#ifndef _IGZIP_H
	if (state->stream_is_init)
		inflateEnd(&state->stream);
#endif
	free(state);
	return return_code;
}

int
seqfrewind(SeqFile file)
{
	if (file == NULL)
		return -1;
	seqf_statep state = (seqf_statep)file;
	if (lseek(state->fd, 0, SEEK_SET) == -1) {
		seqferrno_ = 1;
		return -1;
	}
	state->have = 0;
	state->eof  = false;
#if defined _IGZIP_H
	isal_inflate_reset(&state->stream);
	state->stream.crc_flag = state->compression == GZIP ? ISAL_GZIP : ISAL_ZLIB;
	state->stream.next_in  = state->in_buf;
#else
	if (state->stream_is_init) {
		int ret;
		if (state->compression == GZIP)
			ret = inflateReset2(&state->stream, 16 + MAX_WBITS);
		else if (state->compression == ZLIB)
			ret = inflateReset(&state->stream);
		else
			return -1;
		if (ret != Z_OK)
			return -1;
		state->stream.next_in  = state->in_buf;
		state->stream.avail_in = 0;
	}
#endif
	return 0;
}

bool
seqfeof(SeqFile file)
{
	return ((seqf_statep)file)->eof;
}

int
seqfsetibuf(SeqFile file, size_t bufsize)
{
	if (file == NULL)
		return -1;
	seqf_statep state = (seqf_statep)file;

	unsigned char *t  = realloc(state->in_buf, bufsize);
	if (t == NULL)
		return -1;
	state->in_buf    = t;
	state->in_bufsiz = bufsize;
	return 0;
}

int
seqfsetobuf(SeqFile file, size_t bufsize)
{
	if (file == NULL)
		return -1;
	seqf_statep state = (seqf_statep)file;

	unsigned char *t  = realloc(state->out_buf, bufsize);
	if (t == NULL)
		return -1;
	state->out_buf    = t;
	state->out_bufsiz = bufsize;
	return 0;
}

int
seqfsetbuf(SeqFile file, size_t bufsize)
{
	if (seqfsetibuf(file, bufsize) != 0)
		return -1;
	if (seqfsetobuf(file, bufsize << 1) != 0)
		return -2;
	return 0;
}

int
seqfsettype(SeqFile file, char filetype)
{
	switch (filetype) {
	case 'a':
	case 'q':
	case 's':
	case 'b':
		((seqf_statep)file)->type = filetype;
		return 0;
	default:
		return 1;
	}
}

unsigned char
seqftype(SeqFile file)
{
	return ((seqf_statep)file)->type;
}

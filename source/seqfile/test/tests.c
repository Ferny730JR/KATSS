#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "minunit.h"

#include "seqf_core.h"
#include "seqfile.h"

#define STRINGIZE(arg) #arg
#define TXT2STR(arg)   STRINGIZE(arg)

static UTEST_TYPE
test_seqfopen(void)
{
	init_unit_tests("Testing seqfopen");
	SeqFile file;

#define assert_file(file, exp_compression, exp_type)                                               \
	file != NULL && ((seqf_statep)file)->compression == exp_compression                            \
		&& ((seqf_statep)file)->type == exp_type && ((seqf_statep)file)->mutex_is_init

	file = seqfopen(TXT2STR(EXAMPLE_FASTA), "rfa");
	mu_assert("Open fasta file", assert_file(file, PLAIN, 'a'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTA_GZ), "rfa");
	mu_assert("Open compressed fasta file", assert_file(file, GZIP, 'a'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTQ), "rfq");
	mu_assert("Open fastq file", assert_file(file, PLAIN, 'q'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTQ_GZ), "rfq");
	mu_assert("Open compressed fastq file", assert_file(file, GZIP, 'q'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_READS), "rs");
	mu_assert("Open sequences file", assert_file(file, PLAIN, 's'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_READS_GZ), "rs");
	mu_assert("Open compressed sequences file", assert_file(file, GZIP, 's'));
	seqfclose(file);

	/* Check that return NULL on files that don't exist */
	file = seqfopen("non-existent-dir/non-existent-file", NULL);
	mu_assert("Opening non-existent file", file == NULL);
	seqfclose(file);

	file = seqfopen("example_files/example.reads", "wrong!");
	mu_assert("Open file with unsupported mode", file == NULL);

#undef assert_file
	unit_tests_end;
}

static UTEST_TYPE
test_seqfclose(void)
{
	init_unit_tests("Testing seqfclose");
	mu_assert("Close null file", seqfclose(NULL) == 1);
	mu_assert("Close opened file", seqfclose(seqfopen(TXT2STR(EXAMPLE_READS), "rs")) == 0);

	unit_tests_end;
}

static UTEST_TYPE
test_seqferrno(void)
{
	init_unit_tests("Testing seqferrno");

	mu_assert("seqferrno has not been set", seqferrno == 0);

	seqfopen("non-existent-path/no-existent-file", "r");
	mu_assert("Opening file that does not exist", seqferrno == SEQF_E_ERRNO);

	seqfopen(TXT2STR(EXAMPLE_READS), "wrong mode!");
	mu_assert("Error code for wrong mode", seqferrno == SEQF_E_INVALM);

	unit_tests_end;
}

static UTEST_TYPE
test_seqfgetc(void)
{
	init_unit_tests("Testing seqfgetc");

	/* Prepare variables to use */
	bool passed  = false;
	SeqFile file = seqfopen(TXT2STR(EXAMPLE_READS), "rs");

	passed       = seqfgetc(file) == 'G';
	mu_assert("getc call on uncompressed file", passed);

	seqfrewind(file);
	passed = true;
	char *reads_seq =
		"GCATACGGTGAAAGCTCAGCTTTCCAGCGCTGCTTTACAGTTGGCACGATTAACCCAAGAACGTTATTTCTGTCAAATTTTGAGTTGGTTGTGGGCAAGG";
	char *tmp = reads_seq;
	for (int i = 0; i < 100; i++) {
		if (*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	mu_assert("Many seqfgetc calls", passed);

	seqfclose(file);
	file   = seqfopen(TXT2STR(EXAMPLE_READS_GZ), "rs");
	passed = seqfgetc(file) == 'G';
	mu_assert("getc call on compressed file", passed);

	seqfrewind(file);
	passed = true;
	tmp    = reads_seq;
	for (int i = 0; i < 100; i++) {
		if (*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	seqfclose(file);
	mu_assert("Many seqfgetc calls on compressed file", passed);

	char buf[400];
	tmp       = buf;
	FILE *fp  = fopen(TXT2STR(EXAMPLE_FASTA), "rb");
	file      = seqfopen(TXT2STR(EXAMPLE_FASTA), "rb");
	int nread = fread(buf, 1, sizeof(buf), fp);
	for (int i = 0; i < nread; i++) {
		if (*tmp++ != seqfgetc(file)) {
			passed = false;
			break;
		}
	}
	fclose(fp);
	mu_assert("seqfgetc entire file", passed);

	passed = seqfgetc(file) == EOF;
	mu_assert("seqfgetc when reached end of file", passed);
	seqfclose(file);

	unit_tests_end;
}

static SeqFile
open_tmpfile_with_content(const char *content, size_t len)
{
	char path[] = "/tmp/seqf_test_XXXXXX";
	int fd      = mkstemp(path);
	if (fd == -1)
		return NULL;

	/* Write content then remove the directory entry so the fd owns it */
	size_t written = 0;
	while (written < len) {
		ssize_t n = write(fd, content + written, len - written);
		if (n == -1) {
			close(fd);
			unlink(path);
			return NULL;
		}
		written += (size_t)n;
	}
	unlink(path); /* file lives until fd is closed */
	lseek(fd, 0, SEEK_SET);

	/* "r" mode with no type specifier -> triggers auto-detection */
	return seqfdopen(fd, "r");
}

static UTEST_TYPE
test_determine_type(void)
{
	init_unit_tests("Testing determine_type (auto-detection)");
	SeqFile file;

/* Convenience: check type field of a successfully opened SeqFile */
#define assert_type(sf, expected_type)                                                             \
	((sf) != NULL && ((seqf_statep)(sf))->type == (expected_type))

	/* ---- FASTA -------------------------------------------------------- */
	file = seqfopen(TXT2STR(EXAMPLE_FASTA), "r");
	mu_assert("Auto-detect plain FASTA file", assert_type(file, 'a'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTA_GZ), "r");
	mu_assert("Auto-detect gzip FASTA file", assert_type(file, 'a'));
	seqfclose(file);

	/* Minimal FASTA: single-char sequence */
	file = open_tmpfile_with_content(">s\nA\n", 5);
	mu_assert("Minimal FASTA (1-nt sequence)", assert_type(file, 'a'));
	seqfclose(file);

	/* FASTA with IUPAC ambiguity codes in sequence */
	file = open_tmpfile_with_content(">hdr\nACGTRYSWKMBDHVN\n", 18);
	mu_assert("FASTA with full IUPAC alphabet", assert_type(file, 'a'));
	seqfclose(file);

	/* FASTA with multiple records */
	const char *multi_fasta = ">seq1\nACGT\n>seq2\nTGCA\n";
	file                    = open_tmpfile_with_content(multi_fasta, strlen(multi_fasta));
	mu_assert("FASTA with multiple records", assert_type(file, 'a'));
	seqfclose(file);

	/* FASTA header-only (no sequence body) — still detected as FASTA */
	file = open_tmpfile_with_content(">header_only\n", 13);
	mu_assert("FASTA with header only", assert_type(file, 'a'));
	seqfclose(file);

	/* ---- FASTQ -------------------------------------------------------- */
	file = seqfopen(TXT2STR(EXAMPLE_FASTQ), "r");
	mu_assert("Auto-detect plain FASTQ file", assert_type(file, 'q'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_FASTQ_GZ), "r");
	mu_assert("Auto-detect gzip FASTQ file", assert_type(file, 'q'));
	seqfclose(file);

	/* Minimal FASTQ: one record */
	const char *min_fastq = "@read1\nACGT\n+\nIIII\n";
	file                  = open_tmpfile_with_content(min_fastq, strlen(min_fastq));
	mu_assert("Minimal FASTQ (1 record)", assert_type(file, 'q'));
	seqfclose(file);

	/* FASTQ: One NT record */
	const char *one_fastq = "@r\nA\n+\nI\n";
	file = open_tmpfile_with_content(one_fastq, strlen(one_fastq));
	mu_assert("One NT FASTQ", assert_type(file, 'q'));
	seqfclose(file);

	/* FASTQ: One NT record */
	const char *mult_one_fastq = "@r\nA\n+\nI\n@r\nA\n+\nI\n@r\nA\n+\nI\n@r\nA\n+\nI\n";
	file = open_tmpfile_with_content(one_fastq, strlen(mult_one_fastq));
	mu_assert("Multi Record One NT FASTQ", assert_type(file, 'q'));
	seqfclose(file);

	/* FASTQ with multi-char quality scores */
	const char *multi_fastq = "@r1\nACGTACGT\n+\n!~ABCDEF\n"
							  "@r2\nTGCA\n+\nIIII\n";
	file                    = open_tmpfile_with_content(multi_fastq, strlen(multi_fastq));
	mu_assert("FASTQ with multiple records", assert_type(file, 'q'));
	seqfclose(file);

	/* ---- Sequences ---------------------------------------------------- */
	file = seqfopen(TXT2STR(EXAMPLE_READS), "r");
	mu_assert("Auto-detect plain sequences file", assert_type(file, 's'));
	seqfclose(file);

	file = seqfopen(TXT2STR(EXAMPLE_READS_GZ), "r");
	mu_assert("Auto-detect gzip sequences file", assert_type(file, 's'));
	seqfclose(file);

	/* Single nucleotide — smallest possible sequences file */
	file = open_tmpfile_with_content("G\n", 2);
	mu_assert("Sequences file: single nucleotide", assert_type(file, 's'));
	seqfclose(file);

	/* Single nucleotide, no trailing newline */
	file = open_tmpfile_with_content("A", 1);
	mu_assert("Sequences file: 1 byte, no newline", assert_type(file, 's'));
	seqfclose(file);

	/* Multiple sequences, mixed IUPAC */
	const char *seqs = "ACGTN\nRYSWKM\nBDHV\n";
	file             = open_tmpfile_with_content(seqs, strlen(seqs));
	mu_assert("Sequences file: multiple IUPAC lines", assert_type(file, 's'));
	seqfclose(file);

	/* Lowercase nucleotides */
	file = open_tmpfile_with_content("acgtn\n", 6);
	mu_assert("Sequences file: lowercase nucleotides", assert_type(file, 's'));
	seqfclose(file);

	/* ---- Error cases -------------------------------------------------- */

	/* Empty file -> SEQF_E_EMPTY */
	seqferrno = 0;
	file      = open_tmpfile_with_content("", 0);
	mu_assert("Empty file returns NULL", file == NULL);
	mu_assert("Empty file sets SEQF_E_EMPTY", seqferrno == SEQF_E_EMPTY);
	seqfclose(file);

	/* Only newlines -> SEQF_E_EMPTY */
	seqferrno = 0;
	file      = open_tmpfile_with_content("\n\n\n", 3);
	mu_assert("Newlines-only file returns NULL", file == NULL);
	mu_assert("Newlines-only sets SEQF_E_EMPTY", seqferrno == SEQF_E_EMPTY);
	seqfclose(file);

	/* Garbage content -> SEQF_E_UNKNT */
	seqferrno = 0;
	file      = open_tmpfile_with_content("!invalid content\n", 17);
	mu_assert("Garbage content returns NULL", file == NULL);
	mu_assert("Garbage content sets SEQF_E_UNKNT", seqferrno == SEQF_E_UNKNT);
	seqfclose(file);

	/* FASTQ-like but quality length mismatch -> SEQF_E_UNKNT */
	seqferrno             = 0;
	const char *bad_fastq = "@r\nACGT\n+\nII\n"; /* qual shorter than seq */
	file                  = open_tmpfile_with_content(bad_fastq, strlen(bad_fastq));
	mu_assert("FASTQ qual/seq length mismatch returns NULL", file == NULL);
	mu_assert("FASTQ qual/seq mismatch sets SEQF_E_UNKNT", seqferrno == SEQF_E_UNKNT);
	seqfclose(file);

	/* FASTQ-like but missing '+' separator -> SEQF_E_UNKNT */
	seqferrno           = 0;
	const char *no_plus = "@r\nACGT\nACGT\nIIII\n";
	file                = open_tmpfile_with_content(no_plus, strlen(no_plus));
	mu_assert("FASTQ without '+' separator returns NULL", file == NULL);
	mu_assert("FASTQ without '+' sets SEQF_E_UNKNT", seqferrno == SEQF_E_UNKNT);
	seqfclose(file);

	/* Sequences file with an invalid character mid-file -> SEQF_E_UNKNT */
	seqferrno           = 0;
	const char *bad_seq = "ACGT\nACG!T\nACGT\n";
	file                = open_tmpfile_with_content(bad_seq, strlen(bad_seq));
	mu_assert("Sequences file with invalid char returns NULL", file == NULL);
	mu_assert("Sequences invalid char sets SEQF_E_UNKNT", seqferrno == SEQF_E_UNKNT);
	seqfclose(file);

#undef assert_type
	unit_tests_end;
}

static void
all_tests(void)
{
	init_run_test;

	/* Begin tests */
	mu_run_test(test_seqfopen);
	mu_run_test(test_seqfclose);
	mu_run_test(test_seqferrno);
	mu_run_test(test_determine_type);
	mu_run_test(test_seqfgetc);

	/* End of tests */
	run_test_end;
}

int
main(void)
{
	all_tests();
	return 0;
}

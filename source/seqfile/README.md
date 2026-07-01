![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![GitHub release](https://img.shields.io/github/v/release/Ferny730JR/seqfile)
![GitHub issues](https://img.shields.io/github/issues/Ferny730JR/seqfile)
![CMake >=3.10.0](https://img.shields.io/badge/CMake-%3E=3.10.0-blue)

# SeqFile

**SeqFile** is a lightweight, zlib-inspired C library for efficient reading and
processing of FASTA, FASTQ, and raw sequence file formats.

It is designed to be fast, minimal, thread-safe, and easy to integrate into C
projects. Whether you're building a full-scale bioinformatics pipeline or a
small utility, SeqFile provides a clean and efficient interface to handle
sequence data, mimicking the <stdio.h> and <zlib.h> API's for familiarity.

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [Examples](#examples)
4. [License](#license)
5. [Contact](#contact)

## Installation

### Dependencies

- [zlib](https://zlib.net/) (for compressed input support)
- [CMake ≥ 3.9.0](https://cmake.org/) (build configuration)

### Quick Start

Once you have downloaded the repository and the required dependencies, run the
following commands:

```bash
cd seqfile
cmake -B build && cmake --build build
cmake --build build --target install
```

**Note**: This will install the library into the default system path, which is
typically `/usr/local/lib` and `/usr/local/include`.

### User-dir Installation

One issue you might encounter with the previous installation is that it might 
require root privileges. In case you do not have root privileges in your
computer, you can specify the installation location to a place you have write
access to by using the ```-DCMAKE_INSTALL_PREFIX``` flag as such:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/your/preferred/path -B build && cmake --build build
cmake --build build --target install
```

### Testng

Once you have setup the project, you can run the tests to ensure SeqFile is working as intended.

```bash
./build/test/test-seqf
```

## Usage

To use SeqFile in your own project, include the SeqFile header as follows in
your own code:
```c
#include <seqfile.h>
```

Then, you have two options to link the library:

### 1. Manual Linking

Compile your C file with:

```bash
gcc your_program.c -lseqf -o your_program
```

### 2. Using CMake

You can add SeqFile to your CMake-based project like this:

```cmake
add_subdirectory(seqfile)
target_link_libraries(your_target PRIVATE seqf_static)  # for static linking
# or
target_link_libraries(your_target PRIVATE seqf_shared)  # for shared linking
```

## Examples

Here's a basic example of how to use the library:

```c
#include <seqfile.h>
#include <stdio.h> // for printf
#include <stdlib.h> // for malloc and size_t

int
main(void)
{
	/* Open a fasta file for reading */
	SeqFile sfp = seqfopen("example.fasta", "a");
	if(sfp == NULL) {
		printf("seqfopen: %s", seqfstrerror(seqferrno));
		return 1;
	}

	/* Allocate space for buffer to read into */
	size_t bufsize = 1000;
	char *buffer = malloc(bufsize * sizeof *buffer);

	/* Read all sequences in the file */
	while(seqfgets(sfp, buffer, bufsize) != NULL) {
		printf("sequence: '%s'\n", buffer);
	}

	/* Close the SeqFile handle to release resources */
	if(seqfclose(sfp) != 0) {
		printf("seqfclose: %s\n", seqfstrerror(seqferrno));
		return 1;
	}

	/* Success */
	return 0;
}
```

## License

This project is licensed under the [MIT License](LICENSE).

## Contact
For bug reports, feature requests, or general questions:

- Open an [issue on GitHub](https://github.com/Ferny730JR/seqfile/issues)
- Or reach out via email: **ffcavazos@miners.utep.edu**

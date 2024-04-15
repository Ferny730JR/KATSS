<p align="center">
	<img src="./res/KATSS_logo.png" alt="KATSS Logo of a blue cat with a RNA-shaped tail." width="300">
</p>

# KATSS - K-mer Analysis Tools for Sequence and Structure

A C package containing programs to analyze RNA-protein interaction.

## Introduction

KATSS is a collection of C programs designed for the analysis of RNA sequences, particularly focusing on RNA-protein interactions. The package includes two main tools: `bppPipeline` and `SKA`, each tailored for specific analyses within the realm of RNA sequence analysis.

The tools allow you to:
* Calculate the motif of RNA-binding proteins
* Determine the secondary motifs
* Find the nucleotide and length preference of varying motifs
* Search for specific sequences
* Predict the binding preference of a protein
* TODO: Predict how G4's affects RNA stability

KATSS uses the ViennaRNA Package for the ```bppPipeline``` program, as such it is requires that you have it installed if you plan on using this.

## Table of Contents
1. [Installation](#installation)
2. [Executable Program](#executable-programs)
3. [Example Usage](#example-usage)
4. [License](#license)
## Installation

KATSS uses CMake to install the programs. Alongside, KATSS it is necessary to have the ViennaRNA library installed if you plan on using ```bppPipeline```. The following instructions are to install KATSS from source.

### Quick Start

Once you download the project, you can simply do the following:

```bash
cd KATSS
cmake .
make
sudo make install
```

**Note**: This will install the binaries into the default system path, which is typically `/usr/local/bin`.

### User-dir Installation

One issue you might encounter with the previous installation is that it requires root privileges. In case you do not have root privileges in your computer, you can specify the installation location to a place you do by using the ```-DCMAKE_INSTALL_PREFIX``` flag as such:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/your/preferred/path .
make
make install
```

Ensure that the specified path is included in your system's `PATH` environment variable to execute the binaries conveniently from any location.

## Executable Programs

The KATSS package includes the following executable programs:
| Program       | Description                                                    |
| ------------- | :--------------------------------------------------------------|
| `SKA`         | Compute enriched motifs in RNA-sequences                       |
| `bppPipeline` | Predict enriched base-pair interactions in protein-bound RNA   |

### `bppPipeline`

```bash
bppPipeline [OPTIONS] [<input.fa>] [<bound.fa>]
```

#### Options:

- `-h`, `--help`: Print help and exit.
- `--detailed-help`: Print detailed help, including all options, and exit.
- `-V`, `--version`: Print version and exit.

#### I/O Options:

- `-i`, `--input=filename`: Specify the control sequences file.
- `-b`, `--bound=filename`: Specify the protein-bound RNA sequences file.
- `-o`, `--output=filename`: Set the name of the output files.
- `-k`, `--kmer=INT`: Set the length of k-mers.
- `-d`, `--file-delimiter=CHAR`: Set the delimiter used in the output file.
- `--keepFolds`: Keep the base-pair probability of the folded sequences in a file.
- `--input-fold=filename`: Set the file for the control sequences with pre-calculated base-pair probabilities.
- `--bound-fold=filename`: Set the file for the protein-bound sequences with pre-calculated base-pair probabilities.
- `-j`, `--jobs[=number]`: Process files in parallel using multiple threads.

#### Algorithms:

- `-w`, `--seq-windows[=INT]`: Split the sequence into sliding windows and find the mean probability per position in the window.
- `--bin`: Produce bin file.
- `--frq`: Keep the frequency files.

### `SKA`

```bash
SKA [OPTIONS] [<input.fa>] [<bound.fa>]
```

#### Options:

- `-h`, `--help`: Print help and exit.
- `--detailed-help`: Print detailed help, including all options, and exit.
- `-V`, `--version`: Print version and exit.

#### I/O Options:

- `-i`, `--input=filename`: Specify the control sequences file.
- `-b`, `--bound=filename`: Specify the protein-bound RNA sequences file.
- `-o`, `--output=filename`: Set the name of the output files.
- `-k`, `--kmer=INT`: Set the length of k-mers.
- `-r`, `--iterations=INT`: Set the number of iterations for SKA.
- `-d`, `--file-delimiter=char`: Set the delimiter used in the output file.

#### Algorithms:

- `-p`, `--independent-probs`: Calculate the enrichments without the input reads.
- `--fmotif=string`: Search for a specific fixed motif, rather than all k-mers.
- `--motif=pattern`: Search for a specific motif using regular expression patterns.

## Example Usage

### `bppPipeline`

```bash
bppPipeline -i control_sequences.fa -b bound_sequences.fa -o output.csv -k 3
```

### `SKA`

```bash
SKA -i control_sequences.fa -b bound_sequences.fa -o output.csv -k 5
```

## License

I still need to get a license.
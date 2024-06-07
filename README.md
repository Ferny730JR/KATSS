<p align="center">
	<img src="./res/KATSS_logo.png" alt="KATSS Logo of a blue cat with a RNA-shaped tail." width="300">
</p>

# KATSS - K-mer Analysis Tools for Sequence and Structure

A C package containing programs to analyze RNA-protein interaction.

## Introduction

KATSS is a collection of C programs designed for the analysis of RNA sequences, particularly focusing on RNA-protein interactions. The package includes three main programs: `bppPipeline`, `SKA`, and `aspiire`, each tailored for specific analyses.

The tools allow you to:
* Calculate the motif of RNA-binding proteins
* Determine the secondary motifs
* Find the nucleotide and length preference of varying motifs
* Predict the binding preference of a protein
* Predict for iron regulatory elements

## Table of Contents
1. [Installation](#installation)
2. [Executable Program](#executable-programs)
3. [Example Usage](#example-usage)
4. [License](#license)
## Installation

KATSS uses CMake to install the programs. In order to compile from source, the following dependencies are required:
* ```ViennaRNA```, for RNA structure prediction
* ```zlib```, to read from compressed files

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

One issue you might encounter with the previous installation is that it requires root privileges. In case you do not have root privileges in your computer, you can specify the installation location to a place you have write access to by using the ```-DCMAKE_INSTALL_PREFIX``` flag as such:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/your/preferred/path .
make
make install
```

Ensure that the specified path is included in your system's `PATH` environment variable to execute the binaries conveniently from any location.

## Executable Programs

The KATSS package includes the following executable programs:
| Program       | Description                                                    |
| ------------- | -------------------------------------------------------------- |
| `SKA`         | Compute enriched motifs in RNA-sequences                       |
| `bppPipeline` | Predict enriched base-pair interactions in protein-bound RNA   |
| `aspiire`     | Predict for iron regulatory elements in a given sequence       |

To get more information about a program, run them with the `--help` or `--detailed-help` flag as such:

```bash
SKA --detailed-help
```

## Example Usage

### `bppPipeline`

```bash
bppPipeline -i control_sequences.fastq.gz -b bound_sequences.fastq.gz -o output.csv -k 3
```

### `SKA`

```bash
SKA -i control_sequences.fa -b bound_sequences.fa -o output.csv --kmer=5 --iterations=10
```

### `aspiire`
```bash
aspiire --input="human_utr3.fa" --output="human_utr3_ires"
```

## License

I still need to get a license.
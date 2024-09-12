
# BLAST-Python-Implementation

**BLAST-Python-Implementation** is a Python-based implementation of the Basic Local Alignment Search Tool (BLAST), designed for a Master's 2 Bioinformatics project at **Université Paris-Cité**. This tool facilitates searching for sequence similarity within a sequence database. It supports parallelization for improved speed and efficiency and provides results comparable to BLAST+, the C implementation.

## Features

- **K-mer extraction and indexing**: Efficient k-mer-based sequence search.
- **Double-hit detection**: Identifies hits between sequences on the same diagonal.
- **Smith-Waterman local alignment**: Extends alignments with scoring matrices and gap penalties.
- **BLOSUM62 scoring**: Implements substitution matrix scoring for alignments.
- **Parallel execution**: Supports multi-core CPU parallelism for faster sequence searching.
- **E-value calculation**: Estimates statistical significance of sequence alignments.
- **Tested and validated**: Compares results with BLAST+ using statistical tests (Mann-Whitney).

## Prerequisites

Make sure you have the following installed before running the project:

- **Conda**: A Python environment management tool.
- **Python 3.8+**: The project is built and tested on Python 3.8 and above.
- **Biopython**: Used for sequence parsing and BLOSUM matrix handling.
- **NumPy**: Essential for array and matrix calculations.

## Download this Repository

Clone the repository from GitHub:

```bash
git clone https://github.com/theo-serralta/BLAST-Python-Implementation.git
cd BLAST-Python-Implementation
```

## Install Dependencies

### Conda Environment Setup

1. Install [conda](https://docs.conda.io/en/latest/miniconda.html).
2. Create and activate the environment using the `blast_python.yml` file:

```bash
conda env create -f blast_python.yml
conda activate blast_python
```

### Manual Dependency Installation

Alternatively, you can manually install the required packages:

```bash
conda create -n blast_python python=3.8
conda activate blast_python
pip install biopython numpy
```

## Usage

### Running the BLAST Implementation

1. Place your query and target sequences in FASTA format within the appropriate directory.
2. Use the following command to run the alignment:

```bash
python src/blast_python.py
```

3. To change query sequences or database settings, modify the `src/blast_python.py` file, providing the paths to your query and database sequences.

### Example Usage

For example, to align the following query sequence:

```plaintext
>ADEPILVA
ADEPILVA
```

With a target sequence database, execute the script as described in the previous steps.

### Parallel Processing

To optimize performance, the program supports multi-core processing. You can adjust the number of threads to improve speed by modifying the relevant parameters in the `blast_python.py` script.

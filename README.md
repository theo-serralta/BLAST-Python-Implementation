
# BLAST-Python-Implementation

BLAST-Python-Implementation is a Python implementation of the Basic Local Alignment Search Tool (BLAST) for a Master's 2 project in Bioinformatics at Université Paris-Cité. The program enables the search for similar sequences within a sequence database using a custom BLASTP algorithm, as well as integration with BLASTN and BLASTX through Biopython. It is open-source and available under the MIT license.

## Project Overview

This project implements a custom BLASTP algorithm along with options to run BLASTX and BLASTN using Biopython's integration with the BLAST+ command line tools. The custom BLASTP algorithm includes:
- **K-mer extraction**
- **Double-hit detection**
- **Gapped alignment extension**

Additionally, a GUI version is provided for ease of use, allowing users to configure and run BLAST from a graphical interface, while the CLI version enables command-line execution.

### Filename:
**blast_python.py**

### Description:
This script implements a custom BLASTP algorithm along with options to run BLASTX and BLASTN using Biopython's integration with the BLAST+ command line tools. The custom BLASTP algorithm includes k-mer extraction, double-hit detection, and gapped alignment extension. The GUI version enables users to configure and run BLAST from a graphical interface, while the CLI version allows command-line execution.

### Author:
**SERRALTA Theo**

### Date:
**11/11/2024**

### Dependencies:
- Python 3.x
- Biopython 1.78 (e.g., `pip install biopython==1.78`)
- Numpy (e.g., `pip install numpy`)
- Tkinter (for GUI support, included in standard Python libraries)

For reproducibility, I recommend using the provided Conda environment available on the project's GIT.

## Download this repository

```bash
git clone https://github.com/theo-serralta/BLAST-Python-Implementation.git
cd BLAST-Python-Implementation
```

## Install dependencies

### Conda environment

Install [conda](https://docs.conda.io/en/latest/miniconda.html).

Create the conda environment and install dependencies:

```bash
conda env create -f blast_python.yml 
```

Load the conda environment:

```bash
conda activate blast_python
```

## Usage:

### You can choose to run the script in two modes: **GUI** or **CLI**.

To run the program in GUI mode:
```bash
python blast_python.py GUI
```

To run the program in CLI mode with the custom BLASTP algorithm:
```bash
python blast_python.py CLI -m blastp -q query.fasta -d database.fasta -o output.txt
```

To run BLASTN (using Biopython's BLAST+ integration):
```bash
python blast_python.py CLI -m blastn -q query.fasta -d database.fasta -o output.txt
```

To run BLASTX (using Biopython's BLAST+ integration):
```bash
python blast_python.py CLI -m blastx -q query.fasta -d database.fasta -o output.txt
```

## Arguments (for CLI mode):

- `-m, --method`: The BLAST method to use (`blastp`, `blastn`, `blastx`).
- `-q, --query`: Path to the query sequence file.
- `-d, --database`: Path to the database sequence file.
- `-o, --output`: Path to the output file.
- `-e, --evalue`: E-value threshold (default: 0.001).
- `-f, --outfmt`: Output format (default: 6).
- `-k, --k`: K-mer length for the custom BLASTP algorithm.
- `-w, --max_distance`: Maximum distance for double-hit detection in the custom BLASTP.
- `-b, --bit_score_threshold`: Bit score threshold for the custom BLASTP.
- `-t, --threads`: Number of threads for parallel processing.

## Testing

A `test_blast_python.py` file is available in the `tests/` directory to validate functionality. Run it using:

```bash
python -m unittest discover -s tests
```

## License:
MIT license

"""
Filename: blast_python.py
Description: This script implements a custom BLASTP algorithm along with options to run BLASTX and BLASTN
             using Biopython's integration with the BLAST+ command line tools. The custom BLASTP algorithm 
             includes k-mer extraction, double-hit detection, and gapped alignment extension. The GUI 
             version enables users to configure and run BLAST from a graphical interface, while the CLI 
             version allows command-line execution.

Author: SERRALTA Theo
Date: 11/11/2024

Dependencies:
    - Python 3.x
    - Biopython 1.78 (e.g., `pip install biopython==1.78`)
    - Numpy (e.g., `pip install numpy`)
    - Tkinter (for GUI support, included in standard Python libraries)

    I recommend using my Conda environment available on the project's GIT.

Usage:
    You can choose to run the script in two modes: GUI or CLI.

    To run the program in GUI mode:
        $ python blast_python.py GUI

    To run the program in CLI mode with a custom BLASTP algorithm:
        $ python blast_python.py CLI -m blastp -q query.fasta -d database.fasta -o output.txt
    
    To run BLASTN (using Biopython's BLAST+ integration):
        $ python blast_python.py CLI -m blastn -q query.fasta -d database.fasta -o output.txt
    
    To run BLASTX (using Biopython's BLAST+ integration):
        $ python blast_python.py CLI -m blastx -q query.fasta -d database.fasta -o output.txt
    
    Arguments (for CLI mode):
        -m, --method         : The BLAST method to use (blastp, blastn, blastx).
        -q, --query          : Path to the query sequence file.
        -d, --database       : Path to the database sequence file.
        -o, --output         : Path to the output file.
        -e, --evalue         : E-value threshold (default: 0.001).
        -f, --outfmt         : Output format (default: 6).
        -k, --k              : K-mer length for the custom BLASTP algorithm.
        -w, --max_distance   : Maximum distance for double-hit detection in the custom BLASTP.
        -b, --bit_score_threshold : Bit score threshold for the custom BLASTP.
        -t, --threads        : Number of threads for parallel processing.

License:
    MIT license
"""

import argparse
import os
import math
import sys
import time
import threading
import tkinter as tk
from tkinter import filedialog, simpledialog, ttk
import concurrent.futures
import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastxCommandline
from Bio.Blast import NCBIXML


def parse_arguments():
    """
    Parse command-line arguments for the CLI mode.

    Returns:
        Namespace: A namespace containing the parsed arguments as attributes.
        The attributes include:
        - method (str): The BLAST method to use (blastp, blastn, blastx).
        - query (str): The path to the query sequence file.
        - database (str): The path to the database sequence file.
        - output (str): The path to save the output file.
        - e_value_threshold (float): The E-value threshold for filtering results.
        - outfmt (int): The output format (default is 6).
        - k (int): K-mer length for the custom BLASTP algorithm (default is 3).
        - max_distance (int): Maximum distance for double-hit detection (default is 10).
        - bit_score_threshold (float): The bit score threshold for filtering alignments (default is 22).
        - threads (int): The number of threads for parallel processing (default is 4).
    """
    parser = argparse.ArgumentParser(description="Run BLAST algorithm in CLI mode.")

    # Add the required CLI arguments
    parser.add_argument("-m", "--method", required=True, help="The BLAST method to use (blastp, blastn, blastx)")
    parser.add_argument("-q", "--query", required=True, help="Path to the query sequence file.")
    parser.add_argument("-d", "--database", required=True, help="Path to the database sequence file.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output file.")
    parser.add_argument("-e", "--e_value_threshold", type=float, default=0.001, help="E-value threshold for filtering results.")
    parser.add_argument("-f", "--outfmt", type=int, default=6, help="Output format (default: 6).")
    parser.add_argument("-k", "--k", type=int, default=3, help="K-mer length for the custom BLASTP.")
    parser.add_argument("-w", "--max_distance", type=int, default=10, help="Maximum distance for double-hit detection.")
    parser.add_argument("-b", "--bit_score_threshold", type=float, default=22, help="Bit score threshold.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for parallel processing.")

    # Parse the arguments from the command line
    return parser.parse_args()

#Function to load query sequence fasta
def load_query_sequence(fasta_file):
    """
    Load the query sequence from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        str: The query sequence.
    """
    # Parse the FASTA file and retrieve the first sequence entry
    record = next(SeqIO.parse(fasta_file, "fasta"))
    
    # Convert the sequence into a string and return it
    return str(record.seq)

# Function to load sequences from a FASTA database
def load_fasta_database(fasta_file):
    """
    Load sequences from a FASTA file to use as a database.

    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        list: List of sequences from the file.
    """
    sequences = []  # Initialize an empty list to store sequences
    # Parse the FASTA file using SeqIO and extract each sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))  # Convert the sequence to a string and append to the list
    return sequences  # Return the list of sequences

# Function to obtain the BLOSUM62 score for a pair of amino acids
def get_blosum62_score(a, b, blosum_matrix):
    """
    Get the BLOSUM62 score for a pair of amino acids.

    Args:
        a (str): First amino acid.
        b (str): Second amino acid.
        blosum_matrix (dict): The BLOSUM62 substitution matrix.
    
    Returns:
        int: The BLOSUM62 score of the amino acid pair.
    """
    # Check if the exact pair (a, b) exists in the matrix
    if (a, b) in blosum_matrix:
        return blosum_matrix[(a, b)]
    # If the reverse pair (b, a) exists (since the matrix is symmetric), return its score
    elif (b, a) in blosum_matrix:
        return blosum_matrix[(b, a)]
    # If neither pair exists, return the default penalty for a mismatch or unknown substitution
    else:
        return -4  # default penalty for mismatch/gap

# Function to extract k-mers from a sequence
def extract_kmers(sequence, k):
    """
    Extract k-mers from a sequence.

    Args:
        sequence (str): The sequence to process.
        k (int): The length of the k-mers.
    
    Returns:
        list: A list of tuples containing k-mers and their positions.
    """
    kmers = []  # Initialize an empty list to store the k-mers and their positions

    # Loop through the sequence to extract k-mers of length 'k'
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]  # Extract the k-mer from position i
        kmers.append((kmer, i))  # Store the k-mer along with its starting index
    return kmers  # Return the list of k-mers and their positions

# Function to index k-mers in a target sequence for fast searching
def index_target_sequence(sequence, k):
    """
    Find positions of k-mers in the target sequence using an index.

    Args:
        kmers (list): List of k-mers from the query sequence.
        target_index (dict): Precomputed index of the target sequence.
    
    Returns:
        dict: Dictionary of k-mer positions in both query and target sequences.
    """
    kmer_index = {}  # Initialize an empty dictionary to store k-mers and their positions
    # Iterate through the target sequence to extract and index k-mers of length 'k'
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]  # Extract the k-mer starting at position i
        if kmer not in kmer_index:
            kmer_index[kmer] = []  # Create a new entry for the k-mer if it doesn't exist
        kmer_index[kmer].append(i)  # Append the position to the list of positions for this k-mer
    return kmer_index  # Return the dictionary containing k-mers and their positions

# Function to find positions of k-mers in the query and target sequences
def find_kmer_positions(kmers, target_index):
    """
    Finds the positions of the extracted k-mers in an indexed target sequence.

    Args:
        kmers (list): List of k-mers with their positions.
        target_index (dict): Index of the target sequence.

    Returns:
        dict: Dictionary of k-mers and their respective positions in the query and target sequences.
    """
    positions = {}  # Initialize a dictionary to store matching k-mer positions
    # Iterate through each k-mer from the query along with its position
    for kmer, query_pos in kmers:
        # Check if the k-mer exists in the target sequence index
        if kmer in target_index:
            # For each matching k-mer in the target, store its position in both query and target
            for target_pos in target_index[kmer]:
                if kmer not in positions:
                    positions[kmer] = []  # Initialize a list for new k-mers
                positions[kmer].append((target_pos, query_pos))  # Store target and query positions
    return positions  # Return the dictionary of matching positions

# Function to detect double hits in the sequence comparison
# Double hits are pairs of k-mers that align on the same diagonal without overlap
def find_double_hits(kmer_positions, max_distance):
    """
    Identify double-hits between the query and target sequences.

    Args:
        kmer_positions (dict): The positions of k-mers in the target sequence.
        max_distance (int): Maximum allowable distance between two k-mers for a double-hit.
    
    Returns:
        list: A list of double-hit positions.
    """
    double_hits = []  # Initialize a list to store double hits
    kmer_list = list(kmer_positions.items())  # Convert the dictionary to a list for indexed access
    
    # Loop over all k-mers in the query and target sequences
    for i in range(len(kmer_list)):
        kmer1, positions1 = kmer_list[i]  # Get the first k-mer and its positions
        # Compare with other k-mers in the list
        for pos1, query_pos1 in positions1:
            for j in range(i + 1, len(kmer_list)):
                kmer2, positions2 = kmer_list[j]  # Get the second k-mer and its positions
                for pos2, query_pos2 in positions2:
                    # Check if the k-mers are on the same diagonal
                    if (pos1 - query_pos1) == (pos2 - query_pos2):
                        # Ensure the two hits are close enough but not overlapping
                        if abs(pos2 - pos1) <= max_distance:
                            # Ensure they are not in exactly the same location (strict non-overlap)
                            if pos1 != pos2 or query_pos1 != query_pos2:
                                # Double hit found, append details
                                double_hits.append((kmer1, pos1, query_pos1, kmer2, pos2, query_pos2))
    return double_hits  # Return the list of detected double hits

# Function to calculate the score of a double-hit between two k-mers
def evaluate_double_hit(kmer1, kmer2, blosum_matrix):
    """
    Calculate the score of a double-hit using the BLOSUM62 matrix.

    Args:
        kmer1 (str): First k-mer.
        kmer2 (str): Second k-mer.
        blosum_matrix (dict): The BLOSUM62 matrix.
    
    Returns:
        int: The score of the double-hit.
    """
    score = 0  # Initialize the score to zero
    # Iterate over the amino acids of the two k-mers (assuming both are of the same length)
    for a, b in zip(kmer1, kmer2):
        # Add the BLOSUM62 score for the amino acid pair (a, b)
        score += get_blosum62_score(a, b, blosum_matrix)
    return score  # Return the total score for the double-hit

# Extend a double-hit with gapped alignment
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2):
    """
    Extend a double-hit into a full alignment using Smith-Waterman local alignment.

    Args:
        seq_query (str): The query sequence.
        seq_target (str): The target sequence.
        query_pos (int): Position of the query hit.
        target_pos (int): Position of the target hit.
        blosum_matrix (dict): The BLOSUM62 substitution matrix.
        gap_open_penalty (int): Penalty for opening a gap in the alignment.
        gap_extension_penalty (int): Penalty for extending a gap.
    
    Returns:
        tuple: Alignment score and the final alignment.
    """
    m, n = len(seq_query), len(seq_target)  # Lengths of the query and target sequences
    
    # Initialize the dynamic programming score matrix (Smith-Waterman) and traceback matrix
    previous_row = np.zeros(n + 1)  # Previous row in the dynamic programming matrix
    current_row = np.zeros(n + 1)   # Current row in the dynamic programming matrix
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)  # Traceback matrix for tracking alignment path
    
    # Variables to track the best score and its position
    max_score = 0
    max_pos = (0, 0)  # Position of the maximum score in the alignment matrix

    # Iterate through each position in the query (i) and target (j) sequences
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate match score using the BLOSUM62 matrix for amino acid substitution
            try:
                match = previous_row[j - 1] + get_blosum62_score(seq_query[i - 1], seq_target[j - 1], blosum_matrix)
            except KeyError:
                match = -4  # Default mismatch penalty for unknown amino acid pairs

            # Calculate gap penalties (left: gap in query, up: gap in target)
            gap_query = current_row[j - 1] + (gap_extension_penalty if traceback_matrix[i][j - 1] == 2 else gap_open_penalty)
            gap_target = previous_row[j] + (gap_extension_penalty if traceback_matrix[i - 1][j] == 1 else gap_open_penalty)

            # Choose the best score for the current cell (0 for local alignment termination, match, or gap)
            current_score = max(0, match, gap_query, gap_target)
            current_row[j] = current_score  # Update the current row with the computed score

            # Update the traceback matrix based on the direction of the best score
            if current_score == match:
                traceback_matrix[i][j] = 0  # Diagonal move (match/mismatch)
            elif current_score == gap_query:
                traceback_matrix[i][j] = 2  # Left move (gap in query)
            elif current_score == gap_target:
                traceback_matrix[i][j] = 1  # Up move (gap in target)

            # Track the highest score and its position in the matrix
            if current_score > max_score:
                max_score = current_score
                max_pos = (i, j)

        # Swap rows for the next iteration (use current_row for the next previous_row)
        previous_row, current_row = current_row, np.zeros(n + 1)  # Reset the current row for the next iteration

    # Perform traceback to recover the full alignment
    aligned_query, aligned_target = [], []
    i, j = max_pos  # Start traceback from the position of the maximum score

    while i > 0 and j > 0 and traceback_matrix[i][j] != -1:  # Continue until reaching the edge of the matrix
        if traceback_matrix[i][j] == 0:
            aligned_query.append(seq_query[i - 1])
            aligned_target.append(seq_target[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:
            aligned_query.append(seq_query[i - 1])
            aligned_target.append('-')  # Gap in the target sequence
            i -= 1
        elif traceback_matrix[i][j] == 2:
            aligned_query.append('-')  # Gap in the query sequence
            aligned_target.append(seq_target[j - 1])
            j -= 1

    # Reverse the aligned sequences to return them in the correct order
    aligned_query = aligned_query[::-1]
    aligned_target = aligned_target[::-1]
    final_alignment = list(zip(aligned_query, aligned_target))  # Combine query and target alignments into pairs

    return max_score, final_alignment  # Return the maximum score and the final alignment

# Filter and extend alignments based on double-hits
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, bit_score_threshold):
    """
    Filter and extend alignments based on double-hits.

    Args:
        double_hits (list): The detected double-hits.
        seq_query (str): The query sequence.
        seq_target (str): The target sequence.
        blosum_matrix (dict): The BLOSUM62 matrix.
        bit_score_threshold (float): Bit score threshold for filtering alignments.
    
    Returns:
        list: List of filtered alignments.
    """
    alignments = []

    # Return an empty list if no double-hits were detected
    if not double_hits:
        return alignments

    # Initialize variables to track the best double-hit with the largest distance between k-mers
    best_hit = None
    max_distance = 0

    # Iterate over all detected double-hits to find the one with the largest distance between k-mers
    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        distance = abs(pos2 - pos1)  # Calculate the distance between the two k-mers in the target sequence

        # Update the best_hit if this double-hit has a greater distance
        if distance > max_distance:
            max_distance = distance
            best_hit = hit

    # If a best double-hit is found, extend the alignment
    if best_hit:
        # Extract positions and k-mers from the best double-hit
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = best_hit

        # Extend the alignment starting from the double-hit using Smith-Waterman
        score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)

        # Calculate the bit score of the extended alignment
        bit_score = calculate_bit_score(score)

        # Only keep the alignment if the bit score meets or exceeds the specified threshold
        if bit_score >= bit_score_threshold:
            alignments.append((score, alignment))

    return alignments

#Function to filter duplicate alignments
def filter_duplicate_alignments(alignments):
    """
    Remove duplicate alignments based on positions.

    Args:
        alignments (list): List of alignments.
    
    Returns:
        list: Unique alignments with duplicates removed.
    """
    unique_alignments = []  # Store unique alignments
    seen_positions = set()  # Set to track alignment positions already encountered

    # Iterate through the list of alignments
    for score, alignment in alignments:
        # Create a hashable key from the alignment based on the positions
        alignment_key = tuple(alignment)

        # Only add alignment if its position pattern has not been seen before
        if alignment_key not in seen_positions:
            unique_alignments.append((score, alignment))  # Append unique alignment
            seen_positions.add(alignment_key)  # Mark position as seen

    return unique_alignments

# Function to calculate the bit score from the raw alignment score
def calculate_bit_score(raw_score, lambda_param=0.318, k=0.134):
    """
    Calculate the bit score from a raw score.

    Args:
        raw_score (int): Raw alignment score.
        lambda_param (float): Lambda parameter for scoring.
        K (float): K parameter for scoring.
    
    Returns:
        float: The calculated bit score.
    """
    # Convert the raw score into a normalized bit score using the formula
    return (lambda_param * raw_score - math.log(k)) / math.log(2)

# Function to calculate the E-value from a bit score
def calculate_e_value_from_bitscore(bit_score, m, n):
    """
    Calculate the E-value from a bit score.

    Args:
        bit_score (float): The bit score.
        m (int): Length of the query sequence.
        n (int): Total length of the database.
    
    Returns:
        float: The calculated E-value.
    """
    # Calculate the E-value using the standard formula: E = m * n * 2^(-bit_score)
    return m * n * math.pow(2, -bit_score)

# Calculate E-values for all alignments
def calculate_e_values(alignments, seq_query, len_database):
    """
    Calculate E-values for all alignments.

    Args:
        alignments (list): List of alignments with their scores.
        seq_query (str): The query sequence.
        len_database (int): Total length of the database.
    
    Returns:
        list: List of tuples containing E-values, scores, and alignments.
    """
    e_values = []  # Store the resulting E-values for all alignments
    m = len(seq_query)  # Length of the query sequence
    n = len_database  # Total length of the database

    # Loop through each alignment and its corresponding raw score
    for raw_score, alignment in alignments:
        bit_score = calculate_bit_score(raw_score)  # Calculate the bit score
        e_value = calculate_e_value_from_bitscore(bit_score, m, n)  # Compute E-value using the bit score
        e_values.append((e_value, raw_score, alignment))  # Append the result as a tuple

    return e_values  # Return the list of E-values, scores, and alignments

# Filter alignments by E-value
def filter_by_e_value(e_values, e_value_threshold):
    """
    Filter alignments based on an E-value threshold.

    Args:
        e_values (list): List of alignments with their E-values.
        e_value_threshold (float): E-value threshold.
    
    Returns:
        list: Filtered alignments with E-value below the threshold.
    """
    # Use list comprehension to filter alignments based on E-value threshold
    filtered_alignments = [align for align in e_values if align[0] <= e_value_threshold]
    
    return filtered_alignments  # Return only alignments with significant E-values

# Function to format the alignment for BLAST-like display
def format_alignment(seq_query, seq_target, alignment):
    """
    Format an alignment for BLAST-like display.

    Args:
        seq_query (str): The query sequence.
        seq_target (str): The target sequence.
        alignment (list): The alignment.
    
    Returns:
        tuple: Formatted query, match line, and target sequence strings.
    """
    # Lists to store the formatted query sequence, match line, and target sequence
    query_aligned = []
    target_aligned = []
    match_line = []

    # Iterate over the aligned residues in the query and target sequences
    for q_res, t_res in alignment:
        # Append the query residue or gap
        query_aligned.append(q_res if q_res != '-' else '-')
        # Append the target residue or gap
        target_aligned.append(t_res if t_res != '-' else '-')

        # Determine the type of match and add the appropriate symbol
        if q_res == t_res:
            match_line.append('|')  # Perfect match between query and target
        elif get_blosum62_score(q_res, t_res, MatrixInfo.blosum62) > 0:
            match_line.append(':')  # Reasonable substitution based on BLOSUM62 score
        else:
            match_line.append('.')  # Mismatch or weak substitution

    # Convert the aligned sequences and match line to strings
    query_str = "".join(query_aligned)
    match_str = "".join(match_line)
    target_str = "".join(target_aligned)

    return query_str, match_str, target_str  # Return the formatted alignment

# Display formatted alignment results
def display_blast_like_results(alignments, e_values, seq_query, seq_target, output_file):
    """
    Display alignment results in a BLAST-like format.

    Args:
        alignments (list): List of alignments.
        e_values (list): List of E-values.
        seq_query (str): The query sequence.
        seq_target (str): The target sequence.
        output_file (str): Path to the output file.
    """
    # Parameters for formatting
    line_length = 60  # Number of characters per line in the alignment display

    with open(output_file, "a") as out_f:
        out_f.write(f"{len(e_values)} significant alignments found:\n\n")
        
        for i, (e_value, raw_score, alignment) in enumerate(e_values, 1):
            bit_score = calculate_bit_score(raw_score)
            query_str, match_str, target_str = format_alignment(seq_query, seq_target, alignment)

            # Calculate identities, positives, and gaps
            identities = sum(1 for q, t in alignment if q == t)
            positives = sum(1 for q, t in alignment if get_blosum62_score(q, t, MatrixInfo.blosum62) > 0)
            gaps = sum(1 for q, t in alignment if q == '-' or t == '-')

            # Calculate the percentages
            total_length = len(alignment)
            identity_percentage = (identities / total_length) * 100
            positive_percentage = (positives / total_length) * 100
            gap_percentage = (gaps / total_length) * 100

            # Write summary information
            out_f.write(f"Alignment {i}:\n")
            out_f.write(f" Score = {raw_score}, bits=({bit_score:.2f}), Expect = {e_value:.2e}\n")
            out_f.write(f" Identities = {identities}/{total_length} ({identity_percentage:.0f}%), "
                        f"Positives = {positives}/{total_length} ({positive_percentage:.0f}%), "
                        f"Gaps = {gaps}/{total_length} ({gap_percentage:.0f}%)\n\n")

            # Write alignment block by block
            for block_start in range(0, total_length, line_length):
                block_end = block_start + line_length
                query_block = query_str[block_start:block_end]
                match_block = match_str[block_start:block_end]
                target_block = target_str[block_start:block_end]

                query_label = f"Query {block_start + 1}".ljust(9)
                target_label = f"Sbjct {block_start + 1}".ljust(9)
                
                out_f.write(f"{query_label} {query_block}\n")
                out_f.write(f"{''.ljust(9)} {match_block}\n")
                out_f.write(f"{target_label} {target_block}\n\n")

            out_f.write("\n" + "-"*50 + "\n\n")

def sort_output_file(output_file):
    """
    Sort the output file by E-value and overwrite with sorted results.

    Args:
        output_file (str): The file to be sorted.
    """
    # Open and read the output file
    with open(output_file, "r") as f:
        lines = f.readlines()

    # Placeholder for storing the alignment blocks
    alignments = []
    current_alignment = []

    # Parse the file line by line and collect alignments
    for line in lines:
        if line.startswith("Alignment"):  # New alignment block found
            if current_alignment:  # If an alignment is already being processed, save it
                alignments.append(current_alignment)
            current_alignment = [line]  # Start a new alignment block
        else:
            current_alignment.append(line)  # Continue collecting lines for the current alignment
    
    if current_alignment:  # Ensure the last alignment is added
        alignments.append(current_alignment)

    # Define a helper function to extract E-value from an alignment block
    def extract_evalue(alignment):
        """
        Extracts the E-value from an alignment block.

        Args:
            alignment (list): A list of lines representing an alignment block.

        Returns:
            float: The E-value extracted from the alignment, or infinity if no E-value is found.
        """
        for line in alignment:
            if "Expect" in line:
                try:
                    # Extract the E-value by splitting the line and converting to float
                    evalue_str = line.split('=')[-1].strip()
                    return float(evalue_str.replace("e", "E"))  # Convert scientific notation
                except ValueError:
                    return float('inf')  # Return infinity if the E-value cannot be parsed
        return float('inf')  # Return infinity if no E-value is found

    # Sort the alignment blocks by their extracted E-values (smallest first)
    alignments.sort(key=extract_evalue)

    # Write the sorted alignments back to the file
    with open(output_file, "w") as f:
        for alignment in alignments:
            f.writelines(alignment)  # Write each sorted alignment block

def run_blastx(query_fasta, protein_db, output_file, e_value_threshold=0.001, outfmt=6):
    """
    Run BLASTX using the query sequence against a protein database.

    Args:
        query_fasta (str): Path to the query FASTA file.
        protein_db (str): Path to the protein database.
        output_file (str): Path to save BLASTX results.
        e_value_threshold (float): E-value threshold for filtering results.
        outfmt (int): Output format.
    """
    # Create BLASTX command with the specified parameters for query, database, E-value, output format, and file
    blastx_cline = NcbiblastxCommandline(query=query_fasta, db=protein_db, evalue=e_value_threshold, outfmt=outfmt, out=output_file)
    
    # Display the command being run (useful for debugging or tracking the exact command executed)
    print(f"Running BLASTX command: {blastx_cline}")
    
    # Execute the BLASTX command and capture any standard output or errors
    stdout, stderr = blastx_cline()

def run_blastn(query_fasta, nucleotide_db, output_file, e_value_threshold=0.001, outfmt=6):
    """
    Run BLASTN using the query sequence against a nucleotide database.

    Args:
        query_fasta (str): Path to the query FASTA file.
        nucleotide_db (str): Path to the nucleotide database.
        output_file (str): Path to save BLASTN results.
        e_value_threshold (float): E-value threshold for filtering results.
        outfmt (int): Output format.
    """
    # Create BLASTN command with specified parameters: query file, database, E-value threshold, output format, and output file
    blastn_cline = NcbiblastnCommandline(query=query_fasta, db=nucleotide_db, evalue=e_value_threshold, outfmt=outfmt, out=output_file)
    
    # Print the constructed BLASTN command to confirm what is being executed (useful for debugging and transparency)
    print(f"Running BLASTN command: {blastn_cline}")
    
    # Execute the BLASTN command and capture any standard output or errors generated during the process
    stdout, stderr = blastn_cline()

def parse_blast_results(blast_xml_file):
    """
    Parse BLAST XML results and print alignments.

    Args:
        blast_xml_file (str): Path to BLAST XML file.
    """
    # Open the XML file containing the BLAST results for reading
    with open(blast_xml_file) as result_handle:
        
        # Parse the BLAST results from the XML format using Biopython's NCBIXML parser
        blast_records = NCBIXML.parse(result_handle)
        
        # Iterate over each BLAST record in the results
        for blast_record in blast_records:
            
            # Iterate through each alignment in the current BLAST record
            for alignment in blast_record.alignments:
                
                # Iterate through the high-scoring pairs (HSPs) in each alignment
                for hsp in alignment.hsps:
                    
                    # Print key information about the alignment
                    print(f"****Alignment****")
                    print(f"sequence: {alignment.title}")  # Title of the sequence aligned
                    print(f"length: {alignment.length}")   # Length of the aligned sequence
                    print(f"e-value: {hsp.expect}")       # E-value of the HSP (lower is better)
                    
                    # Print the first 75 characters of the query, match, and subject sequences for clarity
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...\n")

# Function to run alignments in parallel using multiple threads
def run_parallel_alignments(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold, num_threads, e_value_threshold, output_file):
    """
    Run alignments in parallel across multiple threads.

    Args:
        seq_query (str): Query sequence.
        database_sequences (list): List of database sequences.
        k (int): K-mer length.
        max_distance (int): Maximum distance for double-hits.
        blosum62 (dict): BLOSUM62 matrix.
        bit_score_threshold (float): Bit score threshold.
        num_threads (int): Number of threads for parallel processing.
        e_value_threshold (float): E-value threshold.
        output_file (str): Output file to save results.
    """
    # Step 1: Calculate the total length of the database (sum of all target sequence lengths)
    len_database = sum(len(seq) for seq in database_sequences)

    # Step 2: Initialize a process pool for parallel execution, limiting to the number of threads specified
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:

        # Step 3: Submit alignment tasks to the pool
        # Each target sequence is aligned to the query in parallel
        future_to_sequence = {
            executor.submit(align_sequence, seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database, e_value_threshold): seq_target
            for seq_target in database_sequences
        }

        # Step 4: Process the results as they complete
        for future in concurrent.futures.as_completed(future_to_sequence):
            seq_target = future_to_sequence[future]
            try:
                # Retrieve the result from the completed task
                result = future.result()

                # If valid alignments are found, display and save the results
                if result:
                    alignments, e_values = result
                    display_blast_like_results(alignments, e_values, seq_query, seq_target, output_file)
            except Exception as exc:
                # If an exception occurs during alignment, log the error
                print(f"Generated an exception for {seq_target}: {exc}")

# Function to align a target sequence to the query sequence
def align_sequence(seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database, e_value_threshold):
    """
    Perform alignment of a target sequence to the query sequence.

    Args:
        seq_target (str): Target sequence.
        seq_query (str): Query sequence.
        k (int): K-mer length.
        max_distance (int): Maximum allowed distance for double-hits.
        blosum62 (dict): BLOSUM62 matrix.
        bit_score_threshold (float): Bit score threshold.
        len_database (int): Length of the database.
        output_file (str): Path to save the alignment output.
    
    Returns:
        tuple: Alignments and E-values if any alignments are found.
    """
    # Step 1: Extract k-mers from the query sequence
    kmers = extract_kmers(seq_query, k)

    # Step 2: Index the target sequence by k-mers
    target_index = index_target_sequence(seq_target, k)

    # Step 3: Find positions of matching k-mers between query and target
    kmer_positions = find_kmer_positions(kmers, target_index)

    # Step 4: Detect double-hits (pairs of matching k-mers on the same diagonal)
    double_hits = find_double_hits(kmer_positions, max_distance)

    # Step 5: Filter and extend alignments based on detected double-hits
    alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, bit_score_threshold)

    # Step 6: Remove any duplicate alignments
    alignments = filter_duplicate_alignments(alignments)

    # Step 7: Calculate E-values for the alignments
    e_values = calculate_e_values(alignments, seq_query, len_database)

    # Step 8: Filter the alignments based on the E-value threshold
    significant_alignments = filter_by_e_value(e_values, e_value_threshold)

    # Return the significant alignments and E-values if found, otherwise return None
    if significant_alignments:
        return alignments, e_values
    return None


def browse_file(label):
    """
    Open a file dialog to select a file and display the file path.

    Args:
        label (tk.Label): Label widget to update with the selected file path.
    """
    filename = filedialog.askopenfilename()  # Open file dialog
    label.config(text=filename)  # Update the label with the selected file path
    return filename

def ask_db_name_blastx(label):
    """
    Ask the user to input the database name for BLASTX.

    Args:
        label (tk.Label): Label widget to update with the entered database name.
    """
    db_name = simpledialog.askstring("Database Name", "Enter the database name (e.g., 'exemple_db_blastx'):")
    if db_name:
        label.config(text=db_name)
    return db_name

def ask_db_name_blastn(label):
    """
    Ask the user to input the database name for BLASTN.

    Args:
        label (tk.Label): Label widget to update with the entered database name.
    """
    db_name = simpledialog.askstring("Database Name", "Enter the database name (e.g., 'exemple_db_blastn'):")
    if db_name:
        label.config(text=db_name)
    return db_name

def run_blast_thread():
    """
    Launch BLAST execution in a separate thread to keep the GUI responsive.

    This function resets the progress bar, starts it, and spawns a new thread
    to run the `run_blast()` function without blocking the main GUI thread.
    """
    progress_bar['value'] = 0  # Reset progress bar to 0%
    progress_bar.start(10)  # Start progress bar animation
    threading.Thread(target=run_blast).start()  # Start BLAST execution in a separate thread

def run_blast():
    """
    Execute the selected BLAST method based on user input from the GUI.

    This function retrieves the user inputs from the GUI, such as the query and database files,
    and the selected BLAST method, then runs the appropriate BLAST function. It updates the
    progress bar during the execution.
    """
    method = method_var.get()  # Get the selected BLAST method
    query = query_file_label.cget("text")  # Get the query file path
    database = db_file_label.cget("text")  # Get the database file path or name
    output = output_entry.get()  # Get the output file path
    e_value_threshold = float(evalue_entry.get())  # Get the E-value threshold
    outfmt = int(outfmt_entry.get())  # Get the output format
    k = int(k_entry.get())  # Get the k-mer length
    max_distance = int(max_dist_entry.get())  # Get the max distance for double-hit detection
    bit_score_threshold = float(bit_score_entry.get())  # Get the bit score threshold
    num_threads = int(threads_entry.get())  # Get the number of threads

    start = time.time()  # Start the timer

    # Depending on the selected method, run the appropriate BLAST process
    if method == "blastp":
        seq_query = load_query_sequence(query)
        database_sequences = load_fasta_database(database)
        blosum62 = MatrixInfo.blosum62  # Load BLOSUM62 matrix for scoring
        run_parallel_alignments(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold, num_threads, e_value_threshold, output)
        if os.path.exists(output):
            sort_output_file(output)  # Sort the output file if alignments were found
        else:
            print(f"No alignments found. Output file '{output}' not created.")
    elif method == "blastx":
        run_blastx(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)
    elif method == "blastn":
        run_blastn(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)

    end = time.time()  # End the timer
    print(f"Execution time for {method} = {end - start} seconds")  # Print the execution time

    # Update progress bar
    progress_bar.stop()  # Stop the progress bar animation
    progress_bar['value'] = 100  # Set progress bar to 100%

def start_gui():
    """
    Initialize and display the GUI for BLAST tool configuration.

    This function sets up the graphical interface for configuring BLAST options, 
    selecting query and database files, and running the BLAST process.
    """
    global method_var, query_file_label, db_file_label, output_entry, evalue_entry, outfmt_entry
    global k_entry, max_dist_entry, bit_score_entry, threads_entry, progress_bar

    # Create the main window
    root = tk.Tk()
    root.title("Custom BLAST GUI")

    # BLAST method selection
    tk.Label(root, text="Method:").grid(row=0, column=0)
    method_var = tk.StringVar(value="blastp")  # Default to "blastp"
    tk.OptionMenu(root, method_var, "blastp", "blastn", "blastx").grid(row=0, column=1)

    # Query file selection
    tk.Label(root, text="Query file:").grid(row=1, column=0)
    query_file_label = tk.Label(root, text="No file selected", width=40)
    query_file_label.grid(row=1, column=1)
    tk.Button(root, text="Browse", command=lambda: browse_file(query_file_label)).grid(row=1, column=2)

    # Database selection
    tk.Label(root, text="Database:").grid(row=2, column=0)
    db_file_label = tk.Label(root, text="No database selected", width=40)
    db_file_label.grid(row=2, column=1)
    db_button = tk.Button(root, text="Browse", command=lambda: browse_file(db_file_label))
    db_button.grid(row=2, column=2)

    # Update database selection UI based on method
    def update_db_selection(*args):
        method = method_var.get()
        if method == "blastp":
            db_button.config(text="Browse", command=lambda: browse_file(db_file_label))
        elif method == "blastx":
            db_button.config(text="Enter Name", command=lambda: ask_db_name_blastx(db_file_label))
        elif method == "blastn":
            db_button.config(text="Enter Name", command=lambda: ask_db_name_blastn(db_file_label))

    method_var.trace_add("write", update_db_selection)  # Bind method selection to update DB behavior

    # Output file selection
    tk.Label(root, text="Output file:").grid(row=3, column=0)
    output_entry = tk.Entry(root, width=40)
    output_entry.grid(row=3, column=1)

    # E-value threshold input
    tk.Label(root, text="E-value threshold:").grid(row=4, column=0)
    evalue_entry = tk.Entry(root, width=10)
    evalue_entry.grid(row=4, column=1)
    evalue_entry.insert(0, "0.001")

    # Output format input
    tk.Label(root, text="Output format:").grid(row=5, column=0)
    outfmt_entry = tk.Entry(root, width=10)
    outfmt_entry.grid(row=5, column=1)
    outfmt_entry.insert(0, "6")

    # K-mer length input
    tk.Label(root, text="K-mer length:").grid(row=6, column=0)
    k_entry = tk.Entry(root, width=10)
    k_entry.grid(row=6, column=1)
    k_entry.insert(0, "3")

    # Max distance for double-hit detection input
    tk.Label(root, text="Max distance:").grid(row=7, column=0)
    max_dist_entry = tk.Entry(root, width=10)
    max_dist_entry.grid(row=7, column=1)
    max_dist_entry.insert(0, "10")

    # Bit score threshold input
    tk.Label(root, text="Bit score threshold:").grid(row=8, column=0)
    bit_score_entry = tk.Entry(root, width=10)
    bit_score_entry.grid(row=8, column=1)
    bit_score_entry.insert(0, "22")

    # Number of threads input
    tk.Label(root, text="Number of threads:").grid(row=9, column=0)
    threads_entry = tk.Entry(root, width=10)
    threads_entry.grid(row=9, column=1)
    threads_entry.insert(0, "4")

    # Progress bar
    progress_bar = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
    progress_bar.grid(row=10, columnspan=3, pady=10)

    # Run button to trigger BLAST
    tk.Button(root, text="Run BLAST", command=run_blast_thread).grid(row=11, column=1)

    # Start the GUI event loop
    root.mainloop()

def run_cli():
    """
    Execute BLAST based on command-line input.

    This function retrieves the parsed command-line arguments, then runs
    the appropriate BLAST method (blastp, blastx, or blastn).
    """
    # Get parsed arguments
    args = parse_arguments()

    start = time.time()

    # Extract arguments from 'args'
    method = args.method
    query = args.query
    database = args.database
    output = args.output
    e_value_threshold = args.e_value_threshold
    outfmt = args.outfmt
    k = args.k
    max_distance = args.max_distance
    bit_score_threshold = args.bit_score_threshold
    num_threads = args.threads

    # Run the appropriate BLAST method
    if method == "blastp":
        seq_query = load_query_sequence(query)
        database_sequences = load_fasta_database(database)
        blosum62 = MatrixInfo.blosum62
        run_parallel_alignments(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold, num_threads, e_value_threshold, output)
        if os.path.exists(output):
            sort_output_file(output)
        else:
            print(f"No alignments found. Output file '{output}' not created.")
    elif method == "blastx":
        run_blastx(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)
    elif method == "blastn":
        run_blastn(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)
    else:
        print(f"Unknown method: {method}. Please choose from 'blastp', 'blastx', or 'blastn'.")

    end = time.time()
    print(f"Execution time for {method} = {end - start} seconds")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Check if the first argument is "GUI" or "CLI"
        mode = sys.argv[1].upper()

        if mode == "GUI":
            # Remove "GUI" from sys.argv and start the graphical interface
            sys.argv.pop(1)
            start_gui()
        elif mode == "CLI":
            # Remove "CLI" from sys.argv and parse the rest of the arguments for CLI
            sys.argv.pop(1)
            run_cli()  # Call the function to handle CLI logic
        else:
            print("Usage: python blast_python.py [GUI|CLI]")
    else:
        print("Usage: python blast_python.py [GUI|CLI]")

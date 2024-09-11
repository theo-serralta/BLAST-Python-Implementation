"""
Filename: blast_python.py
Description: This script implements a custom version of BLASTP, and utilizes Biopython 
             for running BLASTX and BLASTN. The custom BLASTP algorithm includes k-mer extraction, 
             diagonal double-hit detection, and gapped alignment extension. BLASTX and BLASTN 
             functionalities leverage Biopython's integration with the BLAST+ command line tools.

Author: SERRALTA Theo
Date: 11/11/2024

Dependencies:
    - Python 3.x
    - Biopython 1.78 (e.g., `pip install biopython`)
    - Numpy (e.g., `pip install numpy`)

    I recommend using my Conda environment available on the project's GIT.

Usage:
    To run the custom BLASTP algorithm:
        $ python blast_python.py -m blastp -q query.fasta -d database.fasta -o output.txt
    
    To run BLASTN with a query file and a nucleotide database (using Biopython):
        $ python blast_python.py -m blastn -q query.fasta -d database.fasta -o output.txt
    
    To run BLASTX with a query file and a protein database (using Biopython):
        $ python blast_python.py -m blastx -q query.fasta -d database.fasta -o output.txt
    
    Arguments:
        -m, --method         : The BLAST method to use (blastp, blastn, blastx).
        -q, --query          : Path to the query sequence file.
        -d, --database       : Path to the database sequence file.
        -o, --output         : Path to the output file.
        -e, --evalue         : E-value threshold (default: 0.001).
        -f, --outfmt         : Output format (default: 6).
        -k, --k              : K-mer length for the custom BLASTP.
        -w, --max_distance   : Maximum distance for double hits in the custom BLASTP.
        -b, --bit_score_threshold : Bit score threshold for the custom BLASTP.
        -t, --threads        : Number of threads for parallel processing.

License:
    MIT license
"""

import os
import math
import sys
import time
import getopt
import concurrent.futures
import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML


def parse_arguments():
    """
    Parse command-line arguments.

    This function parses the command-line arguments provided to the script and extracts the values for various parameters.

    Returns:
    --------
    tuple
        Returns a tuple containing:
        - method : str : The method to use (blastp, blastn, blastx)
        - query : str : The query sequence file path
        - database : str : The database file path
        - output : str : The output file path
        - e_value_threshold : float : The E-value threshold for filtering alignments
        - outfmt : int : The output format (default is 6)
        - k : int : K-mer size for custom blastp
        - max_distance : int : Maximum distance between k-mer hits for double-hit detection
        - bit_score_threshold : float : Threshold for filtering alignments based on the bit score
        - num_threads : int : Number of threads for parallel execution
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:q:d:o:e:f:k:w:b:t:", 
                                   ["method=", "query=", "database=", "output=", "e_value_threshold=", "outfmt=", 
                                    "k=", "max_distance=", "bit_score_threshold=", "threads="])
        
        # Initialisation des valeurs par défaut
        method = query = database = output = None
        e_value_threshold = 0.001
        outfmt = 6
        k = 3
        max_distance = 10
        bit_score_threshold = 22
        num_threads = 4

        for opt, arg in opts:
            if opt in ("-m", "--method"):
                method = arg
            elif opt in ("-q", "--query"):
                query = arg
            elif opt in ("-d", "--database"):
                database = arg
            elif opt in ("-o", "--output"):
                output = arg
            elif opt in ("-e", "--e_value_threshold"):
                e_value_threshold = float(arg)
            elif opt in ("-f", "--outfmt"):
                outfmt = int(arg)
            elif opt in ("-k", "--k"):
                k = int(arg)
            elif opt in ("-w", "--max_distance"):
                max_distance = int(arg)
            elif opt in ("-b", "--bit_score_threshold"):
                bit_score_threshold = float(arg)
            elif opt in ("-t", "--threads"):
                num_threads = int(arg)

        if not method or not query or not database or not output:
            print("Missing required arguments. Please specify method, query, database, and output.")
            sys.exit(1)

        return (method, query, database, output, e_value_threshold, outfmt, k, max_distance, bit_score_threshold, num_threads)

    except getopt.GetoptError as err:
        print(f"Error parsing arguments: {err}")
        sys.exit(1)


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
    if (a, b) in blosum_matrix:
        return blosum_matrix[(a, b)]
    elif (b, a) in blosum_matrix:
        return blosum_matrix[(b, a)]
    else:
        return -4  # default penalty for mismatch/gap

# K-mer extraction
def extract_kmers(sequence, k):
    """
    Extract k-mers from a sequence.

    Args:
        sequence (str): The sequence to process.
        k (int): The length of the k-mers.
    
    Returns:
        list: A list of tuples containing k-mers and their positions.
    """
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append((kmer, i))  # Preserve the original index
    return kmers

# Indexing the target sequence for fast k-mer search
def index_target_sequence(sequence, k):
    """
    Find positions of k-mers in the target sequence using an index.

    Args:
        kmers (list): List of k-mers from the query sequence.
        target_index (dict): Precomputed index of the target sequence.
    
    Returns:
        dict: Dictionary of k-mer positions in both query and target sequences.
    """
    kmer_index = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmer_index:
            kmer_index[kmer] = []
        kmer_index[kmer].append(i)  # Preserve the original index
    return kmer_index

# Find k-mer positions in the target sequence
def find_kmer_positions(kmers, target_index):
    """
    Finds the positions of the extracted k-mers in an indexed target sequence.

    Args:
        kmers (list): List of k-mers with their positions.
        target_index (dict): Index of the target sequence.

    Returns:
        dict: Dictionary of k-mers and their respective positions in the query and target sequences.
    """
    positions = {}
    for kmer, query_pos in kmers:
        if kmer in target_index:
            for target_pos in target_index[kmer]:
                if kmer not in positions:
                    positions[kmer] = []
                positions[kmer].append((target_pos, query_pos))  # position in target and in query
    return positions

# Function to find double-hits on the same diagonal and without strict overlap
def find_double_hits(kmer_positions, max_distance):
    """
    Identify double-hits between the query and target sequences.

    Args:
        kmer_positions (dict): The positions of k-mers in the target sequence.
        max_distance (int): Maximum allowable distance between two k-mers for a double-hit.
    
    Returns:
        list: A list of double-hit positions.
    """
    double_hits = []
    kmer_list = list(kmer_positions.items())
    
    for i in range(len(kmer_list)):
        kmer1, positions1 = kmer_list[i]
        for pos1, query_pos1 in positions1:
            for j in range(i + 1, len(kmer_list)):
                kmer2, positions2 = kmer_list[j]
                for pos2, query_pos2 in positions2:
                    # Check diagonal
                    if (pos1 - query_pos1) == (pos2 - query_pos2):
                        # Check for strict non-overlap and distance
                        if abs(pos2 - pos1) <= max_distance:
                            # Avoid hits that are exactly in the same location
                            if pos1 != pos2 or query_pos1 != query_pos2:
                                # Double hit detected
                                double_hits.append((kmer1, pos1, query_pos1, kmer2, pos2, query_pos2))
    return double_hits

# Function to evaluate the initial score of a double-hit
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
    score = 0
    for a, b in zip(kmer1, kmer2):
        score += get_blosum62_score(a, b, blosum_matrix)
    return score

# Extend a double-hit with gapped alignment
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2, X_g=1000000):
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
        X_g (int): Maximum score drop to stop extension.
    
    Returns:
        tuple: Alignment score and the final alignment.
    """
    m, n = len(seq_query), len(seq_target)
    
    # Initialize current and previous score rows
    previous_row = np.zeros(n + 1)
    current_row = np.zeros(n + 1)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)  # Still need full traceback matrix for recovery
    
    # Initialize best score tracking
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Ensure valid BLOSUM62 score calculation
            try:
                match = previous_row[j - 1] + get_blosum62_score(seq_query[i - 1], seq_target[j - 1], blosum_matrix)
            except KeyError:
                match = -4  # Default mismatch penalty for unknown pairs

            # Handle gap penalties (gap in query or gap in target)
            gap_query = current_row[j - 1] + (gap_extension_penalty if traceback_matrix[i][j - 1] == 2 else gap_open_penalty)
            gap_target = previous_row[j] + (gap_extension_penalty if traceback_matrix[i - 1][j] == 1 else gap_open_penalty)

            # Calculate best score for the current cell
            current_score = max(0, match, gap_query, gap_target)

            current_row[j] = current_score

            # Update traceback matrix
            if current_score == match:
                traceback_matrix[i][j] = 0  # Diagonal (match/mismatch)
            elif current_score == gap_query:
                traceback_matrix[i][j] = 2  # Left (gap in query)
            elif current_score == gap_target:
                traceback_matrix[i][j] = 1  # Up (gap in target)

            # Track max score
            if current_score > max_score:
                max_score = current_score
                max_pos = (i, j)

            # Check if current score has dropped more than X_g below max_score
            if max_score - current_score > X_g:
                break  # Stop the extension if score drops too much

        # Swap rows: move current_row to previous_row, and reset current_row
        previous_row, current_row = current_row, np.zeros(n + 1)

    # Traceback to recover the alignment (this part doesn't use the score matrix directly)
    aligned_query, aligned_target = [], []
    i, j = max_pos

    while i > 0 and j > 0 and traceback_matrix[i][j] != -1:  # Ensure valid bounds and traceback
        if traceback_matrix[i][j] == 0:
            aligned_query.append(seq_query[i - 1])
            aligned_target.append(seq_target[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:
            aligned_query.append(seq_query[i - 1])
            aligned_target.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 2:
            aligned_query.append('-')
            aligned_target.append(seq_target[j - 1])
            j -= 1

    aligned_query = aligned_query[::-1]
    aligned_target = aligned_target[::-1]
    final_alignment = list(zip(aligned_query, aligned_target))

    return max_score, final_alignment

# Filter alignments by score with distance selection
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
    
    # Si aucun double-hit n'est détecté, on retourne une liste vide
    if not double_hits:
        return alignments

    # On trouve le double-hit avec la plus grande distance
    best_hit = None
    max_distance = 0

    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        distance = abs(pos2 - pos1)
        
        if distance > max_distance:
            max_distance = distance
            best_hit = hit

    # Si on a trouvé un meilleur double-hit
    if best_hit:
        # On étend l'alignement à partir de ce double-hit
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = best_hit
        score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)
        bit_score = calculate_bit_score(score)
        
        if bit_score >= bit_score_threshold:
            alignments.append((score, alignment))

    return alignments

def filter_duplicate_alignments(alignments):
    """
    Remove duplicate alignments based on positions.

    Args:
        alignments (list): List of alignments.
    
    Returns:
        list: Unique alignments with duplicates removed.
    """
    unique_alignments = []
    seen_positions = set()

    for score, alignment in alignments:
        # Create a unique key based on the alignment positions
        alignment_key = tuple(alignment)
        
        if alignment_key not in seen_positions:
            unique_alignments.append((score, alignment))
            seen_positions.add(alignment_key)

    return unique_alignments

# Fonction pour calculer le bit score
def calculate_bit_score(raw_score, lambda_param=0.318, K=0.134):
    """
    Calculate the bit score from a raw score.

    Args:
        raw_score (int): Raw alignment score.
        lambda_param (float): Lambda parameter for scoring.
        K (float): K parameter for scoring.
    
    Returns:
        float: The calculated bit score.
    """
    return (lambda_param * raw_score - math.log(K)) / math.log(2)

# Calculate E-value from score
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
    e_values = []
    m = len(seq_query)
    n = len_database

    for raw_score, alignment in alignments:
        bit_score = calculate_bit_score(raw_score)  # Calcul du bit score
        e_value = calculate_e_value_from_bitscore(bit_score, m, n)  # Calcul de l'E-value à partir du bit score
        e_values.append((e_value, raw_score, alignment))
    
    return e_values

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
    filtered_alignments = [align for align in e_values if align[0] <= e_value_threshold]
    return filtered_alignments

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
    query_aligned = []
    target_aligned = []
    match_line = []

    for q_res, t_res in alignment:
        query_aligned.append(q_res if q_res != '-' else '-')
        target_aligned.append(t_res if t_res != '-' else '-')

        if q_res == t_res:
            match_line.append('|')  # perfect match
        elif get_blosum62_score(q_res, t_res, MatrixInfo.blosum62) > 0:
            match_line.append(':')  # reasonable match
        else:
            match_line.append('.')  # mismatch or weak match

    query_str = "".join(query_aligned)
    match_str = "".join(match_line)
    target_str = "".join(target_aligned)

    return query_str, match_str, target_str

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
    with open(output_file, "r") as f:
        lines = f.readlines()

    # Placeholder for storing the alignments
    alignments = []
    current_alignment = []

    # Parse the file and collect alignments
    for line in lines:
        if line.startswith("Alignment"):  # New alignment block
            if current_alignment:  # If there's an existing alignment, save it
                alignments.append(current_alignment)
            current_alignment = [line]  # Start a new alignment
        else:
            current_alignment.append(line)  # Append lines to the current alignment
    
    if current_alignment:  # Don't forget the last alignment
        alignments.append(current_alignment)

    # Sort alignments by the E-value, ensuring that we only sort those that have an "Expect =" line
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
                    evalue_str = line.split('=')[-1].strip()
                    return float(evalue_str.replace("e", "E"))
                except ValueError:
                    return float('inf')  # If we can't parse, assume a very high E-value
        return float('inf')  # If no E-value found, assume a very high E-value

    # Sort alignments by the extracted E-value (from smallest to largest E-value)
    alignments.sort(key=extract_evalue)

    # Rewrite the sorted alignments back to the file
    with open(output_file, "w") as f:
        for alignment in alignments:
            f.writelines(alignment)

# Function to load FASTA database
def load_fasta_database(fasta_file):
    """
    Load sequences from a FASTA file to use as a database.

    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        list: List of sequences from the file.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def align_sequence(seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database, output_file):
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
    # Extract k-mers
    kmers = extract_kmers(seq_query, k)
    target_index = index_target_sequence(seq_target, k)
    kmer_positions = find_kmer_positions(kmers, target_index)
    double_hits = find_double_hits(kmer_positions, max_distance)
    
    # Filter and align
    alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, bit_score_threshold)
    alignments = filter_duplicate_alignments(alignments)
    e_values = calculate_e_values(alignments, seq_query, len_database)
    
    # Filter based on E-value
    significant_alignments = filter_by_e_value(e_values, e_value_threshold)
    
    if significant_alignments:
        return alignments, e_values
    return None

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
    len_database = sum(len(seq) for seq in database_sequences)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        # Submit tasks for parallel processing
        future_to_sequence = {
            executor.submit(align_sequence, seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database, e_value_threshold): seq_target
            for seq_target in database_sequences
        }

        for future in concurrent.futures.as_completed(future_to_sequence):
            seq_target = future_to_sequence[future]
            try:
                result = future.result()
                if result:
                    alignments, e_values = result
                    display_blast_like_results(alignments, e_values, seq_query, seq_target, output_file)
            except Exception as exc:
                print(f"Generated an exception for {seq_target}: {exc}")

def load_query_sequence(fasta_file):
    """
    Load the query sequence from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        str: The query sequence.
    """
    # Lire le fichier FASTA et renvoyer la première séquence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq)

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
    # Create BLASTX command
    blastx_cline = NcbiblastxCommandline(query=query_fasta, db=protein_db, evalue=e_value_threshold, outfmt=outfmt, out=output_file)
    print(f"Running BLASTX command: {blastx_cline}")
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
    # Create BLASTN command
    blastn_cline = NcbiblastnCommandline(query=query_fasta, db=nucleotide_db, evalue=e_value_threshold, outfmt=outfmt, out=output_file)
    print(f"Running BLASTN command: {blastn_cline}")
    stdout, stderr = blastn_cline()

def parse_blast_results(blast_xml_file):
    """
    Parse BLAST XML results and print alignments.

    Args:
        blast_xml_file (str): Path to BLAST XML file.
    """
    with open(blast_xml_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    print(f"****Alignment****")
                    print(f"sequence: {alignment.title}")
                    print(f"length: {alignment.length}")
                    print(f"e-value: {hsp.expect}")
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...\n")


# Example usage with FASTA database
if __name__ == "__main__":
    # Récupérer les arguments de ligne de commande
    (method, query, database, output, e_value_threshold, outfmt, k, max_distance, bit_score_threshold, num_threads) = parse_arguments()

    # Chronométrage de l'exécution
    start = time.time()

    if method == "blastp":
        # Charger la séquence de requête et de la base de données
        seq_query = load_query_sequence(query)
        database_sequences = load_fasta_database(database)
        
        # Charger la matrice BLOSUM62
        blosum62 = MatrixInfo.blosum62

        # Exécution de l'algorithme personnalisé BLASTP en parallèle
        run_parallel_alignments(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold, num_threads, e_value_threshold, output)
        
        # Vérifier si le fichier de sortie a été généré
        if os.path.exists(output):
            sort_output_file(output)
        else:
            print(f"No alignments found. Output file '{output}' not created.")

    elif method == "blastx":
        # Exécution de BLASTX avec Biopython
        run_blastx(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)

    elif method == "blastn":
        # Exécution de BLASTN avec Biopython
        run_blastn(query, database, output, e_value_threshold=e_value_threshold, outfmt=outfmt)

    else:
        print(f"Unknown method: {method}. Please choose from 'blastp', 'blastx', or 'blastn'.")

    # Fin du chronomètre
    end = time.time()
    print(f"Execution time for {method} = {end - start} seconds")
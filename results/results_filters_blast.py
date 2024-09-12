import math
import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import concurrent.futures
import time

# Function to obtain the BLOSUM62 score for a pair of amino acids
def get_blosum62_score(a, b, blosum_matrix):
    """
    Returns the BLOSUM62 score for a pair of amino acids.

    Args:
        a (str): First amino acid.
        b (str): Second amino acid.
        blosum_matrix (dict): BLOSUM62 matrix loaded via Biopython.

    Returns:
        int: The score of the amino acid pair or -4 if the pair does not exist.
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
    Extracts k-mers from a given sequence.

    Args:
        sequence (str): The sequence to process.
        k (int): The length of the k-mers.

    Returns:
        list: List of tuples containing the k-mer and its position.
    """
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append((kmer, i))  # Preserve the original index
    return kmers

# Indexing the target sequence for fast k-mer search
def index_target_sequence(sequence, k):
    """
    Indexes a target sequence by k-mers for fast searching.

    Args:
        sequence (str): Target sequence to index.
        k (int): Length of the k-mers.

    Returns:
        dict: Dictionary where keys are k-mers and values are positions in the sequence.
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
    Identifies double-hits in the target sequence.

    Double-hits are pairs of k-mers on the same diagonal, that do not overlap,
    and whose distance is less than a defined threshold.

    Args:
        kmer_positions (dict): Positions of k-mers in the target sequence.
        max_distance (int): Maximum allowed distance between two hits to be considered a double-hit.

    Returns:
        list: List of tuples representing the detected double-hits.
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
    Calculates the initial score of a double-hit based on the BLOSUM62 matrix.

    Args:
        kmer1 (str): First k-mer.
        kmer2 (str): Second k-mer.
        blosum_matrix (dict): BLOSUM62 matrix.

    Returns:
        int: The score of the double-hit.
    """
    score = 0
    for a, b in zip(kmer1, kmer2):
        score += get_blosum62_score(a, b, blosum_matrix)
    return score

# Extend a double-hit with gapped alignment
import time

def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2, X_g=10000):
    """
    Extends a double-hit into an alignment, using Smith-Waterman for local alignment with gap handling.
    The extension is stopped if the score drops more than X_g below the best score.

    Args:
        seq_query (str): Query sequence.
        seq_target (str): Target sequence.
        query_pos (int): Hit position in the query sequence.
        target_pos (int): Hit position in the target sequence.
        blosum_matrix (dict): BLOSUM62 matrix.
        gap_open_penalty (int): Penalty for gap opening.
        gap_extension_penalty (int): Penalty for gap extension.
        X_g (int): Threshold for score drop before stopping the extension.

    Returns:
        tuple: The alignment score and the alignment itself as a list of tuples.
    """
    
    m, n = len(seq_query), len(seq_target)
    
    # Start total timing
    #total_start_time = time.time()
    
    # Initialize current and previous score rows
    #init_start_time = time.time()
    previous_row = np.zeros(n + 1)
    current_row = np.zeros(n + 1)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)  # Still need full traceback matrix for recovery
    #init_end_time = time.time()
    
    # Initialize best score tracking
    max_score = 0
    max_pos = (0, 0)

    # Fill the score matrix using Smith-Waterman logic (only using two rows at a time)
    total_match_time = 0
    total_gap_query_time = 0
    total_gap_target_time = 0
    total_update_time = 0
    total_max_score_time = 0

    #fill_start_time = time.time()
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Ensure valid BLOSUM62 score calculation
            #match_start_time = time.time()
            try:
                match = previous_row[j - 1] + get_blosum62_score(seq_query[i - 1], seq_target[j - 1], blosum_matrix)
            except KeyError:
                match = -4  # Default mismatch penalty for unknown pairs
            #match_end_time = time.time()
            #total_match_time += match_end_time - match_start_time

            # Handle gap penalties (gap in query or gap in target)
            #gap_query_start_time = time.time()
            gap_query = current_row[j - 1] + (gap_extension_penalty if traceback_matrix[i][j - 1] == 2 else gap_open_penalty)
            #gap_query_end_time = time.time()
            #total_gap_query_time += gap_query_end_time - gap_query_start_time

            #gap_target_start_time = time.time()
            gap_target = previous_row[j] + (gap_extension_penalty if traceback_matrix[i - 1][j] == 1 else gap_open_penalty)
            #gap_target_end_time = time.time()
            #total_gap_target_time += gap_target_end_time - gap_target_start_time

            # Calculate best score for the current cell
            #update_start_time = time.time()
            current_score = max(0, match, gap_query, gap_target)
            current_row[j] = current_score
            #update_end_time = time.time()
            #total_update_time += update_end_time - update_start_time

            # Update traceback matrix
            if current_score == match:
                traceback_matrix[i][j] = 0  # Diagonal (match/mismatch)
            elif current_score == gap_query:
                traceback_matrix[i][j] = 2  # Left (gap in query)
            elif current_score == gap_target:
                traceback_matrix[i][j] = 1  # Up (gap in target)

            # Track max score
            #max_score_start_time = time.time()
            if current_score > max_score:
                max_score = current_score
                max_pos = (i, j)
            #max_score_end_time = time.time()
            #total_max_score_time += max_score_end_time - max_score_start_time

            # Check if current score has dropped more than X_g below max_score
            if max_score - current_score > X_g:
                break  # Stop the extension if score drops too much

        # Swap rows: move current_row to previous_row, and reset current_row
        previous_row, current_row = current_row, np.zeros(n + 1)
    #fill_end_time = time.time()

    # Traceback to recover the alignment (this part doesn't use the score matrix directly)
    #traceback_start_time = time.time()
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
    #traceback_end_time = time.time()

    # Total time
    #total_end_time = time.time()

    # Print time taken for different stages
    #print(f"Time for initializing matrices: {init_end_time - init_start_time:.4f} sec")
    #print(f"Time for filling score matrix: {fill_end_time - fill_start_time:.4f} sec")
    #print(f"  Time for match calculations: {total_match_time:.4f} sec")
    #print(f"  Time for gap_query calculations: {total_gap_query_time:.4f} sec")
    #print(f"  Time for gap_target calculations: {total_gap_target_time:.4f} sec")
    #print(f"  Time for score updates: {total_update_time:.4f} sec")
    #print(f"  Time for max score updates: {total_max_score_time:.4f} sec")
    #print(f"Time for traceback: {traceback_end_time - traceback_start_time:.4f} sec")
    #print(f"Total time for alignment: {total_end_time - total_start_time:.4f} sec")

    return max_score, final_alignment

# Filter alignments by score with distance selection
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, bit_score_threshold):
    """
    Filtre et étend les alignements basés sur les double-hits détectés,
    en sélectionnant uniquement le double-hit avec la plus grande distance entre k-mers.

    Args:
        double_hits (list): Liste des double-hits détectés.
        seq_query (str): Séquence query.
        seq_target (str): Séquence target.
        blosum_matrix (dict): Matrice BLOSUM62.
        bit_score_threshold (float): Seuil pour filtrer les alignements basés sur le bit score.

    Returns:
        list: Liste des alignements filtrés (un seul alignement par séquence).
    """
    alignments = []
    rejected_by_bit_score = 0  # Compteur des alignements rejetés par bit score
    total_double_hits = len(double_hits)

    if not double_hits:
        return alignments

    best_hit = None
    max_distance = 0

    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        distance = abs(pos2 - pos1)
        
        if distance > max_distance:
            max_distance = distance
            best_hit = hit

    if best_hit:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = best_hit
        score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)
        
        bit_score = calculate_bit_score(score)

        # Afficher le score brut et le bit score pour analyse
        #print(f"Double-hit entre {kmer1} et {kmer2} | Score brut: {score}, Bit score: {bit_score:.2f}")

        if bit_score >= bit_score_threshold:
            alignments.append((score, alignment))
        else:
            # Si l'alignement est rejeté à cause du bit score
            rejected_by_bit_score += 1

    # Afficher les détails après traitement de tous les double-hits
    #print(f"Double-hits totaux pour cette séquence: {total_double_hits}")
    #print(f"Alignements rejetés par bit score: {rejected_by_bit_score}")

    return alignments

def filter_duplicate_alignments(alignments):
    """
    Filters out duplicate alignments, keeping only unique ones.

    Args:
        alignments (list): List of alignments.

    Returns:
        list: List of unique alignments.
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
    Calcule le bit score à partir du score brut.
    
    Args:
        raw_score (int): Score brut (raw score).
        lambda_param (float): Paramètre lambda pour la distribution des scores.
        K (float): Paramètre K pour la distribution des scores.

    Returns:
        float: Le bit score.
    """
    return (lambda_param * raw_score - math.log(K)) / math.log(2)

# Calculate E-value from score
def calculate_e_value_from_bitscore(bit_score, m, n):
    """
    Calcule l'E-value à partir du bit score.

    Args:
        bit_score (float): Le bit score de l'alignement.
        m (int): Longueur de la séquence query.
        n (int): Longueur totale de la base de données.

    Returns:
        float: E-value calculée.
    """
    return m * n * math.pow(2, -bit_score)

# Calculate E-values for all alignments
def calculate_e_values(alignments, seq_query, len_database):
    """
    Calcule les E-values pour une liste d'alignements en utilisant les bit scores.

    Args:
        alignments (list): Liste des alignements avec leur score brut.
        seq_query (str): Séquence de la query.
        len_database (int): Longueur totale de la base de données.

    Returns:
        list: Liste des E-values associées à chaque alignement.
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
def filter_by_e_value(e_values, threshold=0.01):
    """
    Filters alignments based on an E-value threshold.

    Args:
        e_values (list): List of alignments with their E-values.
        threshold (float): E-value threshold for filtering alignments.

    Returns:
        list: List of significant alignments with E-value below the threshold.
    """
    filtered_alignments = [align for align in e_values if align[0] <= threshold]
    return filtered_alignments

# Function to format the alignment for BLAST-like display
def format_alignment(seq_query, seq_target, alignment):
    """
    Formats an alignment for BLAST-style display.

    Args:
        seq_query (str): Query sequence.
        seq_target (str): Target sequence.
        alignment (list): Alignment to format.

    Returns:
        tuple: Formatted strings for the query, match line, and target.
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
def display_blast_like_results(alignments, e_values, seq_query, seq_target):
    """
    Displays the alignment results in a BLAST-like format, including score, E-value, identities, positives, and gaps.

    Args:
        alignments (list): List of alignments with their score.
        e_values (list): List of E-values associated with the alignments.
        seq_query (str): The query sequence.
        seq_target (str): The target sequence.
    """
    # Parameters for formatting
    line_length = 60  # Number of characters per line in the alignment display

    print(f"{len(e_values)} significant alignments found:\n")
    
    for i, (e_value, raw_score, alignment) in enumerate(e_values, 1):
        bit_score = calculate_bit_score(raw_score)
        query_str, match_str, target_str = format_alignment(seq_query, seq_target, alignment)

        # Calculate identities, positives, and gaps
        identities = sum(1 for q, t in alignment if q == t)
        positives = sum(1 for q, t in alignment if get_blosum62_score(q, t, blosum62) > 0)
        gaps = sum(1 for q, t in alignment if q == '-' or t == '-')

        # Calculate the percentages
        total_length = len(alignment)
        identity_percentage = (identities / total_length) * 100
        positive_percentage = (positives / total_length) * 100
        gap_percentage = (gaps / total_length) * 100

        # Display summary information
        print(f"Alignment {i}:")
        print(f" Score = {raw_score}, bits=({bit_score:.2f}), Expect = {e_value:.2e}")
        print(f" Identities = {identities}/{total_length} ({identity_percentage:.0f}%), "
              f"Positives = {positives}/{total_length} ({positive_percentage:.0f}%), "
              f"Gaps = {gaps}/{total_length} ({gap_percentage:.0f}%)")
        print("\n")

        # Display alignment block by block
        for block_start in range(0, total_length, line_length):
            block_end = block_start + line_length
            query_block = query_str[block_start:block_end]
            match_block = match_str[block_start:block_end]
            target_block = target_str[block_start:block_end]

            query_label = f"Query {block_start + 1}".ljust(9)
            target_label = f"Sbjct {block_start + 1}".ljust(9)
            
            print(f"{query_label} {query_block}")
            print(f"{''.ljust(9)} {match_block}")
            print(f"{target_label} {target_block}")
            print()

        print("\n" + "-"*50 + "\n")

# Function to load FASTA database
def load_fasta_database(fasta_file):
    """
    Loads a FASTA database from a file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        list: List of sequences in the database.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def align_sequence(seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database):
    """
    Aligns a single target sequence to the query sequence by:
    - Extracting k-mers
    - Finding double-hits
    - Extending and filtering alignments

    Args:
        seq_target (str): Target sequence.
        seq_query (str): Query sequence.
        k (int): Length of k-mers.
        max_distance (int): Maximum allowed distance between two k-mers for a double-hit.
        blosum62 (dict): BLOSUM62 scoring matrix.
        bit_score_threshold (float): Threshold for filtering alignments based on the bit score.
        len_database (int): Total length of the database (used to compute E-value).

    Returns:
        tuple: Alignments and their corresponding E-values, or None if no significant alignment is found.
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
    significant_alignments = filter_by_e_value(e_values, threshold=0.00001)
    
    if significant_alignments:
        return alignments, e_values
    return None

# Function to run alignment in parallel using multiple threads
def run_parallel_alignments(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold, num_threads):
    """
    Executes the alignment process in parallel for a set of target sequences using multiple threads.

    Args:
        seq_query (str): Query sequence.
        database_sequences (list): List of target sequences from the database.
        k (int): Length of k-mers.
        max_distance (int): Maximum allowed distance between two k-mers for a double-hit.
        blosum62 (dict): BLOSUM62 scoring matrix.
        bit_score_threshold (float): Threshold for filtering alignments based on the bit score.
        num_threads (int): Number of threads for parallel execution.
    """
    global initial_sequence_count, kmer_filtered_count, double_hits_detected_count, filtered_alignments_count, final_significant_alignments_count
    
    # Initialize counters
    initial_sequence_count = 0
    kmer_filtered_count = 0
    double_hits_detected_count = 0
    filtered_alignments_count = 0
    final_significant_alignments_count = 0

    # Count the number of sequences in the database
    initial_sequence_count = len(database_sequences)
    print(f"Total number of initial sequences: {initial_sequence_count}")

    # Calculate the total length of the database
    len_database = sum(len(seq) for seq in database_sequences)
    
    # Use a ProcessPoolExecutor to run the alignments in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        # Submit the alignment tasks to be run in parallel
        future_to_sequence = {
            executor.submit(align_sequence, seq_target, seq_query, k, max_distance, blosum62, bit_score_threshold, len_database): seq_target
            for seq_target in database_sequences
        }

        # Retrieve and process the results as they complete
        for future in concurrent.futures.as_completed(future_to_sequence):
            seq_target = future_to_sequence[future]
            try:
                result = future.result()
                if result:
                    alignments, e_values = result
                    display_blast_like_results(alignments, e_values, seq_query, seq_target)
            except Exception as exc:
                print(f"Generated an exception for {seq_target}: {exc}")

    # Display the filtering results summary
    print(f"Number of sequences after k-mer extraction: {kmer_filtered_count}")
    print(f"Number of sequences with double-hits detected: {double_hits_detected_count}")
    print(f"Number of sequences with filtered alignments: {filtered_alignments_count}")
    print(f"Number of sequences with significant alignments: {final_significant_alignments_count}")


# --- Global counters ---
initial_sequence_count = 0
kmer_filtered_count = 0
double_hits_detected_count = 0
filtered_alignments_count = 0
final_significant_alignments_count = 0

# Main function that aligns sequences and displays the filtering results
def align_and_count(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold):
    """
    Main function that aligns a query sequence to each target sequence in the database.
    It also tracks and displays the results of k-mer filtering, double-hit detection, 
    alignment filtering, and significant alignments.

    Args:
        seq_query (str): Query sequence.
        database_sequences (list): List of target sequences from the database.
        k (int): Length of k-mers.
        max_distance (int): Maximum allowed distance between two k-mers for a double-hit.
        blosum62 (dict): BLOSUM62 scoring matrix.
        bit_score_threshold (float): Threshold for filtering alignments based on the bit score.
    """
    global initial_sequence_count, kmer_filtered_count, double_hits_detected_count, filtered_alignments_count, final_significant_alignments_count

    # Count initial sequences
    initial_sequence_count = len(database_sequences)
    print(f"Total number of initial sequences: {initial_sequence_count}")

    # Total database length
    len_database = sum(len(seq) for seq in database_sequences)

    total_rejected_by_bit_score = 0  # Count of alignments rejected due to low bit score

    # Iterate over each target sequence
    for seq_target in database_sequences:
        # Extract k-mers
        kmers = extract_kmers(seq_query, k)
        if kmers:
            kmer_filtered_count += 1

        # Index target sequence and find k-mer positions
        target_index = index_target_sequence(seq_target, k)
        kmer_positions = find_kmer_positions(kmers, target_index)
        double_hits = find_double_hits(kmer_positions, max_distance)
        if double_hits:
            double_hits_detected_count += 1

        # Filter and align
        alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, bit_score_threshold)
        alignments = filter_duplicate_alignments(alignments)
        
        # Count alignments rejected by bit score
        rejected_by_bit_score = len(double_hits) - len(alignments)
        total_rejected_by_bit_score += rejected_by_bit_score

        if alignments:
            filtered_alignments_count += 1

        # Calculate E-values
        e_values = calculate_e_values(alignments, seq_query, len_database)
        significant_alignments = filter_by_e_value(e_values, threshold=0.00001)
        if significant_alignments:
            final_significant_alignments_count += 1

    # Display summary of filtering
    print(f"Number of sequences after k-mer extraction: {kmer_filtered_count}")
    print(f"Number of sequences with double-hits detected: {double_hits_detected_count}")
    print(f"Number of sequences with filtered alignments: {filtered_alignments_count}")
    print(f"Number of sequences with significant alignments: {final_significant_alignments_count}")
    print(f"Total number of alignments rejected due to low bit score: {total_rejected_by_bit_score}")

# Example of loading FASTA sequences and running alignment with counts
if __name__ == "__main__":
    # Sample query sequence (example)
    seq_query = "ADEPILVA"
    blosum62 = MatrixInfo.blosum62
    k = 3
    max_distance = 40
    bit_score_threshold = 0

    # Example database sequences
    database_sequences = ["ADEPILVA"]

    # Run the alignment and display filtering results
    align_and_count(seq_query, database_sequences, k, max_distance, blosum62, bit_score_threshold)
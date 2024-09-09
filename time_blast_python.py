import math
import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
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
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2, X_g=1000):
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
    
    start_extend_time = time.time()

    m, n = len(seq_query), len(seq_target)
    
    # Initialize score and traceback matrices
    score_matrix = np.zeros((m + 1, n + 1))
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)
    
    # Initialize best score tracking
    max_score = 0
    best_score = 0
    max_pos = (0, 0)

    # Fill the score matrix using Smith-Waterman logic
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + get_blosum62_score(seq_query[i - 1], seq_target[j - 1], blosum_matrix)
            gap_query = score_matrix[i][j - 1] + (gap_extension_penalty if traceback_matrix[i][j - 1] == 2 else gap_open_penalty)
            gap_target = score_matrix[i - 1][j] + (gap_extension_penalty if traceback_matrix[i - 1][j] == 1 else gap_open_penalty)

            # Calculate best score for the current cell
            current_score = max(0, match, gap_query, gap_target)
            score_matrix[i][j] = current_score

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
            #print(current_score, max_score)
            if max_score - current_score > X_g:
                break  # Stop the extension if score drops too much

    # Traceback to recover the alignment
    aligned_query, aligned_target = [], []
    i, j = max_pos

    while score_matrix[i][j] > 0:
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

    # Since traceback traverses in reverse, we need to reverse the results
    aligned_query = aligned_query[::-1]
    aligned_target = aligned_target[::-1]

    # Combine the results into a final alignment
    final_alignment = list(zip(aligned_query, aligned_target))

    end_extend_time = time.time()
    if (end_extend_time - start_extend_time) >= 1:
        print(f"Time extend_alignment = {end_extend_time - start_extend_time}")

    return max_score, final_alignment


# Filter alignments by score
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, bit_score_threshold):
    """
    Filtre et étend les alignements basés sur les double-hits détectés, avec un seuil basé sur le bit score.

    Args:
        double_hits (list): Liste des double-hits détectés.
        seq_query (str): Séquence query.
        seq_target (str): Séquence target.
        blosum_matrix (dict): Matrice BLOSUM62.
        bit_score_threshold (float): Seuil pour filtrer les alignements basés sur le bit score.

    Returns:
        list: Liste des alignements filtrés.
    """
    alignments = []
    score_value = []
    #print("==================Entering filter_alignments=================")
    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        initial_score = evaluate_double_hit(kmer1, kmer2, blosum_matrix)
        if initial_score not in score_value:
            score_value.append(initial_score)
            if initial_score <= 0:  # Conserver uniquement les alignements avec un score initial negatif
                #print(f"initial_score = {initial_score}")
                score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)
                bit_score = calculate_bit_score(score)
                if bit_score >= bit_score_threshold:
                    alignments.append((score, alignment))
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

# Example usage with FASTA database
if __name__ == "__main__":
    seq_query = "MGRLDGKVIILTAAAQGIGQAAALAFAREGAKVIATDINESKLQELEKYPGIQTRVLDVTKKKQIDQFANEVERLDVLFNVAGFVHHGTVLDCEEKDWDFSMNLNVRSMYLMIKAFLPKMLAQKSGNIINMSSVASSVKGVVNRCVYSTTKAAVIGLTKSVAADFIQQGIRCNCVCPGTVDTPSLQERIQARGNPEEARNDFLKRQKTGRFATAEEIAMLCVYLASDESAYVTGNPVIIDGGWSL"
    
    blosum62 = MatrixInfo.blosum62

    k = 3
    max_distance = 40
    bit_score_threshold = 22  # Seuil basé sur le bit score
    
    fasta_file = "subset_2000_sequences.fasta"
    database_sequences = load_fasta_database(fasta_file)
    #database_sequences = ["MDKVCAVFGGSRGIGRAVAQLMARKGYRLAVIARNLEGAKAAAGDLGGDHLAFSCDVAKEHDVQNTFEELEKHLGRVNFLVNAAGINRDGLLVRTKTEDMVSQLHTNLLGSMLTCKAAMRTMIQQQGGSIVNVGSIVGLKGNSGQSVYSASKGGLVGFSRALAKEVARKKIRVNVVAPGFVHTDMTKDLKEEHLKKNIPLGRFGETIEVAHAVVFLLESPYITGHVLVVDGGLQLIL"]
    #database_sequences = ["MGRLDGKVIILTAAAQGIGQAAALAFAREGAKVIATDINESKLQELEKYPGIQTRVLDVTKKKQIDQFANEVERLDVLFNVAGFVHHGTVLDCEEKDWDFSMNLNVRSMYLMIKAFLPKMLAQKSGNIINMSSVASSVKGVVNRCVYSTTKAAVIGLTKSVAADFIQQGIRCNCVCPGTVDTPSLQERIQARGNPEEARNDFLKRQKTGRFATAEEIAMLCVYLASDESAYVTGNPVIIDGGWSL"]

    # Start measuring time for the whole process
    overall_start_time = time.time()

    # Initialize lists to store timing data for each step
    kmer_extraction_times = []
    target_index_times = []
    kmer_position_times = []
    double_hit_times = []
    filter_alignment_times = []
    total_times = []

    len_database = sum(len(seq) for seq in database_sequences)

    start_kmer_time = time.time()
    kmers = extract_kmers(seq_query, k)
    kmer_extraction_times.append(time.time() - start_kmer_time)

    for i, seq_target in enumerate(database_sequences):
            # Step 1: Check if the query is identical to the target
            if seq_query == seq_target:
                # Directly process identical sequences with perfect alignment
                perfect_alignment = list(zip(seq_query, seq_target))  # Perfect alignment with no gaps
                raw_score = sum(get_blosum62_score(q, t, blosum62) for q, t in perfect_alignment)  # Perfect score
                bit_score = calculate_bit_score(raw_score)  # Convert to bit score
                e_value = calculate_e_value_from_bitscore(bit_score, len(seq_query), len_database)  # Calculate E-value

                print("Perfect alignment detected (query == target).")
                display_blast_like_results([(raw_score, perfect_alignment)], [(e_value, raw_score, perfect_alignment)], seq_query, seq_target)
                continue  # Skip to the next target

            # Step 2: Time to index the target sequence
            start_target_index_time = time.time()
            target_index = index_target_sequence(seq_target, k)
            target_index_times.append(time.time() - start_target_index_time)

            # Step 3: Time to find k-mer positions
            start_kmer_pos_time = time.time()
            kmer_positions = find_kmer_positions(kmers, target_index)
            kmer_position_times.append(time.time() - start_kmer_pos_time)

            # Step 4: Time to find double hits
            start_double_hit_time = time.time()
            double_hits = find_double_hits(kmer_positions, max_distance)
            double_hit_times.append(time.time() - start_double_hit_time)

            # Step 5: Time to filter alignments
            start_filter_alignment_time = time.time()
            alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, bit_score_threshold)
            filter_alignment_times.append(time.time() - start_filter_alignment_time)

            alignments = filter_duplicate_alignments(alignments)
            e_values = calculate_e_values(alignments, seq_query, len_database)
            significant_alignments = filter_by_e_value(e_values, threshold=0.00001)

            if significant_alignments:
                display_blast_like_results(alignments, e_values, seq_query, seq_target)

            # Total time for this sequence
            total_times.append(time.time() - overall_start_time)

    # Total process time for all sequences
    overall_time = time.time() - overall_start_time

    # Generating stats
    print("\n--- Timing Statistics for 2000 Sequences ---")
    print(f"Total time for all sequences: {overall_time:.2f} seconds")
    print(f"Average time per sequence: {sum(total_times)/len(total_times):.2f} seconds")
    print(f"K-mer extraction: avg={sum(kmer_extraction_times)/len(kmer_extraction_times):.2f}s, min={min(kmer_extraction_times):.2f}s, max={max(kmer_extraction_times):.2f}s, total={sum(kmer_extraction_times):.3f}s")
    print(f"Target indexing: avg={sum(target_index_times)/len(target_index_times):.2f}s, min={min(target_index_times):.2f}s, max={max(target_index_times):.2f}s, total={sum(target_index_times):.3f}s")
    print(f"K-mer position finding: avg={sum(kmer_position_times)/len(kmer_position_times):.2f}s, min={min(kmer_position_times):.2f}s, max={max(kmer_position_times):.2f}s, total={sum(kmer_position_times):.3f}s")
    print(f"Double hit finding: avg={sum(double_hit_times)/len(double_hit_times):.2f}s, min={min(double_hit_times):.2f}s, max={max(double_hit_times):.2f}s, total={sum(double_hit_times):.3f}s")
    print(f"Alignment filtering: avg={sum(filter_alignment_times)/len(filter_alignment_times):.2f}s, min={min(filter_alignment_times):.2f}s, max={max(filter_alignment_times):.2f}s, total={sum(filter_alignment_times):.3f}s")
    print(f"Total sequence processing: avg={sum(total_times)/len(total_times):.2f}s, min={min(total_times):.2f}s, max={max(total_times):.2f}s")

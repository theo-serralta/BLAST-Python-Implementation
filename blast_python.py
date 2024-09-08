import math
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

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
                     gap_open_penalty=-11, gap_extension_penalty=-2):
    """
    Extends a double-hit into an alignment with gap penalties.

    Args:
        seq_query (str): Query sequence.
        seq_target (str): Target sequence.
        query_pos (int): Hit position in the query sequence.
        target_pos (int): Hit position in the target sequence.
        blosum_matrix (dict): BLOSUM62 matrix.
        gap_open_penalty (int): Penalty for gap opening.
        gap_extension_penalty (int): Penalty for gap extension.

    Returns:
        tuple: The alignment score and the alignment itself as a list of tuples.
    """
    score = 0
    alignment = []

    # Extend to the left
    left_score, left_alignment = extend_direction(seq_query, seq_target, query_pos - 1, target_pos - 1, 
                                                  blosum_matrix, gap_open_penalty, gap_extension_penalty, direction='left')
    
    # Extend to the right
    right_score, right_alignment = extend_direction(seq_query, seq_target, query_pos, target_pos, 
                                                    blosum_matrix, gap_open_penalty, gap_extension_penalty, direction='right')
    
    # Combine scores and alignments
    score = left_score + right_score
    alignment = left_alignment[::-1] + right_alignment  # Reverse the left alignment and combine with right
    
    return score, alignment

# Extend an alignment in one direction (left or right)
def extend_direction(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2, direction='right', max_dropoff=10):
    """
    Extends an alignment in a given direction, considering gap penalties.

    Args:
        seq_query (str): Query sequence.
        seq_target (str): Target sequence.
        query_pos (int): Position in the query sequence.
        target_pos (int): Position in the target sequence.
        blosum_matrix (dict): BLOSUM62 matrix.
        gap_open_penalty (int): Penalty for gap opening.
        gap_extension_penalty (int): Penalty for gap extension.
        direction (str): Direction of the extension ('left' or 'right').
        max_dropoff (int): Threshold value to stop the extension.

    Returns:
        tuple: The score of the extension and the alignment as a list of tuples.
    """
    score = 0
    alignment = []
    current_score = 0
    i = query_pos
    j = target_pos
    gap_opened_in_query = False
    gap_opened_in_target = False

    while i >= 0 and j >= 0 and i < len(seq_query) and j < len(seq_target):
        match_score = get_blosum62_score(seq_query[i], seq_target[j], blosum_matrix)
        
        # If match score is below the gap opening penalty, consider a gap
        if match_score <= gap_open_penalty:
            if gap_opened_in_query or gap_opened_in_target:
                score += gap_extension_penalty  # Gap extension
            else:
                score += gap_open_penalty  # Gap opening
                gap_opened_in_query = True
                gap_opened_in_target = True
            alignment.append(('-', seq_target[j]) if seq_query[i] == '-' else (seq_query[i], '-'))
        else:
            alignment.append((seq_query[i], seq_target[j]))
            score += match_score
            gap_opened_in_query = False
            gap_opened_in_target = False
        
        current_score += match_score

        if current_score < -max_dropoff:
            break

        # Increment i and j after processing the positions
        if direction == 'left':
            i -= 1
            j -= 1
        else:  # direction 'right'
            i += 1
            j += 1

    return score, alignment

# Filter alignments by score
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, threshold=-30):
    """
    Filters and extends alignments based on the detected double-hits.

    Args:
        double_hits (list): List of detected double-hits.
        seq_query (str): Query sequence.
        seq_target (str): Target sequence.
        blosum_matrix (dict): BLOSUM62 matrix.
        threshold (int): Score threshold for including an alignment.

    Returns:
        list: List of filtered alignments.
    """
    alignments = []
    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        initial_score = evaluate_double_hit(kmer1, kmer2, blosum_matrix)
        if initial_score >= threshold:
            score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)
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

# Calculate E-value from score
def calculate_e_value(score, m, n, lambda_param=0.318, K=0.134):
    """
    Calculates the E-value based on the score and sequence lengths.

    Args:
        score (int): Alignment score.
        m (int): Length of the query sequence.
        n (int): Total length of the database.
        lambda_param (float): Lambda parameter for the score distribution.
        K (float): K parameter for the score distribution.

    Returns:
        float: The calculated E-value.
    """
    e_value = K * m * n * math.exp(-lambda_param * score)
    return e_value

# Calculate E-values for all alignments
def calculate_e_values(alignments, seq_query, len_database):
    """
    Calculates E-values for a list of alignments.

    Args:
        alignments (list): List of alignments with their score.
        seq_query (str): Query sequence.
        len_database (int): Total length of the database.

    Returns:
        list: List of E-values associated with each alignment.
    """
    e_values = []
    m = len(seq_query)
    n = len_database

    for score, alignment in alignments:
        e_value = calculate_e_value(score, m, n)
        e_values.append((e_value, score, alignment))
    
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
def display_results(alignments, e_values):
    """
    Displays the alignment results in a formatted way.

    Args:
        alignments (list): List of alignments.
        e_values (list): List of E-values associated with the alignments.
    """
    print(f"{len(e_values)} significant alignments found:\n")
    
    for i, (e_value, score, alignment) in enumerate(e_values, 1):
        query_str, match_str, target_str = format_alignment(seq_query, seq_target, alignment)

        print(f"Alignment {i}:")
        print(f"Score: {score}, E-value: {e_value:.2e}")
        print(f"Query:  {query_str}")
        print(f"        {match_str}")
        print(f"Target: {target_str}")
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
    score_threshold = 15
    
    fasta_file = "subset_2000_sequences.fasta"
    database_sequences = load_fasta_database(fasta_file)
    
    len_database = sum(len(seq) for seq in database_sequences)
    
    kmers = extract_kmers(seq_query, k)

    for i, seq_target in enumerate(database_sequences):
        target_index = index_target_sequence(seq_target, k)
        kmer_positions = find_kmer_positions(kmers, target_index)
        double_hits = find_double_hits(kmer_positions, max_distance)
        alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, threshold=score_threshold)
        alignments = filter_duplicate_alignments(alignments)
        e_values = calculate_e_values(alignments, seq_query, len_database)
        significant_alignments = filter_by_e_value(e_values, threshold=0.00001)

        if significant_alignments:
            display_results(alignments, significant_alignments)

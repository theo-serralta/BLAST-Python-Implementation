import math
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

# Fonction pour obtenir le score entre deux acides aminés avec BLOSUM62
def get_blosum62_score(a, b, blosum_matrix):
    if (a, b) in blosum_matrix:
        return blosum_matrix[(a, b)]
    elif (b, a) in blosum_matrix:  # matrice symétrique
        return blosum_matrix[(b, a)]
    else:
        return -4  # pénalité par défaut pour mismatch/gap

# Extraction des k-mers
def extract_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append((kmer, i))  # Conserver l'indice d'origine
    return kmers

# Indexation de la séquence cible pour une recherche rapide des k-mers
def index_target_sequence(sequence, k):
    kmer_index = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmer_index:
            kmer_index[kmer] = []
        kmer_index[kmer].append(i)  # Conserver l'indice d'origine
    return kmer_index

# Trouver les positions des k-mers dans la séquence cible
def find_kmer_positions(kmers, target_index):
    positions = {}
    for kmer, query_pos in kmers:
        if kmer in target_index:
            for target_pos in target_index[kmer]:
                if kmer not in positions:
                    positions[kmer] = []
                positions[kmer].append((target_pos, query_pos))  # position dans target et dans query
    return positions

# Fonction pour trouver les double-hits
def find_double_hits(kmer_positions, max_distance):
    double_hits = []
    kmer_list = list(kmer_positions.items())
    
    for i in range(len(kmer_list)):
        kmer1, positions1 = kmer_list[i]
        for pos1, query_pos1 in positions1:
            for j in range(i + 1, len(kmer_list)):
                kmer2, positions2 = kmer_list[j]
                for pos2, query_pos2 in positions2:
                    if abs(pos2 - pos1) <= max_distance:
                        # Double hit détecté uniquement si la distance est respectée
                        double_hits.append((kmer1, pos1, query_pos1, kmer2, pos2, query_pos2))
    return double_hits


# Fonction pour évaluer le score initial d'un double-hit
def evaluate_double_hit(kmer1, kmer2, blosum_matrix):
    score = 0
    for a, b in zip(kmer1, kmer2):
        score += get_blosum62_score(a, b, blosum_matrix)
    return score

# Extension d'un double-hit avec gapped alignment
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2):
    score = 0
    alignment = []

    # Étendre vers la gauche
    left_score, left_alignment = extend_direction(seq_query, seq_target, query_pos - 1, target_pos - 1, 
                                                  blosum_matrix, gap_open_penalty, gap_extension_penalty, direction='left')
    
    # Étendre vers la droite
    right_score, right_alignment = extend_direction(seq_query, seq_target, query_pos, target_pos, 
                                                    blosum_matrix, gap_open_penalty, gap_extension_penalty, direction='right')
    
    # Combiner les scores et les alignements
    score = left_score + right_score
    alignment = left_alignment[::-1] + right_alignment  # Reverse the left alignment and combine with right
    
    return score, alignment



# Fonction pour étendre dans une direction (gauche ou droite)
def extend_direction(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2, direction='right', max_dropoff=10):
    score = 0
    alignment = []
    current_score = 0
    i = query_pos
    j = target_pos
    gap_opened_in_query = False
    gap_opened_in_target = False

    while i >= 0 and j >= 0 and i < len(seq_query) and j < len(seq_target):
        match_score = get_blosum62_score(seq_query[i], seq_target[j], blosum_matrix)
        
        # Si le match score est inférieur à la pénalité d'ouverture, on considère un gap
        if match_score <= gap_open_penalty:
            if gap_opened_in_query or gap_opened_in_target:
                score += gap_extension_penalty  # Extension du gap
            else:
                score += gap_open_penalty  # Ouverture du gap
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

        # Incrémenter i et j après avoir traité les positions
        if direction == 'left':
            i -= 1
            j -= 1
        else:  # direction 'right'
            i += 1
            j += 1

    return score, alignment

# Filtrer les alignements par score
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, threshold=-30):
    alignments = []
    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        initial_score = evaluate_double_hit(kmer1, kmer2, blosum_matrix)
        if initial_score >= threshold:
            score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix, gap_open_penalty=-11, gap_extension_penalty=-2)

            alignments.append((score, alignment))
    return alignments

def filter_duplicate_alignments(alignments):
    unique_alignments = []
    seen_positions = set()

    for score, alignment in alignments:
        # Créer une clé unique basée sur les positions de l'alignement
        alignment_key = tuple(alignment)
        
        if alignment_key not in seen_positions:
            unique_alignments.append((score, alignment))
            seen_positions.add(alignment_key)

    return unique_alignments

# Calcul de l'E-value à partir du score
def calculate_e_value(score, m, n, lambda_param=0.318, K=0.134):
    """
    Calcule l'E-value basée sur un score normalisé et la longueur des séquences query et target
    """
    e_value = K * m * n * math.exp(-lambda_param * score)
    return e_value

# Calcul des E-values pour tous les alignements
def calculate_e_values(alignments, seq_query, len_database):
    e_values = []
    m = len(seq_query)
    n = len_database

    # Boucle pour calculer l'E-value pour chaque alignement
    for score, alignment in alignments:
        e_value = calculate_e_value(score, m, n)  # n reste constant ici
        e_values.append((e_value, score, alignment))
    
    return e_values

# Filtrer les alignements par E-value
def filter_by_e_value(e_values, threshold=0.01):
    filtered_alignments = [align for align in e_values if align[0] <= threshold]
    return filtered_alignments

# Fonction pour formater l'alignement pour un affichage BLAST-like
def format_alignment(seq_query, seq_target, alignment):
    query_aligned = []
    target_aligned = []
    match_line = []

    for q_res, t_res in alignment:
        query_aligned.append(q_res if q_res != '-' else '-')
        target_aligned.append(t_res if t_res != '-' else '-')

        if q_res == t_res:
            match_line.append('|')  # correspondance parfaite
        elif get_blosum62_score(q_res, t_res, MatrixInfo.blosum62) > 0:
            match_line.append(':')  # correspondance raisonnable
        else:
            match_line.append('.')  # mismatch ou faible correspondance

    query_str = "".join(query_aligned)
    match_str = "".join(match_line)
    target_str = "".join(target_aligned)

    return query_str, match_str, target_str

# Affichage des résultats formatés
def display_results(alignments, e_values):
    print(f"{len(e_values)} alignements significatifs trouvés :\n")
    
    for i, (e_value, score, alignment) in enumerate(e_values, 1):
        # Formater l'alignement
        query_str, match_str, target_str = format_alignment(seq_query, seq_target, alignment)

        print(f"Alignement {i}:")
        print(f"Score: {score}, E-value: {e_value:.2e}")
        print(f"Query:  {query_str}")
        print(f"        {match_str}")
        print(f"Target: {target_str}")
        print("\n" + "-"*50 + "\n")

# Fonction pour charger la base de données FASTA
def load_fasta_database(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# Exemple d'utilisation avec la base de données FASTA
if __name__ == "__main__":
    # Charger la séquence query (fournie par l'utilisateur)
    seq_query = "MGRLDGKVIILTAAAQGIGQAAALAFAREGAKVIATDINESKLQELEKYPGIQTRVLDVTKKKQIDQFANEVERLDVLFNVAGFVHHGTVLDCEEKDWDFSMNLNVRSMYLMIKAFLPKMLAQKSGNIINMSSVASSVKGVVNRCVYSTTKAAVIGLTKSVAADFIQQGIRCNCVCPGTVDTPSLQERIQARGNPEEARNDFLKRQKTGRFATAEEIAMLCVYLASDESAYVTGNPVIIDGGWSL"
    
    # Chargement de la matrice BLOSUM62 via Biopython
    blosum62 = MatrixInfo.blosum62  # Biopython charge la matrice blosum62
    
    # Paramètres
    k = 3
    max_distance = 40
    score_threshold = 15
    
    # Charger la base de données de séquences FASTA
    fasta_file = "subset_2000_sequences.fasta"
    database_sequences = load_fasta_database(fasta_file)
    
    len_database = 0
    for seq in database_sequences:
        len_database += len(seq)

    print(f"{len(database_sequences)} séquences chargées depuis la base de données FASTA.")


    # Pour chaque séquence dans la base de données, on effectue l'alignement avec la query
    for seq_target in database_sequences:
        # Extraction des k-mers
        kmers = extract_kmers(seq_query, k)
        
        # Indexation de la séquence cible
        target_index = index_target_sequence(seq_target, k)
        
        # Trouver les positions des k-mers
        kmer_positions = find_kmer_positions(kmers, target_index)
        
        # Trouver les double-hits
        double_hits = find_double_hits(kmer_positions, max_distance)
        
        # Filtrer et étendre les alignements
        alignments = filter_alignments(double_hits, seq_query, seq_target, blosum62, threshold=score_threshold)
        
        alignments = filter_duplicate_alignments(alignments)

        # Calculer les E-values
        e_values = calculate_e_values(alignments, seq_query, len_database)
        
        # Filtrer par E-value
        significant_alignments = filter_by_e_value(e_values, threshold=0.00001)
        
        # Afficher les alignements significatifs s'il y en a
        if significant_alignments:
            display_results(alignments, significant_alignments)

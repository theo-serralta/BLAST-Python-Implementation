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
                        # Double hit détecté
                        double_hits.append((kmer1, pos1, query_pos1, kmer2, pos2, query_pos2))
    return double_hits

# Fonction pour évaluer le score initial d'un double-hit
def evaluate_double_hit(kmer1, kmer2, blosum_matrix):
    score = 0
    for a, b in zip(kmer1, kmer2):
        score += get_blosum62_score(a, b, blosum_matrix)
    return score

# Extension d'un double-hit avec gapped alignment
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, gap_penalty=-4):
    score = 0
    alignment = []

    # Étendre vers la gauche
    left_score, left_alignment = extend_direction(seq_query, seq_target, query_pos, target_pos, blosum_matrix, gap_penalty, direction='left')
    
    # Étendre vers la droite
    right_score, right_alignment = extend_direction(seq_query, seq_target, query_pos, target_pos, blosum_matrix, gap_penalty, direction='right')
    
    # Combiner les scores et les alignements
    score = left_score + right_score
    alignment = left_alignment[::-1] + right_alignment  # Reverse the left alignment and combine with right
    
    return score, alignment

# Fonction pour étendre dans une direction (gauche ou droite)
def extend_direction(seq_query, seq_target, query_pos, target_pos, blosum_matrix, gap_penalty, direction='right', max_dropoff=10):
    score = 0
    alignment = []
    current_score = 0
    i = query_pos
    j = target_pos
    
    # Boucle d'extension
    while i >= 0 and j >= 0 and i < len(seq_query) and j < len(seq_target):
        if direction == 'left':
            i -= 1
            j -= 1
        else:  # direction 'right'
            i += 1
            j += 1
        
        # Stop si on dépasse les limites
        if i < 0 or j < 0 or i >= len(seq_query) or j >= len(seq_target):
            break
        
        # Score de l'alignement actuel
        match_score = get_blosum62_score(seq_query[i], seq_target[j], blosum_matrix)
        
        # Ajouter gap si le match est trop faible
        if match_score <= gap_penalty:
            alignment.append(('-', seq_target[j]))
            score += gap_penalty
        else:
            alignment.append((seq_query[i], seq_target[j]))
            score += match_score
        
        current_score += match_score
        
        # Si la différence entre le score actuel et le score maximum dépasse le max_dropoff, on arrête l'extension
        if current_score < -max_dropoff:
            break

    return score, alignment

# Filtrer les alignements par score
def filter_alignments(double_hits, seq_query, seq_target, blosum_matrix, threshold=15):
    alignments = []
    for hit in double_hits:
        kmer1, pos1, query_pos1, kmer2, pos2, query_pos2 = hit
        initial_score = evaluate_double_hit(kmer1, kmer2, blosum_matrix)
        if initial_score >= threshold:
            score, alignment = extend_alignment(seq_query, seq_target, query_pos1, pos1, blosum_matrix)
            alignments.append((score, alignment))
    return alignments

# Calcul de l'E-value à partir du score
def calculate_e_value(score, m, n, lambda_param=0.318, K=0.134):
    """
    Calcule l'E-value basée sur un score normalisé et la longueur des séquences query et target
    """
    e_value = K * m * n * math.exp(-lambda_param * score)
    return e_value

# Calcul des E-values pour tous les alignements
def calculate_e_values(alignments, seq_query, seq_target):
    e_values = []
    m = len(seq_query)
    n = len(seq_target)
    
    for score, alignment in alignments:
        e_value = calculate_e_value(score, m, n)
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
    max_distance = 10
    score_threshold = 15
    
    # Charger la base de données de séquences FASTA
    fasta_file = "subset_2000_sequences.fasta"
    database_sequences = load_fasta_database(fasta_file)
    
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
        
        # Calculer les E-values
        e_values = calculate_e_values(alignments, seq_query, seq_target)
        
        # Filtrer par E-value
        significant_alignments = filter_by_e_value(e_values, threshold=0.00001)
        
        # Afficher les alignements significatifs s'il y en a
        if significant_alignments:
            display_results(alignments, significant_alignments)

import math
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

# Fonction pour obtenir le score entre deux acides aminés avec BLOSUM62
def get_blosum62_score(a, b, blosum_matrix):
    """
    Retourne le score BLOSUM62 pour une paire d'acides aminés.

    Args:
        a (str): Premier acide aminé.
        b (str): Deuxième acide aminé.
        blosum_matrix (dict): Matrice BLOSUM62 chargée via Biopython.

    Returns:
        int: Le score de la paire d'acides aminés ou -4 si la paire n'existe pas.
    """
    if (a, b) in blosum_matrix:
        return blosum_matrix[(a, b)]
    elif (b, a) in blosum_matrix:
        return blosum_matrix[(b, a)]
    else:
        return -4  # pénalité par défaut pour mismatch/gap

# Extraction des k-mers
def extract_kmers(sequence, k):
    """
    Extrait les k-mers d'une séquence donnée.

    Args:
        sequence (str): La séquence à traiter.
        k (int): La longueur des k-mers.

    Returns:
        list: Liste de tuples contenant le k-mer et sa position.
    """
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append((kmer, i))  # Conserver l'indice d'origine
    return kmers

# Indexation de la séquence cible pour une recherche rapide des k-mers
def index_target_sequence(sequence, k):
    """
    Indexe une séquence cible en fonction des k-mers pour une recherche rapide.

    Args:
        sequence (str): Séquence cible à indexer.
        k (int): Longueur des k-mers.

    Returns:
        dict: Dictionnaire où les clés sont des k-mers et les valeurs sont les positions dans la séquence.
    """
    kmer_index = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmer_index:
            kmer_index[kmer] = []
        kmer_index[kmer].append(i)  # Conserver l'indice d'origine
    return kmer_index

# Trouver les positions des k-mers dans la séquence cible
def find_kmer_positions(kmers, target_index):
    """
    Trouve les positions des k-mers extraits dans une séquence cible indexée.

    Args:
        kmers (list): Liste de k-mers avec leurs positions.
        target_index (dict): Index de la séquence cible.

    Returns:
        dict: Dictionnaire des k-mers et leurs positions respectives dans les séquences query et target.
    """
    positions = {}
    for kmer, query_pos in kmers:
        if kmer in target_index:
            for target_pos in target_index[kmer]:
                if kmer not in positions:
                    positions[kmer] = []
                positions[kmer].append((target_pos, query_pos))  # position dans target et dans query
    return positions

# Fonction pour trouver les double-hits sur la même diagonale et sans superposition stricte
def find_double_hits(kmer_positions, max_distance):
    """
    Identifie les double-hits dans la séquence cible.

    Les double-hits sont des paires de k-mers sur la même diagonale, qui ne se superposent pas et
    dont la distance est inférieure à un seuil défini.

    Args:
        kmer_positions (dict): Positions des k-mers dans la séquence cible.
        max_distance (int): Distance maximale permise entre deux hits pour être considérés comme un double-hit.

    Returns:
        list: Liste de tuples représentant les double-hits détectés.
    """
    double_hits = []
    kmer_list = list(kmer_positions.items())
    
    for i in range(len(kmer_list)):
        kmer1, positions1 = kmer_list[i]
        for pos1, query_pos1 in positions1:
            for j in range(i + 1, len(kmer_list)):
                kmer2, positions2 = kmer_list[j]
                for pos2, query_pos2 in positions2:
                    # Vérification de la diagonale
                    if (pos1 - query_pos1) == (pos2 - query_pos2):
                        # Vérification de la non-superposition stricte et de la distance
                        if abs(pos2 - pos1) <= max_distance:
                            # Éviter uniquement les hits qui sont exactement au même endroit
                            if pos1 != pos2 or query_pos1 != query_pos2:
                                # Double hit détecté
                                double_hits.append((kmer1, pos1, query_pos1, kmer2, pos2, query_pos2))
    return double_hits

# Fonction pour évaluer le score initial d'un double-hit
def evaluate_double_hit(kmer1, kmer2, blosum_matrix):
    """
    Calcule le score initial d'un double-hit en se basant sur la matrice BLOSUM62.

    Args:
        kmer1 (str): Premier k-mer.
        kmer2 (str): Deuxième k-mer.
        blosum_matrix (dict): Matrice BLOSUM62.

    Returns:
        int: Le score du double-hit.
    """
    score = 0
    for a, b in zip(kmer1, kmer2):
        score += get_blosum62_score(a, b, blosum_matrix)
    return score

# Extension d'un double-hit avec gapped alignment
def extend_alignment(seq_query, seq_target, query_pos, target_pos, blosum_matrix, 
                     gap_open_penalty=-11, gap_extension_penalty=-2):
    """
    Étend un double-hit en un alignement avec prise en compte des gaps.

    Args:
        seq_query (str): Séquence query.
        seq_target (str): Séquence cible.
        query_pos (int): Position du hit dans la séquence query.
        target_pos (int): Position du hit dans la séquence cible.
        blosum_matrix (dict): Matrice BLOSUM62.
        gap_open_penalty (int): Pénalité pour l'ouverture d'un gap.
        gap_extension_penalty (int): Pénalité pour l'extension d'un gap.

    Returns:
        tuple: Le score de l'alignement et l'alignement lui-même sous forme de liste de tuples.
    """
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
    """
    Étend un alignement dans une direction donnée avec prise en compte des gaps.

    Args:
        seq_query (str): Séquence query.
        seq_target (str): Séquence cible.
        query_pos (int): Position dans la séquence query.
        target_pos (int): Position dans la séquence cible.
        blosum_matrix (dict): Matrice BLOSUM62.
        gap_open_penalty (int): Pénalité pour l'ouverture d'un gap.
        gap_extension_penalty (int): Pénalité pour l'extension d'un gap.
        direction (str): Direction de l'extension ('left' ou 'right').
        max_dropoff (int): Valeur seuil pour arrêter l'extension.

    Returns:
        tuple: Le score de l'extension et l'alignement sous forme de liste de tuples.
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
    """
    Filtre et étend les alignements basés sur les double-hits détectés.

    Args:
        double_hits (list): Liste des double-hits détectés.
        seq_query (str): Séquence query.
        seq_target (str): Séquence cible.
        blosum_matrix (dict): Matrice BLOSUM62.
        threshold (int): Seuil de score pour l'inclusion d'un alignement.

    Returns:
        list: Liste d'alignements filtrés.
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
    Filtre les alignements dupliqués en conservant uniquement les uniques.

    Args:
        alignments (list): Liste d'alignements.

    Returns:
        list: Liste d'alignements uniques.
    """
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
    Calcule l'E-value basée sur le score et la longueur des séquences.

    Args:
        score (int): Score de l'alignement.
        m (int): Longueur de la séquence query.
        n (int): Longueur totale de la base de données.
        lambda_param (float): Paramètre lambda pour la distribution des scores.
        K (float): Paramètre K pour la distribution des scores.

    Returns:
        float: L'E-value calculée.
    """
    e_value = K * m * n * math.exp(-lambda_param * score)
    return e_value

# Calcul des E-values pour tous les alignements
def calculate_e_values(alignments, seq_query, len_database):
    """
    Calcule les E-values pour une liste d'alignements.

    Args:
        alignments (list): Liste des alignements avec leur score.
        seq_query (str): Séquence query.
        len_database (int): Longueur totale de la base de données.

    Returns:
        list: Liste des E-values associées à chaque alignement.
    """
    e_values = []
    m = len(seq_query)
    n = len_database

    for score, alignment in alignments:
        e_value = calculate_e_value(score, m, n)
        e_values.append((e_value, score, alignment))
    
    return e_values

# Filtrer les alignements par E-value
def filter_by_e_value(e_values, threshold=0.01):
    """
    Filtre les alignements en fonction d'un seuil d'E-value.

    Args:
        e_values (list): Liste des alignements avec leurs E-values.
        threshold (float): Seuil d'E-value pour filtrer les alignements.

    Returns:
        list: Liste des alignements significatifs avec E-value inférieure au seuil.
    """
    filtered_alignments = [align for align in e_values if align[0] <= threshold]
    return filtered_alignments

# Fonction pour formater l'alignement pour un affichage BLAST-like
def format_alignment(seq_query, seq_target, alignment):
    """
    Formate un alignement pour un affichage de style BLAST.

    Args:
        seq_query (str): Séquence query.
        seq_target (str): Séquence cible.
        alignment (list): Alignement à formater.

    Returns:
        tuple: Chaînes formatées pour la query, la ligne de correspondance, et la cible.
    """
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
    """
    Affiche les résultats d'alignement de manière formatée.

    Args:
        alignments (list): Liste des alignements.
        e_values (list): Liste des E-values associées aux alignements.
    """
    print(f"{len(e_values)} alignements significatifs trouvés :\n")
    
    for i, (e_value, score, alignment) in enumerate(e_values, 1):
        query_str, match_str, target_str = format_alignment(seq_query, seq_target, alignment)

        print(f"Alignement {i}:")
        print(f"Score: {score}, E-value: {e_value:.2e}")
        print(f"Query:  {query_str}")
        print(f"        {match_str}")
        print(f"Target: {target_str}")
        print("\n" + "-"*50 + "\n")

# Fonction pour charger la base de données FASTA
def load_fasta_database(fasta_file):
    """
    Charge une base de données FASTA à partir d'un fichier.

    Args:
        fasta_file (str): Chemin vers le fichier FASTA.

    Returns:
        list: Liste des séquences de la base de données.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# Exemple d'utilisation avec la base de données FASTA
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

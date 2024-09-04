#Etape 1 : Recherche des mots dans les séquences protéiques

def extract_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers

def find_kmer_positions(kmers, target_sequence):
    positions = {}
    for kmer in kmers:
        positions[kmer] = []
        start = 0
        while True:
            pos = target_sequence.find(kmer, start)
            if pos == -1:  # Si le k-mer n'est plus trouvé, on arrête la recherche
                break
            positions[kmer].append(pos)
            start = pos + 1  # On continue la recherche après la position trouvée
    return positions

def find_double_hits(kmer_positions, max_distance):
    double_hits = []
    kmer_list = list(kmer_positions.items())
    
    for i in range(len(kmer_list)):
        kmer1, positions1 = kmer_list[i]
        for pos1 in positions1:
            for j in range(i+1, len(kmer_list)):
                kmer2, positions2 = kmer_list[j]
                for pos2 in positions2:
                    if abs(pos2 - pos1) <= max_distance:
                        double_hits.append((kmer1, pos1, kmer2, pos2))
    return double_hits


if __name__ == "__main__":
    # Séquence query
    query_sequence = "MKVYLGI"

    # Séquences cibles (par exemple dans une base de données)
    target_sequences = [
        "ATGMKVYLGIYQWER",
        "MKVYLLVYLGIMKVG"
    ]

    # Longueur du k-mer
    k = 3

    # Distance maximale entre les double-hits
    max_distance = 5

    # Extraire les k-mers de la séquence query
    kmers = extract_kmers(query_sequence, k)
    print(f"K-mers extraits de la séquence query: {kmers}")

    # Pour chaque séquence cible, trouver les positions des k-mers
    for i, target in enumerate(target_sequences):
        positions = find_kmer_positions(kmers, target)
        print(f"\nPositions des k-mers dans la séquence cible {i + 1}: {positions}")
        
        # Détecter les double-hits
        double_hits = find_double_hits(positions, max_distance)
        print(f"Double-hits détectés dans la séquence cible {i + 1}: {double_hits}")

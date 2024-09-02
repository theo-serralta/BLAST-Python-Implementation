###Partie code recherche MSP (Maximal segment pai)
def similarity_score(residue1, residue2, scoring_matrix):
    return scoring_matrix.get((residue1, residue2), 0)

def find_wmer_hits(sequence1, sequence2, w, scoring_matrix, threshold):
    hits = []
    for i in range(len(sequence1) - w + 1):
        wmer1 = sequence1[i:i+w]
        for j in range(len(sequence2) - w + 1):
            wmer2 = sequence2[j:j+w]
            score = sum(similarity_score(wmer1[k], wmer2[k], scoring_matrix) for k in range(w))
            if score >= threshold:
                hits.append((i, j, score))
    return hits

def extend_hit(sequence1, sequence2, start1, start2, scoring_matrix, dropoff):
    max_score = 0
    current_score = 0
    best_coords = (start1, start1, start2, start2)

    # Extend to the right
    i = start1
    j = start2
    while i < len(sequence1) and j < len(sequence2):
        current_score += similarity_score(sequence1[i], sequence2[j], scoring_matrix)
        if current_score > max_score:
            max_score = current_score
            best_coords = (start1, i+1, start2, j+1)
        elif max_score - current_score > dropoff:
            break
        i += 1
        j += 1

    return best_coords, max_score

def find_msp(sequence1, sequence2, scoring_matrix, w=4, threshold=10, dropoff=20):
    hits = find_wmer_hits(sequence1, sequence2, w, scoring_matrix, threshold)
    max_msp = None
    max_score = float('-inf')

    for hit in hits:
        start1, start2, hit_score = hit
        msp, msp_score = extend_hit(sequence1, sequence2, start1, start2, scoring_matrix, dropoff)
        if msp_score > max_score:
            max_score = msp_score
            max_msp = msp

    return max_msp, max_score




scoring_matrix = {
    ('A', 'A'): 5, ('A', 'C'): -4, ('A', 'G'): -4, ('A', 'T'): -4,
    ('C', 'A'): -4, ('C', 'C'): 5, ('C', 'G'): -4, ('C', 'T'): -4,
    ('G', 'A'): -4, ('G', 'C'): -4, ('G', 'G'): 5, ('G', 'T'): -4,
    ('T', 'A'): -4, ('T', 'C'): -4, ('T', 'G'):-4, ('T', 'T'):5
}

sequence1 = "ACTGA"
sequence2 = "GACTG"

msp, max_score = find_msp(sequence1, sequence2, scoring_matrix)
print(f"MSP: {msp}, Score: {max_score}")
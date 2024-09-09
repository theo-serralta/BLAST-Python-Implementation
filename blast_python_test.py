import unittest
import math
from Bio.SubsMat import MatrixInfo
from blast_python import get_blosum62_score, extract_kmers, index_target_sequence, find_kmer_positions
from blast_python import find_double_hits, evaluate_double_hit, extend_alignment, calculate_e_value_from_bitscore, filter_alignments, calculate_e_values, filter_by_e_value

class TestBlastAlgorithm(unittest.TestCase):

    def setUp(self):
        # Initialisation de la matrice BLOSUM62 pour les tests
        self.blosum62 = MatrixInfo.blosum62

    def test_get_blosum62_score(self):
        """ Test the function get_blosum62_score """
        # Test a perfect match
        self.assertEqual(get_blosum62_score('A', 'A', self.blosum62), 4)
        # Test a mismatch
        self.assertEqual(get_blosum62_score('A', 'W', self.blosum62), -3)
        # Test a pair not found in BLOSUM62 (default -4)
        self.assertEqual(get_blosum62_score('A', 'J', self.blosum62), -4)
        # Test a gap situation
        self.assertEqual(get_blosum62_score('-', 'A', self.blosum62), -4)  # Assuming gaps use -4 by default


    def test_extract_kmers(self):
        """ Test de l'extraction des k-mers """
        sequence = "MKTWQ"
        expected_kmers = [('MKT', 0), ('KTW', 1), ('TWQ', 2)]
        self.assertEqual(extract_kmers(sequence, 3), expected_kmers)

    def test_find_kmer_positions(self):
        """ Test de la recherche des positions de k-mers dans une séquence cible """
        sequence_query = "MKTWQ"
        sequence_target = "MKTWQKTWQ"
        kmers = extract_kmers(sequence_query, 3)
        target_index = index_target_sequence(sequence_target, 3)
        expected_positions = {
            'MKT': [(0, 0)],
            'KTW': [(1, 1), (5, 1)],
            'TWQ': [(2, 2), (6, 2)]
        }
        self.assertEqual(find_kmer_positions(kmers, target_index), expected_positions)

    def test_find_double_hits(self):
        """ Test de la détection des double-hits sur la même diagonale et sans superposition """
        kmer_positions = {
            'MKT': [(0, 0)],  # pos1=0, query_pos1=0
            'KTW': [(1, 1), (5, 1)],  # pos1=1, query_pos1=1 et pos2=5, query_pos2=1
            'TWQ': [(2, 2), (6, 2)]   # pos1=2, query_pos1=2 et pos2=6, query_pos2=2
        }
        
        # Les hits doivent être sur la même diagonale, donc on attend des correspondances diagonales :
        max_distance = 2
        expected_hits = [
            ('MKT', 0, 0, 'KTW', 1, 1),  # Même diagonale : pos1 - query_pos1 == pos2 - query_pos2
            ('KTW', 1, 1, 'TWQ', 2, 2),  # Même diagonale
            ('MKT', 0, 0, 'TWQ', 2, 2),  # Même diagonale
            ('KTW', 5, 1, 'TWQ', 6, 2)   # Même diagonale
        ]
        
        # Appel à la fonction et vérification des hits trouvés
        self.assertCountEqual(find_double_hits(kmer_positions, max_distance), expected_hits)

    def test_find_double_hits_not_same_diagonal(self):
        """ Test pour les hits qui ne sont pas sur la même diagonale """
        kmer_positions = {
            'MKT': [(0, 0)],  # pos1=0, query_pos1=0
            'KTW': [(2, 1)],  # pos1=2, query_pos1=1 (différence de diagonale pos1 - query_pos1 != pos2 - query_pos2)
            'TWQ': [(4, 2)]   # pos1=4, query_pos1=2
        }
        
        max_distance = 2
        expected_hits = []  # Aucune paire ne devrait être retenue car les diagonales ne correspondent pas
        
        self.assertCountEqual(find_double_hits(kmer_positions, max_distance), expected_hits)

    def test_find_double_hits_overlapping(self):
        """ Test pour les hits qui se superposent """
        kmer_positions = {
            'MKT': [(0, 0)],  # pos1=0, query_pos1=0
            'KTW': [(0, 0)],  # pos1=0, query_pos1=0 (se superpose complètement avec MKT)
            'TWQ': [(2, 2)]   # pos1=2, query_pos1=2 (ne se superpose pas avec MKT)
        }
        
        max_distance = 2
        expected_hits = [('MKT', 0, 0, 'TWQ', 2, 2), ('KTW', 0, 0, 'TWQ', 2, 2)]  
        
        # 'MKT' et 'KTW' ne devraient pas être trouvés car ils se superposent
        self.assertCountEqual(find_double_hits(kmer_positions, max_distance), expected_hits)

    def test_extension_of_partial_hits(self):
        """Test de l'extension des hits partiels"""
        seq_query = "MKTWQ"
        seq_target = "MKTGQ"
        
        double_hits = [('MKT', 0, 0, 'KTG', 1, 1)]  # Double hit avec mismatch
        alignments = filter_alignments(double_hits, seq_query, seq_target, self.blosum62, bit_score_threshold=-5)
        
        # S'assurer que l'alignement partiel est détecté malgré le mismatch
        self.assertGreater(len(alignments), 0, "Un alignement partiel devrait être trouvé même avec un mismatch.")


    def test_alignment_with_gaps(self):
        """
        Test a case where the alignment should include gaps.
        """
        seq_query = "KVIILTAAAQGIGQAAALAFAREGAK--VIATDINESKLQELEKYPGIQTRVLDVTKKKQIDQFANEVE----RLDVLFNVAGFVHHGTVLDCEEKDWDFSMNLNVRSMYLMIKAFLPKMLAQKSGNIINMSSVASSVKGVVNRCVYSTTKAAVIGLTKSVAADFIQQGIRCNCVCPGTVDTPSLQERIQARGNPEEARNDFLKRQKTGRFATAEEIAMLCVYLASDESAYVTGNPVIIDGGWSL"
        seq_target = "KVCAVFGGSRGIGRAVAQLMARKGYRLAVIARNLEGAKAAAGDLGGDHLAFSCDVAKEHDVQNTFEELEKHLGRVNFLVNAAGINRDGLLVRTKTEDMVSQLHTNLLGSMLTCKAAMRTMIQQQGGSIVNVGSIVGLKGNSGQSVYSASKGGLVGFSRALAKEVARKKIRVNVVAPGFVHTDMTKDLKEEHLKKNIPLGRFGETIEVAHAVVFLLESPYITGHVLVVDGGLQLIL"

        # Test with double-hit and alignment extension
        kmers = extract_kmers(seq_query, 3)
        target_index = index_target_sequence(seq_target, 3)
        kmer_positions = find_kmer_positions(kmers, target_index)
        double_hits = find_double_hits(kmer_positions, max_distance=10)

        # Extend alignment with gap penalties
        alignments = filter_alignments(double_hits, seq_query, seq_target, self.blosum62, bit_score_threshold=-30)

        # Expect an alignment with gaps
        self.assertTrue(len(alignments) > 0, "No alignments with gaps were found")


    def test_gapped_alignment(self):
        """
        Test to check that gaps are correctly inserted when needed.
        """
        seq_query = "KVIILTAAAQGIGQAAALAFAREGAKVIATDINESK"
        seq_target = "KVCAVFGGSRGIGRAVAQLMARKGYRLAVIARNLEGAK"

        # Expected alignment with gaps
        expected_alignment = [
            ('K', 'K'), ('V', 'V'), ('I', 'C'), ('I', 'A'), ('L', 'V'),
            ('T', 'F'), ('A', 'G'), ('A', 'G'), ('A', 'S'), ('Q', 'R'),
            ('G', 'G'), ('I', 'I'), ('G', 'G'), ('Q', 'R'), ('A', 'A'),
            ('A', 'V'), ('A', 'A'), ('L', 'Q'), ('A', 'L'), ('F', 'M'),
            ('A', 'A'), ('R', 'R'), ('E', 'K'), ('G', 'G'), ('A', 'Y'),
            ('K', 'R'), ('-', 'L'), ('-', 'A'), ('V', 'V'), ('I', 'I'),
            ('A', 'A'), ('T', 'R'), ('D', 'N'), ('I', 'L'), ('N', 'E'),
            ('E', 'G'), ('S', 'A'), ('K', 'K')
        ]

        score, alignment = extend_alignment(seq_query, seq_target, 0, 0, self.blosum62)
        #print(f"Generated alignment: {alignment}")
        #print(f"Expected alignment: {expected_alignment}")

        # Check that gaps are inserted correctly
        self.assertEqual(alignment, expected_alignment)

    def test_evaluate_double_hit(self):
        """ Test de la fonction evaluate_double_hit avec 3 cas : correspondance parfaite, mismatch et pénalité par défaut """
        
        # 1. Correspondance parfaite
        kmer1 = "MKT"
        kmer2 = "MKT"
        expected_score = (get_blosum62_score('M', 'M', self.blosum62) +
                        get_blosum62_score('K', 'K', self.blosum62) +
                        get_blosum62_score('T', 'T', self.blosum62))
        self.assertEqual(evaluate_double_hit(kmer1, kmer2, self.blosum62), expected_score, "Le score pour une correspondance parfaite est incorrect.")
        
        # 2. Mismatch (ici 'P' vs 'T')
        kmer1 = "MKT"
        kmer2 = "MKP"  # 'P' est différent de 'T'
        expected_score = (get_blosum62_score('M', 'M', self.blosum62) +
                        get_blosum62_score('K', 'K', self.blosum62) +
                        get_blosum62_score('T', 'P', self.blosum62))  # Pénalité pour mismatch ici
        self.assertEqual(evaluate_double_hit(kmer1, kmer2, self.blosum62), expected_score, "Le score pour un mismatch est incorrect.")

        # 3. Acide aminé inconnu dans BLOSUM62 (par exemple 'Z' qui n'est pas défini)
        kmer1 = "MKJ"  # J est absent de BLOSUM62
        kmer2 = "MKT"
        expected_score = (get_blosum62_score('M', 'M', self.blosum62) +
                        get_blosum62_score('K', 'K', self.blosum62) +
                        -4)  # Pénalité par défaut pour 'J'
        self.assertEqual(evaluate_double_hit(kmer1, kmer2, self.blosum62), expected_score, "Le score pour un acide aminé inconnu est incorrect.")

    def test_partial_alignment_with_gaps(self):
        """ Test a case where the alignment should include gaps """
        seq_query = "MKTW--Q"
        seq_target = "MKTWQQQ"
        
        score, alignment = extend_alignment(seq_query, seq_target, 0, 0, self.blosum62)
        
        # Expecting a score that takes into account the gap penalties
        expected_score = (
            get_blosum62_score('M', 'M', self.blosum62) +
            get_blosum62_score('K', 'K', self.blosum62) +
            get_blosum62_score('T', 'T', self.blosum62) +
            get_blosum62_score('W', 'W', self.blosum62) +
            -11 + -2 +  # Gap penalties
            get_blosum62_score('Q', 'Q', self.blosum62)
        )
        
        self.assertEqual(score, expected_score)
        self.assertEqual(alignment, [
            ('M', 'M'), ('K', 'K'), ('T', 'T'), ('W', 'W'),
            ('-', 'Q'), ('-', 'Q'), ('Q', 'Q')
        ])



    def test_extend_alignment(self):
        seq_query = "MKTWQ"
        seq_target = "MKTWQ"
        score, alignment = extend_alignment(seq_query, seq_target, 0, 0, self.blosum62)
        
        #print(f"Alignment: {alignment}, Score: {score}")

        # Vérifier que la longueur de l'alignement est correcte
        self.assertEqual(len(alignment), len(seq_query))  # Doit être 5
        
        # Vérifier que le score de l'alignement est correct
        self.assertEqual(score, 31)  # Le score total attendu est 26

    def calculate_e_value_from_bitscore(self):
        """ Test du calcul des E-values """
        score = 50
        m = 100  # Longueur de la séquence query
        n = 200  # Longueur de la séquence target
        lambda_param = 0.318
        K = 0.134
        expected_e_value = K * m * n * math.exp(-lambda_param * score)
        self.assertAlmostEqual(calculate_e_value_from_bitscore(score, m, n), expected_e_value, places=5)

    def test_filter_alignments(self):
        """ Test du filtrage des alignements par E-value """
        
        # Création d'alignements fictifs avec leurs scores
        alignments = [
            (100, [('MKT', 'MKT')]),  # Alignement avec un score élevé
            (80, [('KTW', 'KTW')]),   # Alignement avec un score moyen
            (50, [('TWQ', 'TWQ')]),   # Alignement avec un score plus faible
            (20, [('QW', 'QW')])      # Alignement avec un score très faible
        ]
        
        # On considère que la query et la target ont une longueur de 5
        seq_query = "MKTWQ"
        seq_target = "MKTWQ"
        len_database = len(seq_target)

        # Calculer les E-values pour ces alignements
        e_values = calculate_e_values(alignments, seq_query, len_database)
        
        # Afficher les E-values calculées pour vérification
        for e_value, score, alignment in e_values:
            print(f"Alignement: {alignment}, Score: {score}, E-value: {e_value}")
        
        # Test du filtrage des alignements avec un seuil d'E-value = 0.01
        filtered = filter_by_e_value(e_values, threshold=0.00001)
        print(f"Alignements retenus : {filtered}")
        
        # Assurer que des alignements ont été filtrés
        self.assertTrue(len(filtered) > 0, "Aucun alignement n'a été retenu après filtrage.")
        
        # Vérifier que seuls les alignements avec des E-values <= 0.01 sont présents
        for e_value, score, alignment in filtered:
            self.assertLessEqual(e_value, 0.01, f"E-value trop élevée pour {alignment} : {e_value}")
        
        # Test du filtrage avec un seuil plus élevé, par exemple 1.0 (pour s'assurer que tous les alignements sont inclus)
        filtered_higher_threshold = filter_by_e_value(e_values, threshold=1.0)
        self.assertEqual(len(filtered_higher_threshold), len(e_values), "Tous les alignements devraient être inclus avec un seuil de 1.0.")


    def test_integration(self):
        seq_query = "MKTWQ"
        seq_target = "MKTWQKTWQ"
        len_database = len(seq_target)

        # Exécution de l'algorithme complet
        kmers = extract_kmers(seq_query, 3)
        target_index = index_target_sequence(seq_target, 3)
        #print(f"K-mers: {kmers}")
        
        kmer_positions = find_kmer_positions(kmers, target_index)
        #print(f"K-mer positions: {kmer_positions}")
        
        double_hits = find_double_hits(kmer_positions, 2)
        #print(f"Double-hits: {double_hits}")
        
        alignments = filter_alignments(double_hits, seq_query, seq_target, self.blosum62, bit_score_threshold=-5)
        #print(f"Alignments: {alignments}")
        
        e_values = calculate_e_values(alignments, seq_query, len_database)
        #print(f"E-values: {e_values}")
        
        filtered_alignments = filter_by_e_value(e_values, threshold=0.01)
        #print(f"Filtered Alignments: {filtered_alignments}")
        
        # Vérification que l'alignement produit est correct et non vide
        self.assertTrue(len(filtered_alignments) > 0)


if __name__ == '__main__':
    unittest.main(argv=[''], verbosity=2, exit=False)

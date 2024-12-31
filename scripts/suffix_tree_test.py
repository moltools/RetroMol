class SuffixTreeNode:
    def __init__(self):
        self.children = {}
        self.indices = []  # Tracks indices of strings in which this node appears

class SuffixTree:
    def __init__(self):
        self.root = SuffixTreeNode()

    def add_kmers(self, sequence, index, k):
        """Add kmers of length `k` from a sequence to the suffix tree."""
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            current_node = self.root
            for item in kmer:
                if item not in current_node.children:
                    current_node.children[item] = SuffixTreeNode()
                current_node = current_node.children[item]
                current_node.indices.append(index)

    def find_kmers(self, kmer, scoring_matrix, threshold):
        """Find all sequence indices containing a sufficiently similar kmer based on scoring matrix."""
        current_node = self.root
        matching_indices = set()

        def is_match(node, remaining_kmer, current_score):
            if not remaining_kmer:
                if current_score / len(kmer) >= threshold:  # Normalize the score
                    matching_indices.update(node.indices)
                return

            item = remaining_kmer[0]
            for child_item, child_node in node.children.items():
                score = scoring_matrix.get(item, {}).get(child_item, 0)
                is_match(child_node, remaining_kmer[1:], current_score + score)

        is_match(current_node, kmer, 0)
        return list(matching_indices)

# Example sequences with fruit names (for testing), now as lists of strings
sequences = [
    ["apple", "banana", "cherry"],
    ["banana", "apple", "grape", "cherry"],
    ["grape", "apple", "cherry"],
    ["apple", "banana", "grape"],
    ["apple", "banana", "grape", "cherry"],
]

# Example scoring matrix for determining matches between fruits
fruit_scoring_matrix = {
    "apple": {"apple": 1.0, "banana": 0.8, "cherry": 0.5, "grape": 0.6},
    "banana": {"apple": 0.8, "banana": 1.0, "cherry": 0.4, "grape": 0.7},
    "cherry": {"apple": 0.5, "banana": 0.4, "cherry": 1.0, "grape": 0.6},
    "grape": {"apple": 0.6, "banana": 0.7, "cherry": 0.6, "grape": 1.0},
}

# Step 1: Construct a suffix tree for all kmers of sequences
k = 3  # Length of kmers
threshold = 0.8
suffix_tree = SuffixTree()
for i, sequence in enumerate(sequences):
    suffix_tree.add_kmers(sequence, i, k)

# Step 2: Query sequences using kmers and retrieve matches

def find_similar_sequences_kmers(suffix_tree, query_sequence, scoring_matrix, k, threshold=0.5):
    """Find sequences sufficiently similar to the query sequence based on kmers."""
    print(query_sequence)
    kmer_scores = {}
    for i in range(len(query_sequence) - k + 1):
        kmer = query_sequence[i:i + k]
        print(kmer)
        potential_matches = suffix_tree.find_kmers(kmer, scoring_matrix, threshold)
        print(potential_matches)
        for idx in potential_matches:
            if idx not in kmer_scores:
                kmer_scores[idx] = 0
            kmer_scores[idx] += 1  # Count matching k-mers
    print(kmer_scores)
    exit()

    similar_indices = [idx for idx, count in kmer_scores.items() if count > 0]
    return similar_indices

# Step 3: Create a comparison matrix based on kmer matching

def create_comparison_matrix_kmers(sequences, suffix_tree, scoring_matrix, k, threshold=0.5):
    """Generate a binary matrix indicating pairs of sequences to compare based on kmer scoring."""
    n = len(sequences)
    comparison_matrix = [[0] * n for _ in range(n)]

    for i, seq1 in enumerate(sequences):
        similar_indices = find_similar_sequences_kmers(suffix_tree, seq1, scoring_matrix, k, threshold)
        for j in similar_indices:
            if i != j:
                comparison_matrix[i][j] = 1  # Mark sequences to be compared

    return comparison_matrix

# Generate the comparison matrix
comparison_matrix = create_comparison_matrix_kmers(sequences, suffix_tree, fruit_scoring_matrix, k, threshold)

# Display the comparison matrix
for row in comparison_matrix:
    print(row)

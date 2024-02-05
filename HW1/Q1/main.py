### An implementation of Algorithms 1-5 from Richard Durbin's paper
### "Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT)"

from typing import List, Tuple

class GenomeStringGenerator:
    def __init__(self, sequences: List[str]):
        self.sequences = sequences
        self.num_sequences = len(sequences)
        self.sequence_length = len(sequences[0]) if sequences else 0
        self.c_k = self.calculate_ck()  # Calculate c_k based on sequences

    def calculate_ck(self) -> List[int]:
        # Calculate the count of '0's at each position k in all sequences.
        c_k = [0] * self.sequence_length
        for sequence in self.sequences:
            for k, nucleotide in enumerate(sequence):
                if nucleotide == '0':
                    c_k[k] += 1
        return c_k

# Algorithms 1 & 2
class PrefixAndDivergenceArraysBuilder:
    def __init__(self, genome_generator: GenomeStringGenerator):
        self.genome_generator = genome_generator

    def build_prefix_and_divergence_arrays(self, k: int, prefix_order: List[int], divergence_start: List[int]) -> Tuple[List[int], List[int]]:
        M = self.genome_generator.num_sequences
        a_next, b_next, divergence_next, divergence_extension = [], [], [], []
        for i in range(M):
            sequence_prefix = self.genome_generator.sequences[prefix_order[i]][k]  # Extract k-th value based on prefix_order
            next_divergence_value = divergence_start[i] if divergence_start[i] != k else k + 1

            if sequence_prefix == '0':
                a_next.append(prefix_order[i])
                divergence_next.append(next_divergence_value)
            else:
                b_next.append(prefix_order[i])
                divergence_extension.append(next_divergence_value)

        prefix_order_next = a_next + b_next  # Concatenate a_next and b_next to form prefix_order_next
        divergence_start_next = divergence_next + divergence_extension  # Concatenate divergence_next and divergence_extension

        print(f"Updated prefix_order at position {k+1}: {prefix_order_next}")
        print(f"Updated divergence_start at position {k+1}: {divergence_start_next}")
        return prefix_order_next, divergence_start_next

# Algorithm 3
class MatchReporter:
    MIN_MATCH_LENGTH = 2  # Minimum length for a match

    def __init__(self, genome_generator: GenomeStringGenerator):
        self.genome_generator = genome_generator

    def report_long_matches(self, prefix_order: List[int], divergence_start: List[int], k: int) -> None:
        M = self.genome_generator.num_sequences
        a_group, b_group = [], []
        for i in range(M):
            if divergence_start[i] <= k - self.MIN_MATCH_LENGTH:
                a_group.append(prefix_order[i]) if self.genome_generator.sequences[prefix_order[i]][k] == '0' else b_group.append(prefix_order[i])

        if len(a_group) > 0 and len(b_group) > 0:
            for index_a in a_group:
                for index_b in b_group:
                    print(f"Match from {index_a} to {index_b} ending at position {k}, longer than {self.MIN_MATCH_LENGTH}")

# Algorithm 4
class SetMaximalMatchReporter:
    def __init__(self, genome_generator: GenomeStringGenerator):
        self.genome_generator = genome_generator

    def report_set_maximal_matches(self, prefix_order: List[int], divergence_start: List[int]) -> None:
        N = self.genome_generator.sequence_length
        M = self.genome_generator.num_sequences
        # Extend divergence_start with sentinels to prevent out-of-range errors
        divergence_start_extended = [N + 1] + divergence_start + [N + 1]

        for k in range(N):
            # Process each sequence in the order determined for position k
            for i in range(1, M + 1):  # Adjust for the indexing with sentinels
                m, n = i - 1, i + 1
                while m > 0 and divergence_start_extended[m] <= divergence_start_extended[i]:
                    m -= 1
                while n <= M and divergence_start_extended[n] <= divergence_start_extended[i + 1]:
                    n += 1

                # Report matches considering adjusted indices for prefix_order
                for j in range(m + 1, i):
                    print(f"Match of {prefix_order[i-1]} to {prefix_order[j-1]} from {divergence_start_extended[i]} to {k}")
                for j in range(i + 1, n):
                    print(f"Match of {prefix_order[i-1]} to {prefix_order[j-1]} from {divergence_start_extended[i + 1]} to {k}")

# Algorithm 5
class ZMatchUpdater:
    def __init__(self, genome_generator: GenomeStringGenerator, z: str):
        self.genome_generator = genome_generator
        self.z = z
        self.e_k = [0] * (genome_generator.sequence_length + 1)  # Start of the longest match ending at k
        self.f_k = [0] * (genome_generator.sequence_length + 1)  # Start index of interval in a_k with the match
        self.g_k = [0] * (genome_generator.sequence_length + 1)  # End index of interval in a_k with the match

    def update_z_matches(self, k: int, prefix_order: List[int], divergence_start: List[int], sorted_sequences: List[str]):
        # Simplified calculation for f' and g' based on z[k]
        f_prime, g_prime = self.calculate_w(k, self.z[k])
        
        # Determine if the matches can be extended
        if f_prime < g_prime:
            e_prime = self.e_k[k]
        else:
            print(f"Reporting set-maximal matches ending at {k}")
            e_prime = divergence_start[f_prime] - 1  # Adjusted calculation for e'
            
            # Adjust f' and g' based on the new sequence's base at e_prime
            f_prime, g_prime, e_prime = self.adjust_boundaries(k, f_prime, g_prime, e_prime, sorted_sequences, divergence_start)

        # Update for the next position (k + 1)
        self.f_k[k + 1], self.g_k[k + 1], self.e_k[k + 1] = f_prime, g_prime, e_prime

    def calculate_w(self, k: int, z_value: str) -> Tuple[int, int]:
        # Logic to update f and g based on z[k]
        f_prime = 0 if z_value == '0' else self.genome_generator.c_k[k]
        g_prime = self.genome_generator.c_k[k] if z_value == '0' else self.genome_generator.num_sequences
        return f_prime, g_prime

    def adjust_boundaries(self, k: int, f_prime: int, g_prime: int, e_prime: int, sorted_sequences: List[str]) -> Tuple[int, int, int]:
        # Adjust f_prime and g_prime based on the base at e_prime and the sequences in X
        if e_prime >= 0 and f_prime > 0:
            if self.z[e_prime] == '0':
                # Assuming z matches better with sequences having '0' at e_prime
                f_prime -= 1
            else:
                # Assuming z matches better with sequences having '1' at e_prime
                f_prime += 1

            # Simplified logic to adjust e_prime; real logic would require comparing z to sequences in X
            while e_prime > 0 and self.z[e_prime - 1] == sorted_sequences[f_prime][e_prime - 1]:
                e_prime -= 1

        return f_prime, g_prime, e_prime

def main():
    sequences = ["001", "110", "100", "011"]  # Example haplotype sequences for bi-allelic data
    z = "010"
    print(f"Input sequences: {sequences}\n")
    genome_generator = GenomeStringGenerator(sequences)
    prefix_divergence_builder = PrefixAndDivergenceArraysBuilder(genome_generator)

    # Assume initial prefix_order and divergence_start for k=0
    initial_order = list(range(genome_generator.num_sequences))
    initial_divergence = [0] * genome_generator.num_sequences  # Assume all start positions for divergence are 0
    prefix_order_next, divergence_start_next = prefix_divergence_builder.build_prefix_and_divergence_arrays(0, initial_order, initial_divergence)

    # SetMaximalMatchReporter reports set maximal matches within the sequences
    set_maximal_match_reporter = SetMaximalMatchReporter(genome_generator)
    set_maximal_match_reporter.report_set_maximal_matches(prefix_order_next, divergence_start_next)

    z_match_updater = ZMatchUpdater(genome_generator, z)
    sorted_sequences = sorted(sequences)

    # Simulate c_k for demonstration; actual values would depend on the content of sequences in X
    z_match_updater.c_k = [2, 1, 1, 2]  # Mock counts of '0's at each position

    for k in range(len(z)):
        f_prime, g_prime = z_match_updater.calculate_w(k, z[k])
        f_prime, g_prime, e_prime = z_match_updater.adjust_boundaries(k, f_prime, g_prime, z_match_updater.e_k[k], sorted_sequences)
        print(f"At position {k}: f'={f_prime}, g'={g_prime}, e'={e_prime}")

if __name__ == "__main__":
    main()
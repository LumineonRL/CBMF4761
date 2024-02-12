import random
from typing import List
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

class GenomeGenerator:
    """A class to generate random genomes of specified length and introduce mutations."""
    
    def __init__(self, length: int) -> None:
        """Initialize the GenomeGenerator with a specified genome length."""
        self.length = length
        self.nucleotides = ['A', 'C', 'G', 'T']
    
    def generate_genome(self) -> str:
        """Generate a random genome of the specified length."""
        return ''.join(random.choices(self.nucleotides, k=self.length))
    
    def insert_mutation(self, genome: str, p: float = 0.001) -> str:
        """Insert a random nucleotide into the genome with probability p at each position."""
        mutated_genome = ""
        for nucleotide in genome:
            if random.random() < p:
                inserted_nucleotide = random.choice(self.nucleotides)
                mutated_genome += inserted_nucleotide
                print(f"Gap of length 1 introduced: inserted {inserted_nucleotide}")
            mutated_genome += nucleotide
        return mutated_genome

    def delete_mutation(self, genome: str, p: float = 0.001) -> str:
        """Delete a nucleotide from the genome with probability p at each position."""
        mutated_genome = ""
        for i, nucleotide in enumerate(genome):
            if random.random() < p:
                print("Deletion occurred.")
                continue  # Skip adding this nucleotide to mutated_genome
            mutated_genome += nucleotide
        return mutated_genome

    def substitution_mutation(self, genome: str, p: float = 0.001) -> str:
        """Substitute a nucleotide in the genome with another, with probability p at each position."""
        mutated_genome = ""
        for nucleotide in genome:
            if random.random() < p:
                new_nucleotide = random.choice([n for n in self.nucleotides if n != nucleotide])
                mutated_genome += new_nucleotide
                print(f"Substitution occurred: {nucleotide} -> {new_nucleotide}")
            else:
                mutated_genome += nucleotide
        return mutated_genome
    
    def insert_triplet_mutation(self, genome: str, p: float = 0.001) -> str:
        """Insert a triplet of the same nucleotide into the genome with probability p."""
        if random.random() < p:
            position = random.randint(0, len(genome))
            nucleotide = random.choice(self.nucleotides) * 3
            genome = genome[:position] + nucleotide + genome[position:]
            print(f"Gap of length 3 introduced: inserted {nucleotide}")
        return genome

    def delete_triplet_mutation(self, genome: str, p: float = 0.001) -> str:
        """Delete a triplet of the same nucleotide from the genome with probability p."""
        mutated_genome = genome
        for triplet in [nuc * 3 for nuc in self.nucleotides]:
            start_pos = mutated_genome.find(triplet)
            if start_pos != -1 and random.random() < p:
                mutated_genome = mutated_genome[:start_pos] + mutated_genome[start_pos + 3:]
                print(f"Deletion of triplet occurred: deleted {triplet}")
        return mutated_genome

    @staticmethod
    def display_genome(name: str, genome: str) -> None:
        """Display a genome with a given name, adjusting display based on genome length."""
        display_text = genome[:50] + "..." + genome[-50:] if len(genome) > 100 else genome
        print(f"{name}: {display_text} (Length: {len(genome)})")

class GenomeAlignment:
    def __init__(self, ref_genome: str, read_genome: str,
                 match_score: int, mismatch_penalty: int, gap_opening_penalty: int,
                gap_extension_penalty: int, triplet_opening_penalty: int, triplet_extension_penalty: int):
        self.ref_genome = ref_genome
        self.read_genome = read_genome
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_opening_penalty = gap_opening_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.triplet_opening_penalty = triplet_opening_penalty
        self.triplet_extension_penalty = triplet_extension_penalty
        self.alignment_matrix = [[0 for _ in range(len(read_genome) + 1)] for _ in range(len(ref_genome) + 1)]

    def gap_penalty(self, k: int) -> int:
        if k == 3:
            return self.triplet_opening_penalty
        elif k > 3:
            return self.triplet_opening_penalty + self.triplet_extension_penalty * (k - 3)
        else:
            return self.gap_opening_penalty + self.gap_extension_penalty * (k - 1)

    def align_genomes(self):
        # Initialize the first row and column of the matrix
        for i in range(1, len(self.ref_genome) + 1):
            self.alignment_matrix[i][0] = self.gap_penalty(i)
        for j in range(1, len(self.read_genome) + 1):
            self.alignment_matrix[0][j] = self.gap_penalty(j)
        
        # Fill in the rest of the alignment matrix
        for i in range(1, len(self.ref_genome) + 1):
            for j in range(1, len(self.read_genome) + 1):
                match_mismatch = self.alignment_matrix[i-1][j-1] + (self.match_score if self.ref_genome[i-1] == self.read_genome[j-1] else self.mismatch_penalty)
                gap_in_ref = max(self.alignment_matrix[i-k][j] - self.gap_penalty(k) for k in range(1, i+1))
                gap_in_read = max(self.alignment_matrix[i][j-k] - self.gap_penalty(k) for k in range(1, j+1))
                self.alignment_matrix[i][j] = max(match_mismatch, gap_in_ref, gap_in_read, 0)  # Ensure no negative scores

    
    def display_alignment_matrix(self) -> None:
        """Display the alignment matrix as a heatmap with genome characters as axis labels and cell values."""
        plt.figure(figsize=(10, 8))
        
        # Convert the alignment matrix to a NumPy array for advanced slicing
        alignment_matrix_np = np.array(self.alignment_matrix)
        
        # Use NumPy slicing to exclude the first row and column
        heatmap_data = alignment_matrix_np[1:, 1:]  # Now works thanks to NumPy
        
        # Use genome characters for x and y axes labels
        x_labels = list(self.read_genome)
        y_labels = list(self.ref_genome)
        
        # Plot heatmap with annotations and color bar for the adjusted data
        sns.heatmap(heatmap_data, annot=True, cmap="viridis",
                    xticklabels=x_labels, yticklabels=y_labels, cbar=True, fmt="d")
        
        plt.title("Genome Alignment Matrix")
        plt.xlabel("Read Genome Characters")
        plt.ylabel("Reference Genome Characters")
        
        # Rotate the x-axis labels for better visibility
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        
        plt.show()


def main() -> None:
    SEED = 17

    GENOME_LENGTH = 15
    MUTATION_PROBABILITY = 0.1

    MATCH_SCORE = -10
    MISMATCH_PENALTY = 5
    GAP_OPENING_PENALTY = 14
    GAP_EXTENSION_PENALTY = 8
    TRIPLET_OPENING_PENALTY = 20
    TRIPLET_EXTENSION_PENALTY = 12

    random.seed(SEED)

    genome_generator = GenomeGenerator(GENOME_LENGTH)
    
    # Generate the reference genome and introduce mutations to create the read genome
    genome_reference = genome_generator.generate_genome()
    genome_read = genome_generator.insert_mutation(genome_reference, MUTATION_PROBABILITY)
    genome_read = genome_generator.delete_mutation(genome_read, MUTATION_PROBABILITY)
    genome_read = genome_generator.substitution_mutation(genome_read, MUTATION_PROBABILITY)
    genome_read = genome_generator.insert_triplet_mutation(genome_read, MUTATION_PROBABILITY)
    genome_read = genome_generator.delete_triplet_mutation(genome_read, MUTATION_PROBABILITY)
    
    # Perform genome alignment
    genome_alignment = GenomeAlignment(genome_reference, genome_read,
                                       match_score=MATCH_SCORE,
                                       mismatch_penalty=MISMATCH_PENALTY,
                                       gap_opening_penalty=GAP_OPENING_PENALTY,
                                       gap_extension_penalty=GAP_EXTENSION_PENALTY,
                                       triplet_opening_penalty=TRIPLET_OPENING_PENALTY,
                                       triplet_extension_penalty=TRIPLET_EXTENSION_PENALTY)
    genome_alignment.align_genomes()
    
    # Display the generated genomes and the alignment matrix
    GenomeGenerator.display_genome("Genome Reference", genome_reference)
    GenomeGenerator.display_genome("Mutated Genome Read", genome_read)
    print("\nAlignment Matrix:")

    genome_alignment.display_alignment_matrix()

if __name__ == "__main__":
    main()

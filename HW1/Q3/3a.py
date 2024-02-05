import numpy as np

# Set seed for reproducibility
np.random.seed(42)

def simulate_genome(length=250000):
    """Simulate a random genome sequence of a given length."""
    return ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=length))

def introduce_mutations(genome, mutation_rate=1/250000):
    """Introduce mutations into the genome."""
    mutated_genome = ''
    for char in genome:
        if np.random.rand() < mutation_rate:
            mutated_genome += np.random.choice([c for c in 'ACGT' if c != char])
        else:
            mutated_genome += char
    return mutated_genome

def simulate_misreading(genome, error_rate=0.008):
    """Simulate misreading the genome with a given error rate."""
    misread_genome = ''
    for char in genome:
        if np.random.rand() < error_rate:
            misread_genome += np.random.choice([c for c in 'ACGT' if c != char])
        else:
            misread_genome += char
    return misread_genome

def compare_stretches(original_genome, read_genome, stretch_length=50):
    """Compare 50-letter stretches between the original and read genomes."""
    differences = 0
    for i in range(0, len(original_genome), stretch_length):
        original_stretch = original_genome[i:i+stretch_length]
        read_stretch = read_genome[i:i+stretch_length]
        if original_stretch != read_stretch:
            differences += 1
    return differences

# Simulate the genome and introduce errors
GENOME_LENGTH = 25000000
MUTATION_RATE = 1/250000
ERROR_RATE = 0.008
STRETCH_LENGTH = 50

original_genome = simulate_genome(GENOME_LENGTH)
mutated_genome = introduce_mutations(original_genome, MUTATION_RATE)
read_genome = simulate_misreading(mutated_genome, ERROR_RATE)

# Compare the stretches
segments_evaluated = GENOME_LENGTH / STRETCH_LENGTH
differences = compare_stretches(original_genome, read_genome, STRETCH_LENGTH)
error_free_percentage = 1 - (differences / segments_evaluated)
print(f"Number of differing 50-letter stretches: {differences}")
print(f"Percentage of error-free 50-letter stretches: {error_free_percentage}")
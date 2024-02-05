from collections import Counter
import random

def burrows_wheeler_transform(s: str) -> str:
    """Apply Burrows-Wheeler Transform to a given string."""
    s += '$'  # Append EOS character
    rotations = [s[i:] + s[:i] for i in range(len(s))]
    sorted_rotations = sorted(rotations)
    bwt = ''.join(rotation[-1] for rotation in sorted_rotations)
    return bwt

def simulate_genome_with_cag_repeats(p: int, q: int, genome_size: int) -> str:
    """Simulate a genome sequence of a given size with specified CAG repeats."""
    # Ensure p > q + 1 and both are positive
    assert p > q + 1 and p > 0 and q > 0
        
    # Create the repeat sequences
    repeat_p = "CAG" * p
    repeat_q = "CAG" * q

    repeats = [repeat_p, repeat_q]
    
    # Fill the rest of the genome with random nucleotides to reach the desired size
    NUCLEOTIDES = ['A', 'C', 'G', 'T']
    remaining_size = genome_size - sum(len(repeat) for repeat in repeats)
    if remaining_size < 0:
        raise ValueError("Genome size too small to fit specified CAG repeats.")
    
    filler = ''.join(random.choices(NUCLEOTIDES, k=remaining_size))
    
    # Construct the genome with the two CAG repeats and filler
    insert_positions = sorted(random.sample(range(remaining_size + 1), 2))
    genome = filler[:insert_positions[0]] + repeats[0] + filler[insert_positions[0]:insert_positions[1]] + repeats[1] + filler[insert_positions[1]:]
    
    return genome

def analyze_bwt_of_genome_with_cag_repeats(p: int, q: int, genome_size: int):
    genome = simulate_genome_with_cag_repeats(p, q, genome_size)
    bwt = burrows_wheeler_transform(genome)
    
    # Find the longest run of 'C's in the BWT'd genome
    longest_run_length = 0
    longest_run_start = -1
    current_run_length = 0
    for i, char in enumerate(bwt):
        if char == 'C':
            current_run_length += 1
            if current_run_length > longest_run_length:
                longest_run_length = current_run_length
                longest_run_start = i - longest_run_length + 1
        else:
            current_run_length = 0
    
    nucleotide_count = Counter(bwt[:3*(p-q)-1])
    
    print(f"genome: {genome}")
    print(f"BWT'd genome: {bwt}")
    print(f"First 3(p-q)-1 nucleotides of (CAG)^p within the BWT: {bwt[:3*(p-q)-1]}")
    print(f"Nucleotide distribution in this segment: {nucleotide_count}")
    print(f"Longest run of 'C's in the BWT'd genome starts at position {longest_run_start} with length {longest_run_length}.")

# Example parameters
p = 5  # Length of the longest CAG repeat
q = 3  # Length of the second longest CAG repeat
genome_size = 100  # Total size of the simulated genome

analyze_bwt_of_genome_with_cag_repeats(p, q, genome_size)

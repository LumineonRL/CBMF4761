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

def calculate_occurrences(bwt):
    """Calculate the number of occurrences of each character in BWT up to each position."""
    occ = {'$': [0], 'A': [0], 'C': [0], 'G': [0], 'T': [0]}
    for char in bwt:
        for key in occ:
            occ[key].append(occ[key][-1] + (1 if char == key else 0))
    return occ

def get_cumulative_count(bwt):
    """Calculate the starting index of each character in the sorted BWT (first column)."""
    count = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for char in bwt:
        for c in count.keys():
            count[c] += (char < c)
    return count

def fm_mapping(bwt, pattern, occ, C):
    """Perform FM-mapping to check if a pattern exists in the BWT."""
    # Start with the last character of the pattern
    char = pattern[-1]
    # Get the range of rows in BWT starting with the last character of the pattern
    if char in C:
        l, r = C[char], C[char] + occ[char][-1] - 1
    else:
        return False  # Character not in BWT, pattern doesn't exist
    
    for char in reversed(pattern[:-1]):
        # Update the range using Occurrence and C
        l = C[char] + occ[char][l]
        r = C[char] + occ[char][r + 1] - 1
        if l > r:
            return False  # Pattern does not exist
    
    return True  # Pattern exists

def main():
    p, q, genome_size = 5, 3, 1000
    pattern = "CAGCAGCAGCAGCAG"
    
    genome = simulate_genome_with_cag_repeats(p, q, genome_size)
    bwt = burrows_wheeler_transform(genome)
    occ = calculate_occurrences(bwt)
    C = get_cumulative_count(bwt)

    pattern_exists = fm_mapping(bwt, pattern, occ, C)
    
    print(f"Pattern '{pattern}' exists in the genome: {pattern_exists}")

if __name__ == "__main__":
    main()
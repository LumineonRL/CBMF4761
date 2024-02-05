import random

def nucleotide_to_bits(nucleotide):
    """Convert a nucleotide to its 2-bit representation."""
    return {
        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11',
    }[nucleotide]

def read_to_bits(read):
    """Convert a DNA read into a bit-packed representation."""
    return ''.join(nucleotide_to_bits(n) for n in read)

def generate_random_reads(num_reads, read_length=50):
    """Generate a specified number of random reads of a given length."""
    NUCLEOTIDES = ['A', 'C', 'G', 'T']
    return [''.join(random.choices(NUCLEOTIDES, k=read_length)) for _ in range(num_reads)]

def store_reads_in_hashmap(reads):
    """Store bit-packed reads in a hashmap."""
    read_hashmap = {}
    for read in reads:
        read_bits = read_to_bits(read)
        read_hashmap[read_bits] = True
    return read_hashmap

# Simulating storage constraint
max_reads_within_1MB = 838860
reads = generate_random_reads(max_reads_within_1MB)
read_hashmap = store_reads_in_hashmap(reads)

print(f"Number of reads stored: {len(read_hashmap)}")

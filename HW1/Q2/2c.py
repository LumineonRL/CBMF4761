import bz2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def generate_genome_string(length):
    """Generate a random genome string of specified length."""
    return ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=length))

def compress_with_bz2(input_string):
    """Compress a string using bz2 and return the compressed data."""
    return bz2.compress(input_string.encode())

def calculate_bits(input_data):
    """Calculate the number of bits of a bytes object."""
    return len(input_data) * 8

def main():
    genome_sizes = np.logspace(1, 8, num=10, base=10).astype(int)  # From 10 to 100,000,000
    original_bits_list = []
    compressed_bits_list = []

    for size in genome_sizes:
        genome = generate_genome_string(size)
        original_bits = calculate_bits(genome.encode())  # Convert string to bytes for bit calculation
        compressed_genome = compress_with_bz2(genome)
        compressed_bits = calculate_bits(compressed_genome)
        
        original_bits_list.append(original_bits)
        compressed_bits_list.append(compressed_bits)

    # Plotting the results
    plt.figure(figsize=(10, 6))
    sns.lineplot(x=genome_sizes, y=original_bits_list, label='Original Bits', marker='o', color='blue')
    sns.lineplot(x=genome_sizes, y=compressed_bits_list, label='Compressed Bits', marker='o', color='red')
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Genome Size')
    plt.ylabel('Number of Bits (log(n), base 10)')
    plt.title('Compression Performance using bz2')
    plt.legend()
    plt.grid(True)
    plt.show()

main()

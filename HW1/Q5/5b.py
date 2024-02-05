import numpy as np

def simulate_read_with_minimizer_impact(read_length: int, kmer_length: int, error_rate: float) -> bool:
    critical_minimizer_start = np.random.randint(0, read_length - kmer_length + 1)
    errors = np.random.rand(read_length) < error_rate
    minimizer_affected = np.any(errors[critical_minimizer_start:critical_minimizer_start+kmer_length])
    return minimizer_affected

def monte_carlo_simulation_adjusted(num_reads: int, read_length: int, kmer_length: int, error_rate: float) -> float:
    missed_reads = sum(simulate_read_with_minimizer_impact(read_length, kmer_length, error_rate) for _ in range(num_reads))
    return missed_reads / num_reads

def find_optimal_kmer_length(num_reads: int, read_length: int, initial_kmer_length: int, error_rate: float, target_miss_rate: float) -> int:
    kmer_length = initial_kmer_length
    while True:
        fraction_missed = monte_carlo_simulation_adjusted(num_reads, read_length, kmer_length, error_rate)
        print(f"Trying k-mer length {kmer_length}, fraction missed: {fraction_missed:.4f}")
        
        if fraction_missed < target_miss_rate:
            break
        else:
            kmer_length -= 1
        
        if kmer_length == 0:
            raise ValueError("k-mer length reduced to zero without achieving target miss rate.")
    
    return kmer_length

NUM_READS = 10000
READ_LENGTH = 330
INITIAL_KMER_LENGTH = 31
ERROR_RATE = 1 / 31
TARGET_MISS_RATE = 0.4

optimal_kmer_length = find_optimal_kmer_length(NUM_READS, READ_LENGTH, INITIAL_KMER_LENGTH, ERROR_RATE, TARGET_MISS_RATE)
print(f"Optimal k-mer length to achieve a miss rate below {TARGET_MISS_RATE}: {optimal_kmer_length}")

import numpy as np

def simulate_read_with_minimizer_impact(read_length: int, kmer_length: int, error_rate: float) -> bool:
    """
    Simulate a single read and determine if an error impacts the selected minimizer.
    """
    # Assume the position of a critical minimizer within the read
    critical_minimizer_start = np.random.randint(0, read_length - kmer_length + 1)
    
    # Simulate errors across the read
    errors = np.random.rand(read_length) < error_rate
    
    # Check if any error occurs within the critical minimizer region
    minimizer_affected = np.any(errors[critical_minimizer_start:critical_minimizer_start+kmer_length])
    
    return minimizer_affected

def monte_carlo_simulation_adjusted(num_reads: int, read_length: int, kmer_length: int, error_rate: float) -> float:
    """
    Perform a Monte Carlo simulation to estimate the fraction of reads missed due to errors impacting minimizers.
    """
    missed_reads = sum(simulate_read_with_minimizer_impact(read_length, kmer_length, error_rate) for _ in range(num_reads))
    return missed_reads / num_reads

# Adjusted simulation parameters
np.random.seed(42)
NUM_READS = 10000000
READ_LENGTH = 330
KMER_LENGTH = 31
ERROR_RATE = 1 / 31

# Perform the adjusted Monte Carlo simulation
fraction_missed_adjusted = monte_carlo_simulation_adjusted(NUM_READS, READ_LENGTH, KMER_LENGTH, ERROR_RATE)

print(f"Adjusted fraction of reads missed: {fraction_missed_adjusted:.4f}")

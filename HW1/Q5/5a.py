import numpy as np

def simulate_read_with_errors(read_length: int, kmer_length: int, error_rate: float) -> bool:
    """
    Simulate a single read and determine if it contains an erroneous k-mer.
    :param read_length: Length of the read.
    :param kmer_length: Length of the k-mers.
    :param error_rate: Probability of an error in a single base.
    :return: True if the read contains at least one erroneous k-mer, False otherwise.
    """
    # Total number of k-mers in the read
    num_kmers = read_length - kmer_length + 1
    
    # Probability of a k-mer being error-free
    prob_error_free_kmer = (1 - error_rate) ** kmer_length
    
    # Using binomial distribution to determine if at least one k-mer is erroneous
    # If no k-mers are error-free, the read is considered missed
    num_error_free_kmers = np.random.binomial(n=num_kmers, p=prob_error_free_kmer)
    
    return num_error_free_kmers == 0

def monte_carlo_simulation(num_reads: int, read_length: int, kmer_length: int, error_rate: float) -> float:
    """
    Perform a Monte Carlo simulation to estimate the fraction of reads missed.
    :param num_reads: Number of reads to simulate.
    :param read_length: Length of each read.
    :param kmer_length: Length of the k-mers.
    :param error_rate: Probability of an error in a single base.
    :return: Fraction of reads missed due to errors.
    """
    missed_reads = sum(simulate_read_with_errors(read_length, kmer_length, error_rate) for _ in range(num_reads))
    return missed_reads / num_reads

# Parameters
np.random.seed(42)
NUM_READS = 1000000
READ_LENGTH = 330  # Length of the reads
KMER_LENGTH = 31   # Length of the k-mers
ERROR_RATE = 1 / 31  # Average error rate

# Perform the Monte Carlo simulation
fraction_missed = monte_carlo_simulation(NUM_READS, READ_LENGTH, KMER_LENGTH, ERROR_RATE)

print(f"Fraction of reads missed: {fraction_missed:.9f}")


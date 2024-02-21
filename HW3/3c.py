import numpy as np
from scipy.stats import norm, binomtest
import matplotlib.pyplot as plt

def simulate_transcript_reads(mean_fragment_length=250, std_fragment_length=25, read_length=60, target_probability=0.99, total_simulations=10000):
    # Definitions of transcript lengths according to the example
    exon1_length = 900  # Fixed length for exon 1
    exon2_short = 200   # Short version of exon 2
    exon2_long = 300    # Long version of exon 2
    exon3_short = 200   # Short version of exon 3
    exon3_long = 300    # Long version of exon 3
    poly_a_tail = 10    # Poly-A tail length
    
    # Dictionary to hold the full length of each transcript type
    transcripts_lengths = {
        'short_short': exon1_length + exon2_short + exon3_short + poly_a_tail,
        'short_long': exon1_length + exon2_short + exon3_long + poly_a_tail,
        'long_short': exon1_length + exon2_long + exon3_short + poly_a_tail,
        'long_long': exon1_length + exon2_long + exon3_long + poly_a_tail,
    }
    
    # The number of fully observed transcripts after simulations
    fully_observed_counts = {key: 0 for key in transcripts_lengths.keys()}
    
    # For calculating the simulation's required parameters
    mean_detected_transcripts = 0
    success_criteria = target_probability
    mean_reads_to_success = 0
    
    # Loop over simulations
    for simulation in range(total_simulations):
        observed_transcripts = {key: False for key in transcripts_lengths.keys()}  # Reset for each simulation
        reads_counter = 0
        while not all(observed_transcripts.values()):
            # 1. Choose a random transcript according to their proportion in the sample
            chosen_transcript = np.random.choice(list(transcripts_lengths.keys()))
            # 2. Sample a normal distributed fragment length
            fragment_length = int(round(np.random.normal(mean_fragment_length, std_fragment_length)))
            # 3. Attempt to simulate reads covering exons' switches (unique identifiers for each transcript)
            if fragment_length >= read_length:  # Check if the fragment is large enough for reads from both ends
                read_from_start = read_length
                read_from_end = transcripts_lengths[chosen_transcript] - read_length  # Get a read from the other end
                if chosen_transcript == 'short_short' and read_from_end <= exon1_length + exon2_short + exon3_short:
                    observed_transcripts['short_short'] = True
                elif chosen_transcript == 'short_long' and read_from_end <= exon1_length + exon2_short + exon3_long:
                    observed_transcripts['short_long'] = True
                elif chosen_transcript == 'long_short' and read_from_end <= exon1_length + exon2_long + exon3_short:
                    observed_transcripts['long_short'] = True
                elif chosen_transcript == 'long_long' and read_from_end <= exon1_length + exon2_long + exon3_long:
                    observed_transcripts['long_long'] = True
            # 4. Counter for reads to ensure we capture how many needed to cover all types
            reads_counter += 1
            if all(observed_transcripts.values()):  # If all detected at this stage, track and exit while loop
                mean_reads_to_success += reads_counter
                break
                
        # Loop to see how successful the detection was
        for t in observed_transcripts.keys():
            if observed_transcripts[t]:
                fully_observed_counts[t] += 1
    
    # Print how many times out of the simulations all were fully observed
    for transcript, observed_count in fully_observed_counts.items():
        mean_detected_transcripts = observed_count / total_simulations
        print(f"Transcript {transcript} was observed in {observed_count}/{total_simulations} simulations. (Mean detected probability: {mean_detected_transcripts})")
    
    # Average number of reads needed to meet the goal
    average_reads_needed = mean_reads_to_success / total_simulations
    
    return average_reads_needed

# Running the function
average_reads_needed = simulate_transcript_reads()
print(f"Average reads needed to observe all 4 transcripts with >99% probability: {average_reads_needed}")

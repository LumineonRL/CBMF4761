import pandas as pd
import numpy as np
from scipy.stats import geom
from typing import Tuple

class StopCodonAnalysis:
    @staticmethod
    def calculate_p_stop_codon(q: float) -> float:
        p_tag_tga = 2 * (1/2 - q)**2 * q  # Combined probability for TAG and TGA
        p_taa = (1/2 - q)**3  # Probability for TAA
        return p_tag_tga + p_taa

    @staticmethod
    def find_additional_trinucleotides(q: float, confidence_level: float) -> Tuple[int, float]:
        p_stop_codon = StopCodonAnalysis.calculate_p_stop_codon(q)
        target_p_value = 1 - confidence_level

        n_initial = 0
        while True:
            n_initial += 1
            cumulative_p_value = geom.sf(n_initial, p_stop_codon)
            if cumulative_p_value <= target_p_value:
                break

        return n_initial, cumulative_p_value

def main():
    q_values = np.arange(0, 0.50, 0.01)  # Q values from 0 to 0.5 with a step of 0.01
    confidence_level = 0.995  # Defined in main as per instruction

    results = []
    for q in q_values:
        n_additional, final_p_value = StopCodonAnalysis.find_additional_trinucleotides(q, confidence_level)
        results.append({
            'q': q,
            'Confidence Level': confidence_level,
            'N (Additional Trinucleotides)': n_additional,
            'Final p-value': final_p_value
        })
    
    results_df = pd.DataFrame(results)
    csv_path = 'stop_codon_analysis_results.csv'
    results_df.to_csv(csv_path, index=False)
    print(f"Results saved to '{csv_path}'.")

main()
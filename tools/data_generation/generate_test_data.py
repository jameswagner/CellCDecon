#!/usr/bin/env python3
"""
Generate synthetic methylation data for CellCDecon testing.

This script creates a test file with known cell type compositions and 
methylation values for each probe/tissue combination.
"""

import numpy as np
import pandas as pd
import argparse
import os

def generate_methylation_data(
    num_samples=50,
    num_probes=200,
    num_cell_types=3,
    output_file="synthetic_methylation.txt",
    add_noise=True,
    noise_level=0.03,
    generate_truth_files=True,
    concentration_param=0.5  # Lower values create more biased cell compositions
):
    """
    Generate synthetic methylation data with known cell type compositions.
    
    Parameters:
    - num_samples: Number of samples (individuals)
    - num_probes: Number of CpG probes to generate
    - num_cell_types: Number of distinct cell types
    - output_file: Path to save the generated methylation matrix
    - add_noise: Whether to add random noise to the pure signals
    - noise_level: Standard deviation of Gaussian noise
    - generate_truth_files: Whether to output true compositions and means
    - concentration_param: Dirichlet concentration parameter (lower = more biased toward one cell type)
    """
    print(f"Generating synthetic methylation data with {num_samples} samples, "
          f"{num_probes} probes, and {num_cell_types} cell types...")
    
    # Generate probe IDs
    probe_ids = [f"cg{str(i).zfill(6)}" for i in range(num_probes)]
    
    # Generate sample IDs
    sample_ids = [f"Sample_{i+1}" for i in range(num_samples)]
    
    # Generate true cell type compositions for each sample (sum to 1)
    # Using Dirichlet distribution with LOW alpha to generate compositions with clear dominance
    true_compositions = np.random.dirichlet(
        np.ones(num_cell_types) * concentration_param, size=num_samples
    )
    
    # Create some samples with extreme bias (>80% one cell type)
    extreme_samples = int(num_samples * 0.3)  # 30% of samples will be extreme
    for i in range(extreme_samples):
        dominant_type = i % num_cell_types
        true_compositions[i] = np.zeros(num_cell_types)
        true_compositions[i][dominant_type] = 0.8 + np.random.random() * 0.15  # 80-95% dominant
        
        # Distribute remaining proportion randomly
        remaining = 1.0 - true_compositions[i][dominant_type]
        other_types = list(range(num_cell_types))
        other_types.remove(dominant_type)
        
        for j in range(len(other_types) - 1):
            prop = remaining * np.random.random()
            true_compositions[i][other_types[j]] = prop
            remaining -= prop
        
        true_compositions[i][other_types[-1]] = remaining
    
    # Generate true methylation values for each cell type and probe
    # We'll create very distinct patterns for different cell types
    true_methylation = np.zeros((num_cell_types, num_probes))
    
    # Divide probes into categories:
    # 1. Cell-type specific (high in one, low in others)
    # 2. Gradually varying across cell types
    # 3. Random variation (but still distinct)
    
    num_type_specific = int(num_probes * 0.4)  # 40% of probes are cell-type specific
    num_gradual = int(num_probes * 0.3)  # 30% show gradual variation
    num_random = num_probes - num_type_specific - num_gradual  # The rest are random
    
    # Cell-type specific probes
    for cell_type in range(num_cell_types):
        # Assign probes specific to this cell type
        probes_per_type = num_type_specific // num_cell_types
        start_idx = cell_type * probes_per_type
        end_idx = start_idx + probes_per_type
        
        # Set methylation values: high for this cell type, low for others
        for probe_idx in range(start_idx, end_idx):
            for ct in range(num_cell_types):
                if ct == cell_type:
                    # High methylation (0.7-0.95) for this cell type
                    true_methylation[ct, probe_idx] = 0.7 + np.random.random() * 0.25
                else:
                    # Low methylation (0.05-0.3) for other cell types
                    true_methylation[ct, probe_idx] = 0.05 + np.random.random() * 0.25
    
    # Gradual variation probes
    start_idx = num_type_specific
    end_idx = start_idx + num_gradual
    
    for probe_idx in range(start_idx, end_idx):
        # Create a gradient effect across cell types
        base_value = np.random.random() * 0.3 + 0.1  # Random base (0.1-0.4)
        step = (0.9 - base_value) / (num_cell_types - 1)  # Step to reach ~0.9
        
        for cell_type in range(num_cell_types):
            true_methylation[cell_type, probe_idx] = base_value + cell_type * step
            
    # Random distinct variation probes
    start_idx = num_type_specific + num_gradual
    
    for probe_idx in range(start_idx, num_probes):
        # Generate random values that are still distinct between cell types
        values = np.random.random(num_cell_types) * 0.8 + 0.1  # Values between 0.1-0.9
        
        # Ensure at least 0.2 difference between any two cell types
        for i in range(num_cell_types):
            for j in range(i+1, num_cell_types):
                while abs(values[i] - values[j]) < 0.2:
                    # Adjust one of the values
                    if np.random.random() < 0.5:
                        values[i] = min(0.9, max(0.1, values[i] + np.random.choice([-0.2, 0.2])))
                    else:
                        values[j] = min(0.9, max(0.1, values[j] + np.random.choice([-0.2, 0.2])))
        
        # Assign the distinct values
        for cell_type in range(num_cell_types):
            true_methylation[cell_type, probe_idx] = values[cell_type]
    
    # Generate the observed methylation matrix
    observed_methylation = np.zeros((num_samples, num_probes))
    
    # For each sample, mix the cell type methylation patterns according to composition
    for i in range(num_samples):
        for j in range(num_probes):
            # Calculate methylation as weighted sum of cell type methylation values
            observed_methylation[i, j] = np.sum(
                true_compositions[i] * true_methylation[:, j]
            )
            
            # Add some noise if requested
            if add_noise:
                noise = np.random.normal(0, noise_level)
                observed_methylation[i, j] += noise
                observed_methylation[i, j] = np.clip(observed_methylation[i, j], 0.001, 0.999)
                observed_methylation[i, j] = round(observed_methylation[i, j], 3)  # Round to 3 decimal places
    
    # Create DataFrame for output
    methylation_df = pd.DataFrame(observed_methylation, columns=probe_ids)
    methylation_df.insert(0, "id", sample_ids)
    
    # Transpose the DataFrame to match CellCDecon input format (probes as rows)
    methylation_df = methylation_df.set_index("id").T.reset_index()
    methylation_df = methylation_df.rename(columns={"index": "id"})
    
    # Save methylation matrix to file
    methylation_df.to_csv(output_file, sep=" ", index=False)
    print(f"Methylation matrix saved to {output_file}")
    
    # Generate truth files if requested
    if generate_truth_files:
        # Save true cell compositions
        compositions_df = pd.DataFrame(true_compositions, 
                                     columns=[f"CellType_{i+1}" for i in range(num_cell_types)])
        compositions_df.insert(0, "Sample", sample_ids)
        compositions_file = output_file + ".cellcomp"
        compositions_df.to_csv(compositions_file, sep=" ", index=False)
        print(f"True cell compositions saved to {compositions_file}")
        
        # Save true methylation values
        means_df = pd.DataFrame(true_methylation.T, 
                              columns=[f"CellType_{i+1}" for i in range(num_cell_types)])
        means_df.insert(0, "Probe", probe_ids)
        means_file = output_file + ".means"
        means_df.to_csv(means_file, sep=" ", index=False)
        print(f"True methylation means saved to {means_file}")
    
    return output_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic methylation data for CellCDecon testing")
    parser.add_argument("--samples", type=int, default=50, help="Number of samples")
    parser.add_argument("--probes", type=int, default=200, help="Number of probes")
    parser.add_argument("--cell-types", type=int, default=3, help="Number of cell types")
    parser.add_argument("--output", type=str, default="synthetic_methylation.txt", help="Output file path")
    parser.add_argument("--noise", type=float, default=0.03, help="Noise level (standard deviation)")
    parser.add_argument("--concentration", type=float, default=0.5, 
                        help="Dirichlet concentration parameter (lower = more biased compositions)")
    parser.add_argument("--no-truth-files", action="store_false", dest="truth_files", 
                        help="Don't generate truth files")
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate the data
    output_file = generate_methylation_data(
        num_samples=args.samples,
        num_probes=args.probes,
        num_cell_types=args.cell_types,
        output_file=args.output,
        noise_level=args.noise,
        generate_truth_files=args.truth_files,
        concentration_param=args.concentration
    )
    
    print(f"\nTo use this file with CellCDecon, run:")
    print(f"./CellCDecon -k {args.cell_types} -n {args.samples} -f {args.output} -c 1") 
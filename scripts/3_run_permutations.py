import os
import itertools
import yaml
import subprocess
from multiprocessing import Pool

# Define the possible permutations
marker_types = ["DNase", "H3K4me3"]
rna_seq_files = {
    "adipocytes": "./data/rna_seq_processed/GSE249195_adipocytes.txt",
    "imac": "./data/rna_seq_processed/GSE220450_imac.txt",
    "k7": "./data/rna_seq_processed/GSE220450_k7.txt"
}

# Paths for essential files
essential_genes_file = "./data/essential_genes_ensemble.csv"

# Marker files dictionary to handle different marker file variants
marker_files = {
    "dnase": ["./results/2_basket_generation/dnase/dnase_15_cell_lines.csv"],
    "h3k4me3": [
        "./results/2_basket_generation/h3k4me3/h3k4me3_15_cell_lines.csv",
        "./results/2_basket_generation/h3k4me3/h3k4me3_random_5_cell_lines.csv",
        "./results/2_basket_generation/h3k4me3/h3k4me3_random_10_cell_lines.csv",
        "./results/2_basket_generation/h3k4me3/h3k4me3_random_12_cell_lines.csv"
    ]
}

# Directory to store generated config files
os.makedirs("config_files", exist_ok=True)

# Function to create directories recursively if they don't exist
def create_directory(path):
    """Creates directory recursively if it doesn't exist."""
    os.makedirs(path, exist_ok=True)

# Define the function to create a YAML config and run the Rscript
def run_rscript(config):
    # Adjust the config filename to include the variant structure
    config_filename = (
        f"./run_perm_configs/config_{config['variables']['marker_type']}_"
        f"{config['marker_variant']}_{config['rna_seq_cell_type']}.yaml"
    )

    # Save the configuration as a YAML file
    with open(config_filename, "w") as f:
        yaml.dump(config, f)

    # Ensure the output and results directories exist
    create_directory(config["paths"]["figure_output"])
    create_directory(config["paths"]["results_folder"])

    # Run the Rscript with the generated config
    command = f"Rscript scripts/3_cco_results.r {config_filename}"
    print(f"Running: {command}")
    subprocess.run(command, shell=True)

# Generate all permutations and create config dictionaries
configs = []
for marker_type, (cell_type, rna_seq_file) in itertools.product(marker_types, rna_seq_files.items()):
    for marker_file in marker_files[marker_type.lower()]:
        # Extract the marker variant (e.g., "random_5" or "full_15")
        if "random" in marker_file:
            marker_variant = "_".join(marker_file.split("_")[-3:-1])
        else:
            marker_variant = "full_15"

        # Define output paths based on the marker type, variant, and RNA-seq cell type
        figure_output = f"./results/3_cco_results/{marker_type.lower()}/{marker_variant}/{cell_type}"
        results_folder = figure_output  # Same path for results folder

        # Create the config dictionary for each permutation
        config = {
            "paths": {
                "essential_genes_file": essential_genes_file,
                "marker_merged_file": marker_file,
                "rnaseq_file": rna_seq_file,
                "figure_output": figure_output,
                "results_folder": results_folder,
            },
            "variables": {
                "marker_type": marker_type.lower(),
                "output_file_name": "gene_expression_plots"
            },
            "marker_variant": marker_variant,  # Store marker variant for filename/logging
            "rna_seq_cell_type": cell_type  # Store RNA-seq cell type for logging
        }
        configs.append(config)

# Run the configurations in parallel
if __name__ == "__main__":
    with Pool(processes=os.cpu_count()) as pool:  # Use all available CPUs
        pool.map(run_rscript, configs)

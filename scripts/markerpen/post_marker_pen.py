import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import yaml
import argparse

# Set up argument parser to accept config file as a command-line argument
parser = argparse.ArgumentParser(description='Run the Python part of the workflow.')
parser.add_argument('--config', required=True, help='Path to the configuration file (YAML format)')
args = parser.parse_args()

# Load the config file
with open(args.config, 'r') as file:
    config = yaml.safe_load(file)

# Function to parse the log file and extract the relevant data
def parse_cc_output(file_path):
    data = {
        "cco_level": [],
        "original": [],
        "out": [],
        "retention_percentage": []
    }

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for i in range(len(lines)):
            if "original =" in lines[i]:
                values = re.findall(r"\d+", lines[i])
                original = int(values[0])
                out = int(values[1])
                
                cco_level = int(re.search(r"CCO level (\d+)", lines[i + 1]).group(1))

                retention_percentage = (1 - (out / original)) * 100 if original != 0 else 100

                data["cco_level"].append(cco_level)
                data["original"].append(original)
                data["out"].append(out)
                data["retention_percentage"].append(retention_percentage)

    df = pd.DataFrame(data)
    return df.sort_values(by='cco_level')

# File paths from config
log_file_path = os.path.join(config['r_script']['output_folder'], config['python_script']['log_file'])
pdf_file_path = os.path.join(config['r_script']['output_folder'], config['python_script']['pdf_file'])

# Load and sort the data into a DataFrame
df = parse_cc_output(log_file_path)

with PdfPages(pdf_file_path) as pdf:
    plt.figure(figsize=(10, 5))
    plt.plot(df['cco_level'], df['original'], marker='o')
    plt.plot(df['cco_level'], df['out'], marker='o')
    plt.xlabel('CCO Level')
    plt.ylabel('Absolute Values')
    plt.title('Original vs Out (Absolute Values)')
    plt.grid(True)
    pdf.savefig()
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.plot(df['cco_level'], df['retention_percentage'], marker='o', color='purple')
    plt.xlabel('CCO Level')
    plt.ylabel('Retention Percentage')
    plt.title('Retention Percentage of Genes by MarkerPen')
    plt.grid(True)
    pdf.savefig()
    plt.close()

print(f"Plots saved to {pdf_file_path}")

# PART 3: Calculate Jaccard Index
folder_path = config['r_script']['output_folder']

def read_gene_set(file_path):
    df = pd.read_csv(file_path)
    return set(df['x'])

gene_sets = []
for file_name in os.listdir(folder_path):
    if file_name.endswith(".csv"):
        file_path = os.path.join(folder_path, file_name)
        gene_set = read_gene_set(file_path)
        gene_sets.append(gene_set)

intersection = set.intersection(*gene_sets)
union = set.union(*gene_sets)

jaccard_index = len(intersection) / len(union) if union else 0

print(f"Jaccard Index: {jaccard_index}")

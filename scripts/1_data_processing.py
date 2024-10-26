import os
import pandas as pd
import numpy as np
import pyarrow as pa
import matplotlib.pyplot as plt
import seaborn as sns
import pyarrow.csv as csv
import logging
import yaml
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import time
import subprocess
import argparse

class DataProcessor:
    def __init__(self, config):
        self.config = config
        self.slop_length = config['slop_length']
        self.markertype = config['markertype']
        self.tss_filepath = config['tss_filepath']
        self.humanchrom_filepath = config['humanchrom_filepath']
        self.bedfile_folder = config['bedfile_folder']
        self.indiv_bedfile = config['indiv_bedfile']
        self.percentiles = config['percentiles']
        self.output_paths = config['output_paths']
        self.log_file_path = config['log_filepath']
        self.setup_logging()

    def setup_logging(self):
        logging.basicConfig(filename=self.log_file_path, level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')

    def file_exists(self, filepath):
        return os.path.isfile(filepath)

    def filter_bed_file(self, file):
        filename = self.get_filename_without_extension(file)
        tmp_bed = os.path.join(self.output_paths['tmp_bed'], f"{filename}_{self.slop_length}.bed")
        if not self.file_exists(tmp_bed):
            logging.info(f"Filtering noise for {file}")
            self.filter_noise(file, tmp_bed)
        else:
            logging.info(f"Filtered BED file already exists: {tmp_bed}")
        return tmp_bed

    def filter_noise(self, file, output_file):
        header_names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']
        bed = pd.read_csv(file, sep='\t', header=None)
        bed.columns = header_names
        bed['length'] = bed['chromEnd'] - bed['chromStart']
        bed = bed[bed['length'] < np.percentile(bed['length'], self.percentiles['upper'])]
        bed = bed[bed['length'] > np.percentile(bed['length'], self.percentiles['lower'])]
        bed.drop('length', axis=1, inplace=True)
        bed.to_csv(output_file, sep='\t', header=None, index=None)

    def generate_coverage(self, tmp_bed):
        filename = self.get_filename_without_extension(tmp_bed)
        tmp_outname = os.path.join(self.output_paths['tmp_outname'], f"{filename}.txt")
        if not self.file_exists(tmp_outname):
            logging.info(f"Generating coverage for {tmp_bed}")
            self.create_coverage(tmp_bed, tmp_outname)
        else:
            logging.info(f"Coverage file already exists: {tmp_outname}")
        return tmp_outname

    def create_coverage(self, tmp_bed, output_file):
        cmd = (f"bedtools slop -s -l {self.slop_length} -r {self.slop_length} -i {self.tss_filepath} -g {self.humanchrom_filepath} | "
               f"bedtools coverage -d -a - -b {tmp_bed} > {output_file}")
        subprocess.run(cmd, shell=True, check=True)

    def convert_to_pileup(self, tmp_outname):
        filename = self.get_filename_without_extension(tmp_outname)
        pileup_res_outname = os.path.join(self.output_paths['pileup_res_outname'], f"{filename}.csv")
        if not self.file_exists(pileup_res_outname):
            logging.info(f"Converting coverage to pileup for {tmp_outname}")
            self.wrangle_coverage(tmp_outname, pileup_res_outname)
        else:
            logging.info(f"Pileup file already exists: {pileup_res_outname}")
        return pileup_res_outname

    def wrangle_coverage(self, infile, outfile):
        gene_id_index = 3
        gene_name_index = 13
        marker = self.markertype

        pileup_dict = {}
        with open(infile, "r") as ifh:
            for line in ifh:
                split_line = line.strip().split("\t")
                line_ele = [split_line[0], split_line[1], split_line[2], split_line[gene_id_index], split_line[gene_name_index], split_line[-1]]
                group_key = (line_ele[3], line_ele[4], marker)
                if group_key not in pileup_dict:
                    pileup_dict[group_key] = []
                try:
                    depth = int(line_ele[-1])
                except ValueError:
                    logging.error(f"Cannot convert depth at line into valid integer. Exiting.")
                    raise
                pileup_dict[group_key].append(depth)

        expected_length = self.slop_length * 2 + 1
        pileup_array_dict = {}
        for key, value in pileup_dict.items():
            if len(value) == expected_length:
                pileup_array_dict[key] = np.array(value).reshape(-1, 1)

        self.publish_depth_matrix(pileup_array_dict, expected_length, outfile)

    def publish_depth_matrix(self, pileup_array_dict, expected_length, outfile):
        gene_id_data = []
        gene_name_data = []
        marker_data = []
        arrays_to_stack = []

        for key, array in pileup_array_dict.items():
            gene_id, gene_name, marker = key
            gene_id_data.append(gene_id)
            gene_name_data.append(gene_name)
            marker_data.append(marker)
            arrays_to_stack.append(array)

        pileup_matrix = np.hstack(arrays_to_stack)
        base_data = [f"B{i+1}" for i in range(expected_length)]

        df = pd.concat([pd.DataFrame({"TRANSCRIPT": gene_id_data,
                                     "GENE": gene_name_data,
                                     "MARKER": marker_data}),
                        pd.DataFrame(data=pileup_matrix.T, columns=base_data)],
                       axis=1)

        df_pa_table = pa.Table.from_pandas(df)
        csv.write_csv(df_pa_table, outfile)

    def plot_pileup(self, pileup_res_outname):
        plot_res_outname = os.path.join(self.output_paths['plot_res_outname'], os.path.basename(pileup_res_outname).replace('.csv', ''))
        if not self.file_exists(f"{plot_res_outname}.png"):
            logging.info(f"Plotting pileup generated from {pileup_res_outname}")
            self.create_plot(pileup_res_outname, plot_res_outname)
        else:
            logging.info(f"Plot already exists: {plot_res_outname}.png")

    def create_plot(self, pileup_res_outname, plot_res_outname):
        bn = pd.read_csv(pileup_res_outname)
        sliced = bn.iloc[:, 3:]
        sliced.columns = [i - self.slop_length for i in range(1 + 2 * self.slop_length)]
        alpha = sliced.T.astype(bool).sum(axis=1)

        plt.figure(figsize=(16, 12))
        plt.rcParams.update({'font.size': 18})
        sns.lineplot(alpha)
        plt.title(f'Plot of {self.markertype} pileup within a window of {self.slop_length} base pairs up and downstream from TSS\n')
        plt.xlabel('Distance away from TSS')
        plt.ylabel('Number of markers present')
        plt.savefig(f"{plot_res_outname}.png")
        plt.close()

    def get_filename_without_extension(self, filepath):
        return Path(filepath).stem

    def process_file(self, file):
        start_time = time.time()
        logging.info(f"Processing file: {file}")

        tmp_bed = self.filter_bed_file(file)
        tmp_outname = self.generate_coverage(tmp_bed)
        pileup_res_outname = self.convert_to_pileup(tmp_outname)
        self.plot_pileup(pileup_res_outname)

        end_time = time.time()
        logging.info(f"Completed processing of {file}. Time taken: {end_time - start_time:.2f} seconds.")

    def run(self):
        files = []
        if self.bedfile_folder:
            files += [os.path.join(self.bedfile_folder, f) for f in os.listdir(self.bedfile_folder) if f.endswith('.bed')]
        if self.indiv_bedfile:
            files += [self.indiv_bedfile]
        with ProcessPoolExecutor() as executor:
            executor.map(self.process_file, files)

def main():
    print("Script 1_data_processing.py has started")
    parser = argparse.ArgumentParser(description="Run the Basket Analysis")
    parser.add_argument('-c', "--config", required=True, help="Path to the YAML configuration file")
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    if not Path(config['output_paths']['pileup_res_outname']).exists():
        Path(config['output_paths']['pileup_res_outname']).mkdir(parents=True, exist_ok=True)
    if not Path(config['output_paths']['plot_res_outname']).exists():
        Path(config['output_paths']['pileup_res_outname']).mkdir(parents=True, exist_ok=True)
    processor = DataProcessor(config)
    processor.run()

if __name__ == "__main__":
    main()

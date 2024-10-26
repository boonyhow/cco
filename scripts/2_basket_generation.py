import os
import pandas as pd
import numpy as np
import pathlib
import argparse
import yaml
import logging
from tqdm import tqdm
import concurrent.futures
import random

class BasketAnalysis:
    def __init__(self, config_path, sample_size=None, sample_method='random', indices=None):
        with open(config_path, 'r') as file:
            self.config = yaml.safe_load(file)
        
        # Setup logging
        logging.basicConfig(filename=self.config['log_file'], 
                            level=logging.INFO, 
                            format='%(asctime)s - %(levelname)s - %(message)s')

        self.sample_size = sample_size  # How many files to sample
        self.sample_method = sample_method  # 'random' or 'index'
        self.indices = indices  # List of indices if sampling by index

    def generate_vectors(self, row, tss_pos, boundary, proximal_boundary=None):
        len_of_pileups = len(row)
        positions = {}

        if proximal_boundary:
            end_of_enhancer_left = tss_pos - boundary - 1
            end_of_proximal_promoters_left = tss_pos - proximal_boundary 
            start_of_proximal_promoters_right = tss_pos + proximal_boundary 
            start_of_enhancer_right = tss_pos + boundary
            positions = {
                0: (0, end_of_enhancer_left),
                1: (end_of_enhancer_left, end_of_proximal_promoters_left),
                2: (end_of_proximal_promoters_left, tss_pos),
                3: (tss_pos, start_of_proximal_promoters_right),
                4: (start_of_proximal_promoters_right, start_of_enhancer_right),
                5: (start_of_enhancer_right, len_of_pileups)
            }
        else:
            end_of_first = tss_pos - boundary - 2
            end_of_second = tss_pos + boundary + 1
            positions = {
                0: (0, end_of_first),
                1: (end_of_first, tss_pos),
                2: (tss_pos, end_of_second),
                3: (end_of_second, len_of_pileups)
            }

        out = np.zeros(len(positions))
        non_zero_indices = np.nonzero(row.values)[0]

        for i in non_zero_indices:
            for j in positions:
                if positions[j][0] <= i < positions[j][1]:
                    out[j] += 1
                    break
                if i == len_of_pileups - 1:
                    out[-1] += 1
                    break

        return out.tolist()

    def split_into_regions(self, file, tss_path, out_folder, boundary, proximal_boundary=None):
        tss = pd.read_csv(tss_path, sep='\t', header=None)
        tss.rename(columns={3: 'gene_id', 5:'strand'}, inplace=True)
        tss.index = tss['gene_id'].apply(lambda x: x.split('.')[0])
        tss = tss.iloc[:, 5]
        tss = tss[~tss.index.duplicated(keep='first')]

        pileup = pd.read_csv(file)
        pileup.index = pileup['TRANSCRIPT'].apply(lambda x: x.split('.')[0])
        pileup = pileup.iloc[:, 3:]

        results = pd.DataFrame()
        results['raw'] = pileup.apply(self.generate_vectors, axis=1, 
                                      tss_pos=1 + ((len(pileup.columns) - 1) // 2),
                                      boundary=boundary,
                                      proximal_boundary=proximal_boundary)
        results = results.merge(tss, left_index=True, right_index=True)

        if proximal_boundary:
            results[['enhancer_left', 'proximal_promoter_left', 'core_promoter_left', 
                     'core_promoter_right', 'proximal_promoter_right','enhancer_right']] = results['raw'].to_list()
        else:
            results[['enhancer_left', 'promoter_left', 'promoter_right', 'enhancer_right']] = results['raw'].to_list()
            
        if not pathlib.Path(out_folder).exists():
            pathlib.Path(out_folder).mkdir(parents=True, exist_ok=True)
        
        output_path = f"{out_folder}/{os.path.basename(file).replace('.csv', '_split_data.csv')}"
        results.to_csv(output_path)

        logging.info(f"File {file} split into regions and saved to {output_path}")
        return output_path

    def generate_baskets(self, split_data_fp, region, out_folder, fixedL, downstream):
        results = pd.read_csv(split_data_fp, index_col=0)
        results.index.name = 'gene_id'

        if fixedL:
            if downstream == True:
                results['in_basket'] = results.apply(lambda x: 
                                                     (x['strand'] == '+' and x[f'{region}_right'] > fixedL) or 
                                                     (x['strand'] == '-' and x[f'{region}_left'] > fixedL), axis=1)
            else:
                results['in_basket'] = results.apply(lambda x: 
                                                     (x['strand'] == '-' and x[f'{region}_right'] > fixedL) or 
                                                     (x['strand'] == '+' and x[f'{region}_left'] > fixedL), axis=1)

        if not pathlib.Path(out_folder).exists():
            pathlib.Path(out_folder).mkdir(parents=True, exist_ok=True)

        direction = 'downstream' if downstream == True else 'upstream'
        output_path = f"{out_folder}/{region}_fixedL_{fixedL}_{direction}_{os.path.basename(split_data_fp).replace('_split_data.csv', '_basket.csv')}"
        results.to_csv(output_path)

        logging.info(f"Baskets generated for {split_data_fp} and saved to {output_path}")

    # def combine_baskets(self, dir_path, out_folder, output_name):
    #     gencode = pd.read_csv(self.config['gencode_path'], sep='\t', header=None)
    #     df = pd.DataFrame()
    #     df['gene_id'] = gencode.iloc[:, 3].apply(lambda x: x.split('.')[0]).drop_duplicates()

    #     for file in os.listdir(dir_path):
    #         try:
    #             tmp_df = pd.read_csv(os.path.join(dir_path, file))
    #             subdf = tmp_df.loc[:, ['gene_id', 'in_basket']].drop_duplicates(subset='gene_id')
    #             subdf.rename(columns={'in_basket': f"{'_'.join(file.split('_')[4:])}"}, inplace=True)
    #             df = subdf.merge(df, on='gene_id', how='left')
    #         except IsADirectoryError:
    #             pass

    #     df = df.applymap(lambda x: str(x).upper())
    #     output_path = out_folder
    #     output_path = pathlib.Path(f"{out_folder}/{output_name}.csv")
    #     df.to_csv(output_path)
    #     logging.info(f"Combined baskets saved to {output_path}")
    
    def combine_baskets(self, basket_files, out_folder, output_name):
        gencode = pd.read_csv(self.config['gencode_path'], sep='\t', header=None)
        df = pd.DataFrame()
        df['gene_id'] = gencode.iloc[:, 3].apply(lambda x: x.split('.')[0]).drop_duplicates()

        for file in basket_files:
            try:
                tmp_df = pd.read_csv(os.path.join(self.config['basket_out'], file))
                subdf = tmp_df.loc[:, ['gene_id', 'in_basket']].drop_duplicates(subset='gene_id')
                subdf.rename(columns={'in_basket': f"{'_'.join(file.split('_')[4:])}"}, inplace=True)
                df = subdf.merge(df, on='gene_id', how='left')
            except IsADirectoryError:
                pass

        df = df.applymap(lambda x: str(x).upper())
        output_path = pathlib.Path(f"{out_folder}/{output_name}.csv")
        df.to_csv(output_path)
        logging.info(f"Combined baskets saved to {output_path}")

        
    def file_exists(self, filepath):
        return os.path.isfile(filepath)
    
    def process_file(self, file):
        file_path = os.path.join(self.config['pileup_folder'], file)
        if not self.file_exists(f"{self.config['split_out']}/{os.path.basename(file).replace('.csv', '_split_data.csv')}"):
            split_data_fp = self.split_into_regions(file_path, 
                                                    self.config['tss_path'], 
                                                    self.config['split_out'], 
                                                    self.config['enhancer_promoter_boundary'], 
                                                    self.config['proximal_core_boundary'] if self.config['proximal'] else None)
        else:
            logging.warning(f'Earlier split of {file} has been found and the code will be using that version.')
            split_data_fp = f"{self.config['split_out']}/{os.path.basename(file).replace('.csv', '_split_data.csv')}"
            
        direction = 'downstream' if self.config['downstream'] == True else 'upstream'
        if not self.file_exists(f"{self.config['basket_out']}/{self.config['region_for_basket']}_fixedL_{self.config['length_for_basket']}_{direction}_{os.path.basename(split_data_fp).replace('_split_data.csv', '_basket.csv')}"):
            self.generate_baskets(split_data_fp, 
                                self.config['region_for_basket'], 
                                self.config['basket_out'], 
                                self.config['length_for_basket'], 
                                self.config['downstream'])
        else:
            logging.warning(f'Earlier baskets of {file} has been found, no new baskets will be generated')
            
    def run(self):
        # Get all files in the directory
        files_to_process = [f for f in os.listdir(self.config['pileup_folder']) if f.endswith('.csv')]
        original_length = len(files_to_process)
        output_name = self.config['marker_type']
        
        # Track basket files to pass to combine_baskets
        basket_files = []
        
        # Apply sampling based on method
        if self.sample_method == 'random':
            if self.sample_size and self.sample_size <= original_length:
                files_to_process = random.sample(files_to_process, self.sample_size)
                output_name += f"_random_{self.sample_size}_cell_lines"
            else:
                logging.warning(f"Sample size {self.sample_size} is larger than available files. Using all files.")
                output_name += f"_{original_length}_cell_lines"
        elif self.sample_method == 'index':
            if self.indices:
                files_to_process = [files_to_process[i] for i in self.indices if i < original_length]
                index_string = '_'.join(map(str, self.indices))
                output_name += f"_index_{index_string}_length_{len(self.indices)}_cell_lines"
            else:
                logging.warning(f"No valid indices provided, using all files.")
                output_name += f"_{original_length}_cell_lines"
        else:
            output_name += f"_{original_length}_cell_lines"
        
        logging.info(f"Number of files to be analyzed: {len(files_to_process)}")

        # Process the files using parallel execution
        with concurrent.futures.ProcessPoolExecutor() as executor:
            list(tqdm(executor.map(self.process_file, files_to_process), total=len(files_to_process), unit="file"))

        # Collect basket files that were generated during processing
        for file in files_to_process:
            # Dynamically generate the basket file name using the YAML config values
            direction = 'downstream' if self.config['downstream'] else 'upstream'
            basket_file = f"{self.config['region_for_basket']}_fixedL_{self.config['length_for_basket']}_{direction}_{os.path.basename(file).replace('.csv', '_basket.csv')}"
            basket_files.append(basket_file)

        # Combine the results with the appropriate file name
        self.combine_baskets(basket_files, self.config['merged_out'], output_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Basket Analysis")
    parser.add_argument('-c',"--config", required=True, help="Path to the YAML configuration file")
    parser.add_argument("--sample-size", type=int, help="Number of files to sample from the input folder (default is all files)", default=None)
    parser.add_argument("--sample-method", type=str, choices=['random', 'index'], default='random', help="Method to select files: 'random' or 'index'")
    parser.add_argument("--indices", type=int, nargs='+', help="Indices of files to select if using 'index' method", default=None)
    args = parser.parse_args()

    analysis = BasketAnalysis(args.config, args.sample_size, args.sample_method, args.indices)
    analysis.run()

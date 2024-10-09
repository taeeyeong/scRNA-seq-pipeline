import os
import sys
import argparse
import logging

class ScRNAseqPipeline:
    def __init__(self, data_path, output_dir, n_dims=6, random_state=42):
        self.data_path = data_path
        self.output_dir = output_dir
        self.n_dims = n_dims
        self.random_state = random_state
        self.adata = None
        
        # logging
        self.logger = logging.getLogger('ScRNAseqPipeline')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # stream handler
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(stream_handler)
        
        # file handler
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        file_handler = logging.FileHandler(os.path.join(output_dir, 'sc_rna_pipeline.log'))
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        
        self.logger.info('Pipeline initialized')
    
    def load_data(self):
        """_summary_
            load data from data_path
            support files format: mtx, h5ad
        """
        try:
            self.logger.info('Loading data from {}'.format(self.data_path))
            
            # if file type is not offered, detect file type automatically
            if self.file_type is None:
                if os.path.isdir(self.data_path):
                    self.file_type = '10x_mtx'
                elif self.data_path.endswith('.h5ad'):
                    self.file_type = 'h5ad'
                else:
                    self.logger.error('Unsupported data format')
                    sys.exit(1)
        
        except Exception as e:
            self.logger.error('Error loading data: {}'.format(e))
            sys.exit(1)
        
def parse_argments():
    parser = argparse.ArgumentParser(description="Single cell RNA-seq pipeline")
    parser.add_argument('--data_path', type=str, required=True, help='Path to Cell Ranger output directory (filtered_feature_bc_matrix) or h5ad data format')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save outputs')
    parser.add_argument('--n_dims', type=int, default=6, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    return parser.parse_args()
    
def main():
    args = parse_argments()
    
    pipeline = ScRNAseqPipeline()

if __name__ == "__main__":
    main()
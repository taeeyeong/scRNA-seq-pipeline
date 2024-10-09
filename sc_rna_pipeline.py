import os
import sys
import argparse
import logging
import scanpy as sc

from multiprocessing import Pool, cpu_count
from scipy import sparse
from scipy.stats import pearsonr

class ScRNAseqPipeline:
    def __init__(self, data_path, output_dir, gene, n_dims=6, random_state=42):
        self.data_path = data_path
        self.output_dir = output_dir
        self.gene = gene
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
            
            if os.path.isdir(self.data_path):
                self.file_type = '10x_mtx'
            elif self.data_path.endswith('.h5ad'):
                self.file_type = 'h5ad'
            else:
                self.logger.error('Unsupported data format')
                sys.exit(1)
            
            if self.file_type == '10x_mtx':
                self.adata = sc.read_10x_mtx(
                    self.data_path,
                    var_names='gene_symbols',
                    cache=True
                )
            elif self.file_type == 'h5ad':
                self.adata = sc.read(self.data_path)
            else:
                self.logger.error('Unsupported data format')
                sys.exit(1)
        
        except Exception as e:
            self.logger.error('Error loading data: {}'.format(e))
            sys.exit(1)
    
    def quality_control(self):
        """_summary_
            perform quality control on the data
        """
        self.logger.info('Starting quality control')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during quality control: {}'.format(e))
            sys.exit(1)
    
    def preprocess_data(self):
        """_summary_
            Normalization and scaling of the data
        """
        self.logger.info('Starting data preprocessing')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during data preprocessing: {}'.format(e))
            sys.exit(1)
    
    def run_pca(self):
        """_summary_
            Run PCA on the data and visualize the results
        """
        self.logger.info('Starting PCA')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during PCA: {}'.format(e))
            sys.exit(1)
    
    def clustering(self):
        """_summary_
            Perform finding neighbors and clustering on the data
        """
        self.logger.info('Starting clustering')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during clustering: {}'.format(e))
            sys.exit(1)
    
    def visualize_clusters(self):
        """_summary_
            Visualize the clustering results using UMAP
        """
        self.logger.info('Starting visualization of clusters')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during visualization of clusters: {}'.format(e))
            sys.exit(1)
    
    def annotate_cell_types(self):
        """_summary_
            Annotate cell types based on Celltypist
        """
        self.logger.info('Starting cell type annotation')
        try:
            pass
        except Exception as e:
            self.logger.error('Error during cell type annotation: {}'.format(e))
            sys.exit(1)
    
    def compute_correlation(self, gene):
        """_summary_
            Compute correlation between gene of interest and the other genes
            
           _note_
           	•벡터 연산 활용: 루프를 제거하고 Numpy 배열 연산으로 대체하여 계산 속도를 높였습니다.
	        •메모리 관리: toarray()를 사용하여 희소 행렬을 밀집 행렬로 변환하므로, 메모리 사용량이 증가할 수 있습니다. 
                데이터가 매우 큰 경우 아래의 배치 처리 또는 희소 행렬 연산 방법을 고려하세요.
            - multi processing: 여러 프로세스에서 병렬로 상관계수 계산 
        """
        gene_exp = self.adata[:, self.gene].X
        gene_data = self.adata[:, gene].X
        if sparse.issparse(gene_exp):
            gene_exp = gene_exp.toarray().flatten()
        else:
            gene_exp = gene_exp.flatten()
        if sparse.issparse(gene_data):
            gene_data = gene_data.toarray().flatten()
        else:
            gene_data = gene_data.flatten()
        corr, p_value = pearsonr(gene_exp, gene_data)
        return (gene, corr, p_value)
    
    def save_results(self):
        """_summary_
            Save the results of the pipeline
        """
        self.logger.info('Saving results')
        try:
            pass
        except Exception as e:
            self.logger.error('Error saving results: {}'.format(e))
            sys.exit(1)
    
    def run_pipeline(self):
        """_summary_
            Run the entire pipeline
        """
        self.load_data()
        self.quality_control()
        self.preprocess_data()
        # self.run_pca()
        # self.clustering()
        # self.visualize_clusters()
        self.annotate_cell_types()
        self.save_results()
        self.logger.info('Pipeline execution completed successfully')
        
def parse_argments():
    parser = argparse.ArgumentParser(description="Single cell RNA-seq pipeline")
    parser.add_argument('--data_path', type=str, required=True, help='Path to Cell Ranger output directory (filtered_feature_bc_matrix) or h5ad data format')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save outputs')
    parser.add_argument('--gene', type=str, help='Gene of interest to compute correlation')
    parser.add_argument('--n_dims', type=int, default=6, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    return parser.parse_args()
    
def main():
    args = parse_argments()
    
    pipeline = ScRNAseqPipeline(
        data_path=args.data_path,
        output_dir=args.output_dir,
        gene=args.gene,
        n_dims=args.n_dims,
        random_state=args.random_state
    )
    pipeline.run_pipeline()
    
    all_genes = pipeline.adata.var_names.tolist()
    
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(pipeline.compute_correlation, all_genes)
    
    correlations = results
    significant_genes = [item for item in correlations if item[2] < 0.05]
    significant_genes.sort(key=lambda x: abs(x[1]), reverse=True)
    
    for gene, corr, p_value in significant_genes[:30]:
        print(f'{gene}: corr={corr}, p-value={p_value}')
    
if __name__ == "__main__":
    main()
import os
import sys
import argparse
import logging
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import celltypist
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
    
    # TODO: 시각화 skip 기능 추가 
    def quality_control(self):
        """_summary_
            perform quality control on the data
        """
        self.logger.info('Starting quality control')
        try:
            # Mitochondrial gene ratio
            self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(
                self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
            )
            # visualize quality control metrics
            sc.pl.violin(
                self.adata, 
                ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                jitter=0.4, 
                multi_panel=True, 
                save='_qc_violin.png'
            )
            # visualize n_genes_by_counts distribution
            sns.histplot(self.adata.obs['n_genes_by_counts'], bins=50)
            plt.xlabel('n_genes_by_counts')
            plt.ylabel('Number of cells')
            plt.title('Distribution of n_genes_by_counts')
            plt.savefig(os.path.join(self.output_dir, 'n_genes_by_counts_distribution.png'))
            plt.close()
            
            # lowest threshold for number of genes and counts
            X = np.percentile(self.adata.obs['n_genes_by_counts'], 5) # lower threshold: 5%
            Y = np.percentile(self.adata.obs['total_counts'], 99) # upper threshold: 1%
            
            self.logger.info(f"Calculated lower threshold (X): {X}")
            self.logger.info(f"Calculated upper threshold (Y): {Y}")
            
            # visualize pct_counts_mt distribution
            sns.histplot(self.adata.obs['pct_counts_mt'], bins=50)
            plt.xlabel('pct_counts_mt')
            plt.ylabel('Number of cells')
            plt.title('Distribution of pct_counts_mt')
            plt.savefig(os.path.join(self.output_dir, 'pct_counts_mt_distribution.png'))
            plt.close()
            
            Z = np.percentile(self.adata.obs['pct_counts_mt'], 90) # upper threshold: 90%
            self.logger.info(f"Calculated upper threshold for pct_counts_mt (Z): {Z}")
            
            # Filtering
            initial_cell_count = self.adata.n_obs
            self.adata = self.adata[self.adata.obs.n_genes_by_counts > X, :]
            self.adata = self.adata[self.adata.obs.total_counts < Y, :]
            self.adata = self.adata[self.adata.obs.pct_counts_mt < Z, :]
            filtered_cell_count = self.adata.n_obs
            
            self.logger.info(f"Cells before QC: {initial_cell_count}")
            self.logger.info(f"Cells after QC: {filtered_cell_count}")
            self.logger.info('Quliaity control completed')
        except Exception as e:
            self.logger.error('Error during quality control: {}'.format(e))
            sys.exit(1)
    
    def preprocess_data(self):
        """_summary_
            Normalization and scaling of the data
        """
        self.logger.info('Starting data preprocessing')
        try:
            # Normalization
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            # Logarithmize the data
            sc.pp.log1p(self.adata)
            
            # TODO: Scaling method 추가 (Celltypist 할땐 skip 해야함)
            self.logger.info('Data preprocessing completed')
            # Scaling
            
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
            celltypist.models.download_models()
            # TODO: 모델 선택 기능 추가
            model = celltypist.models.Model.load(model='Immune_All_Low.pkl')
            predictions = celltypist.annotate(
                self.adata,
                model=model,
                majority_voting=True,
            )
            self.adata.obs['cell_type'] = predictions.predicted_labels.predicted_labels
            # visualization
            sc.pl.umap(self.adata, color='cell_type', save='_celltype_annotation.png')
            self.logger.info('Cell type annotation completed')
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

            Args:
                gene: gene compute correlation with gene of interest
    
            Note:
           	    - 벡터 연산 활용: 루프를 제거하고 Numpy 배열 연산으로 대체하여 계산 속도를 높임
	            - 메모리 관리: toarray()를 사용하여 희소 행렬을 밀집 행렬로 변환하므로, 메모리 사용량이 증가할 수 있음
                    데이터가 매우 큰 경우 아래의 배치 처리 또는 희소 행렬 연산 방법을 고려
                - multi processing: 여러 프로세스에서 병렬로 상관계수 계산 
        """
        try: 
            # TODO: multiprocessing 버전으로 수정 
            self.logger.info(f'Computing correlation between {self.gene} and {gene}')
            gene_exp = self.adata[:, self.gene].X
            if sparse.issparse(gene_exp):
                gene_exp = gene_exp.toarray().flatten()
            else:
                gene_exp = gene_exp.flatten()
            gene_data = self.adata[:, gene].X
            if sparse.issparse(gene_data):
                gene_data = gene_data.toarray().flatten()
            else:
                gene_data = gene_data.flatten()
            corr, p_value = pearsonr(gene_exp, gene_data)
            return (gene, corr, p_value)
        except Exception as e:
            self.logger.error(f'Error computing correlation between {self.gene} and {gene}: {e}')
            return (gene, np.nan, np.nan)
    
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
        self.run_pca()
        self.clustering()
        self.visualize_clusters()
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
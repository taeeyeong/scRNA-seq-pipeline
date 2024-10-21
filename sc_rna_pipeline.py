# TODO: sphinx setting
import os
import sys
import argparse
import logging
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import celltypist
from multiprocessing import current_process
from scipy.stats import pearsonr
from logging.handlers import QueueHandler

class ScRNAseqPipeline:
    def __init__(
        self, data_path, output_dir, gene, n_dims=4, random_state=42, skip_qc=False,
        only_highly_variable_genes=False, skip_preprocessing=False, skip_pca=False, 
        skip_clustering=False, skip_visualization=False, skip_annotation=False, 
        skip_correlation=False
        ):
        self.data_path = data_path
        self.output_dir = output_dir
        self.gene = gene
        self.n_dims = n_dims
        self.random_state = random_state
        self.adata = None
        
        self.only_highly_variable_genes = only_highly_variable_genes
        self.skip_qc = skip_qc
        self.skip_preprocessing = skip_preprocessing
        self.skip_pca = skip_pca
        self.skip_clustering = skip_clustering
        self.skip_visualizaton = skip_visualization
        self.skip_annotation = skip_annotation
        self.skip_correlation = skip_correlation
        
        self.logger = logging.getLogger('ScRNAseqPipeline')
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(stream_handler)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        file_handler = logging.FileHandler(
            os.path.join(output_dir, 'sc_rna_pipeline.log')
        )
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
            
            if self.only_highly_variable_genes:
                tg_mask = self.adata.var_names.isin([self.gene])
                hvg_mask = self.adata.var.highly_variable
                combined_mask = tg_mask | hvg_mask 
                self.adata = self.adata[:, combined_mask].copy()
        #TODO: loaded data 출력 (cells, genes)
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
            
            # find highly variable genes
            sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
            
            # TODO: Scaling method 추가 (Celltypist 할땐 skip 해야함)
            self.logger.info('Data preprocessing completed')
            # Scaling
            
        except Exception as e:
            self.logger.error('Error during data preprocessing: {}'.format(e))
            sys.exit(1)

    def annotate_cell_types(self):
        """_summary_
            Annotate cell types based on Celltypist
        """
        self.logger.info('Starting cell type annotation')
        try:
            celltypist.models.download_models()
            # TODO: 모델 선택 기능 추가
            model = celltypist.models.Model.load(model='Immune_All_High.pkl')
            predictions = celltypist.annotate(
                self.adata,
                model=model,
                majority_voting=True,
            )
            self.adata.obs['cell_type'] = predictions.predicted_labels.predicted_labels
        except Exception as e:
            self.logger.error('Error during cell type annotation: {}'.format(e))
            sys.exit(1)
    
    def run_pca(self):
        """_summary_
            Run PCA on the data and visualize the results
        """
        self.logger.info('Starting PCA')
        try:
            sc.pp.scale(self.adata, max_value=10)
            sc.pp.pca(self.adata, svd_solver='arpack')
            sc.pl.pca_variance_ratio(self.adata, log=False, save='_pc_elbow_plot.png')
        except Exception as e:
            self.logger.error('Error during PCA: {}'.format(e))
            sys.exit(1)
    
    def add_gene_expression_to_obs(self, genes):
        """_summary_

        Args:
            genes (List[str]): List of gene names

        Returns:
            _type_: _description_
        """
        if isinstance(genes, str):
            genes = [genes]
        
        status_columns = []
        for gene in genes:
            if gene not in self.adata.var_names:
                self.logger.error(f'Gene {gene} not found in the dataset.')
                sys.exit(1)
            gene_exp = self.adata[:, gene].X.toarray().flatten()
            gene_status = np.where(gene_exp > 0, 'Expressed', 'Not expressed')
            gene_label = [f'{gene}: {status}' for status in gene_status]
            column_name = f'{gene}_status'
            self.adata.obs[column_name] = gene_label
            status_columns.append(column_name)
            self.logger.info(f'Added {column_name} to adata.obs')
        combined_status = self.adata.obs[status_columns].agg(' | '.join, axis=1)
        combined_label_name = '_'.join(genes) + '_combined_status'
        self.adata.obs[combined_label_name] = combined_status
        self.logger.info(f'Added {combined_label_name} to adata.obs')
        
    # TODO: tSNE pca 과정 필요
    def run_tsne(self, color_by, color_map=None, title=None):
        """_summary_
            Run tSNE and plot the results with specified labels and colors.

            Args:
                color_by (str): the column name in adata.obs. the standard to color in the plot
                color_map (dict or str, optiional): color pallette, or dictionary of color mapping. 
                    Default: None
        """
        self.logger.info('Starting tSNE colored by {color_by}')
        if title is None:
            title = f't-SNE plot Colored by {color_by}'
        try:
            sc.tl.tsne(self.adata, n_pcs=self.n_dims)
            cmap=plt.colormaps['viridis']
            sc.pl.tsne(
                self.adata,
                color=color_by,
                cmap=cmap,
                size=5,
                title=title,
                legend_loc='right margin',
                save='_{color_by}_tSNE.png' 
            )
        except Exception as e:
            self.logger.error('Error during tSNE: {}'.format(e))
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
            # visualization
            # sc.pl.umap(self.adata, color='cell_type', save='_celltype_annotation.png')
            # self.logger.info('Cell type annotation completed')
        except Exception as e:
            self.logger.error('Error during visualization of clusters: {}'.format(e))
            sys.exit(1)
    
    def compute_correlation(self, gene, top_n=30):
        """_summary_
            Compute correlation between gene of interest and the other genes

            Args:
                gene: gene compute correlation with gene of interest
        """
        try: 
            self.logger.info(f'Process {current_process().name}: Computing correlation between {self.gene} and {gene}')
            if self.only_highly_variable_genes:
                tg_mask = self.adata.var_names.isin(self.gene)
                hvg_mask = self.adata.var.highly_variable
                combined_mask = tg_mask | hvg_mask 
                self.adata = self.adata[:, combined_mask].copy()
            gene_exp = self.adata[:, self.gene].X.toarray().flatten()
            all_genes = self.adata.var_names
            correlations = []
            for gene in all_genes:
                gene_data = self.adata[:, gene].X.toarray().flatten()
                corr, p_value = pearsonr(gene_exp, gene_data)
                correlations.append(gene, corr, p_value)
                print(f'{gene}: corr={corr}, p-value={p_value}')
            significant_genes = [item for item in correlations if item[2] < 0.05]
            top_genes = sorted(significant_genes, key=lambda x: -abs(x[1]))[:top_n]
            return top_genes
        
        except Exception as e:
            self.logger.error(f'Error computing correlation between {self.gene} and {gene}: {e}')
            return (gene, np.nan, np.nan)
    
    def visualize_correlated_genes(self):
        try:
            self.logger.info(f'Process {current_process().name}: Visualize the top 30 highly correlated genes with {self.gene}')
            
        except Exception as e:
            self.logger.error(f'Error visualizing correlated genes: {e}')           
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
    
    # TODO: 일부 단계는 이전단계가 필요. 의존성 고려해서 사용자에게 경고 메시지 출력 | 강제로 이전 단계 실행
    def run_pipeline(self):
        """_summary_
            Run the entire pipeline
        """
        self.load_data()
        if not self.skip_qc:
            self.quality_control()
        else:
            self.logger.info('Skipping quality control step')
            
        if not self.skip_preprocessing:
            self.preprocess_data()
        else:
            self.logger.info('Skipping preprocessing step')
        
        if not self.skip_annotation:
            self.annotate_cell_types()
        else:
            self.logger.info('Skipping cell type annotation step')
                  
        if not self.skip_pca:
            self.run_pca()
        else:
            self.logger.info('Skipping PCA step')
        
        if not self.skip_clustering:
            self.clustering()
        else:
            self.logger.info('Skipping clustering step')
        self.add_gene_expression_to_obs(['AR', 'METTL3']) # TODO: 수정
        #TODO: add skip step 
        # self.run_tsne(color_by='AR_status') # TODO: 수정
        self.run_tsne(color_by='AR_METTL3_combined_status') # TODO: 수정
        
        
        if not self.skip_visualizaton:
            self.visualize_clusters()
        else:
            self.logger.info('Skipping visualization step')
        
        if not self.skip_correlation:
            self.compute_correlation()
        else:
            self.logger.info('Skipping correlation computation step')
        self.save_results()
        self.logger.info('Pipeline execution completed successfully')

def init_process(log_queue):
    logger = logging.getLogger('ScRNAseqPipeline')
    logger.handlers = []
    logger.setLevel(logging.INFO)
    queue_handler = QueueHandler(log_queue)
    logger.addHandler(queue_handler)  
    
def parse_arguments():
    parser = argparse.ArgumentParser(description="Single cell RNA-seq pipeline")
    parser.add_argument('--data_path', type=str, required=True, help='Path to Cell Ranger output directory (filtered_feature_bc_matrix) or h5ad data format')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save outputs')
    parser.add_argument('--gene', type=str, help='Gene of interest to compute correlation')
    parser.add_argument('--n_dims', type=int, default=4, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    
    parser.add_argument('--only_highly_variable_genes', action='store_true', help='Using only highly variable genes in all analysis steps')
    parser.add_argument('--skip_qc', action='store_true', help='Skip the quality control step')
    parser.add_argument('--skip_preprocessing', action='store_true', help='Skip the preprocessing step')
    parser.add_argument('--skip_pca', action='store_true', help='Skip the PCA step')
    parser.add_argument('--skip_clustering', action='store_true', help='Skip the clustering step')
    parser.add_argument('--skip_visualization', action='store_true', help='Skip the visualization step')
    parser.add_argument('--skip_annotation', action='store_true', help='Skip the cell type annotation step')
    parser.add_argument('--skip_correlation', action='store_true', help='Skip the correlation computation step')
    
    return parser.parse_args()
    
def main():
    args = parse_arguments()
       
    pipeline = ScRNAseqPipeline(
        data_path=args.data_path,
        output_dir=args.output_dir,
        gene=args.gene,
        n_dims=args.n_dims,
        random_state=args.random_state,
        only_highly_variable_genes=args.only_highly_variable_genes,
        skip_qc=args.skip_qc,
        skip_preprocessing=args.skip_preprocessing,
        skip_pca=args.skip_pca,
        skip_clustering=args.skip_clustering,
        skip_visualization=args.skip_visualization,
        skip_annotation=args.skip_annotation,
        skip_correlation=args.skip_correlation
    )
    pipeline.run_pipeline()
    
if __name__ == "__main__":
    main()
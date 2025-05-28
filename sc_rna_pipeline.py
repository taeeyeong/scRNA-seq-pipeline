# TODO: sphinx setting
import os
import sys
import argparse
import logging
import numpy as np
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from typing import Dict, List

import celltypist


class ScRNAseqPipeline:
    def __init__(
        self, data_paths, output_dir, genes_file=None, meta_file=None,target_gene=None, reference_gene=None, 
        n_dims=4, random_state=42, only_highly_variable_genes=False, steps=[]
        ):
        self.data_paths = data_paths
        self.genes_file = genes_file
        self.meta_file = meta_file
        self.output_dir = output_dir
        self.target_gene = target_gene
        self.reference_gene = reference_gene
        self.n_dims = n_dims
        self.random_state = random_state
        self.adata = None
        
        self.only_highly_variable_genes = only_highly_variable_genes
        self.steps = steps
        
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
        
        sc.settings.figdir = self.output_dir
        
    def load_data(self):
        """_summary_
            load data from data_paths
            support files format: mtx, h5ad
        """
        try:
            self.logger.info('Loading data from provided paths')
            adata_list = []
            for idx, data_path in enumerate(self.data_paths):
                self.logger.info('Loading data from {}'.format(data_path))
            
                if os.path.isdir(data_path):
                    file_type = '10x_mtx'
                elif data_path.endswith('.h5ad'):
                    file_type = 'h5ad'
                elif data_path.endswith('txt.gz'):
                    file_type = 'txt_gz'
                elif data_path.endswith('.tsv'):
                    file_type = 'tsv'
                elif data_path.endswith('.csv.gz'):
                    file_type = 'csv_gz'
                elif data_path.endswith('.mtx.gz'):
                    file_type = 'mtx_gz'
                else:
                    self.logger.error(f'Unsupported data format for {data_path}')
                    sys.exit(1)
            
                if file_type == '10x_mtx':
                    adata = sc.read_10x_mtx(
                        data_path,
                        var_names='gene_symbols',
                        cache=True
                    )
                elif file_type == 'h5ad':
                    adata = sc.read(data_path)
                elif file_type == 'txt_gz':
                    adata = sc.read_csv(
                        data_path, 
                        first_column_names=True, 
                        delimiter='\t'
                    ).T
                elif file_type == 'tsv':
                    adata = sc.read_csv(
                        data_path, 
                        first_column_names=True,
                        delimiter='\t'
                    )
                elif file_type == 'csv_gz':
                    adata = sc.read_csv(
                        data_path, 
                        first_column_names=True
                    ).T
                elif file_type == 'mtx_gz':
                    if self.genes_file is None or self.meta_file is None:
                        self.logger.error('genes_file and meta_file must be specified for mtx.gz format.')
                        continue

                    self.logger.info(f'Reading mtx.gz file from {data_path}')

                    try:
                        import gzip
                        from scipy import io
                        import anndata as ad
                        with gzip.open(data_path, 'rb') as f:
                            X = io.mmread(f).tocsr()

                        genes_df = pd.read_csv(self.genes_file, compression='gzip')
                        cells_df = pd.read_csv(self.meta_file, compression='gzip')
                        genes = genes_df['gene_name'].tolist()
                        barcodes = cells_df['bc_wells'].tolist()
                   
                        if X.shape[0] == len(genes) and X.shape[1] == len(barcodes):
                            adata = ad.AnnData(X=X.T, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
                        elif X.shape[0] == len(barcodes) and X.shape[1] == len(genes):
                            adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
                        else:
                            self.logger.error(f'Dimension mismatch: X{X.shape}, genes({len(genes)}), barcodes({len(barcodes)})')
                            adata = None
                            continue
                    except Exception as inner_e:
                        self.logger.error(f'Exception reading mtx.gz files: {inner_e}')
                        adata = None
                else:
                    self.logger.error(f'Unsupported data format for {data_path}')
                    sys.exit(1)
                if adata is None:
                    self.logger.error(f'Failed to load data from {data_path}')
                    continue
                if adata.n_obs == 0:
                    self.logger.warning(f'Dataset loaded from {data_path} is empty (no cells)')
                    continue
                
                if not 'sample' in adata.obs.columns:
                    sample_name = os.path.basename(os.path.normpath(data_path))
                    adata.obs['sample'] = sample_name
                adata_list.append(adata)
            if not adata_list:
                self.logger.error('No datasets loaded successfully. Exiting pipeline.')
                sys.exit(1)
            self.logger.info('All datasets loaded successfully')
            self.adata_list = adata_list
            if len(self.adata_list) == 1:
                self.adata = adata_list[0]
            for idx, adata in enumerate(self.adata_list):
                self.logger.info(f'Dataset {idx} contains {adata.n_obs} cells and {adata.n_vars} genes')
        except Exception as e:
            self.logger.error('Error loading data: {}'.format(e))
            sys.exit(1)
        
    def quality_control(self):
        """_summary_
            perform quality control on the data
        """
        self.logger.info('Starting quality control')
        try:
            qc_adata_list = []
            for adata in self.adata_list:
                adata = adata.copy()
                adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
                if 'batch' in adata.obs.columns:
                    batch_name = adata.obs['batch'].unique()[0]
                else:
                    batch_name = 'unknown_batch'
                # Mitochondrial gene ratio
                sc.pp.calculate_qc_metrics(
                    adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
                )
                plt.figure(figsize=(10, 8))
                plt.hist(adata.obs['n_genes_by_counts'], bins=50, color='skyblue', edgecolor='black')
                
                plt.axvline(x=200, color='red', linestyle='--', linewidth=2, label='200 genes')
                plt.axvline(x=5000, color='darkred', linestyle='--', linewidth=2, label='5000 genes')
                
                plt.xlabel('Number of genes per cell', fontsize=14)
                plt.ylabel('Number of cells', fontsize=14)
                plt.title(f'Distribution of Genes per Cell', fontsize=16)
                plt.xlim(0, np.percentile(adata.obs['n_genes_by_counts'], 100))
                plt.legend(fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, f"n_genes_by_counts_distribution.png"))
                plt.close()

                plt.figure(figsize=(10, 7))
                plt.hist(adata.obs['pct_counts_mt'], bins=50, color='lightgreen', edgecolor='black')
                
                plt.axvline(x=10, color='purple', linestyle='--', linewidth=2, label='10% MT')
                
                plt.xlabel('Mitochondrial gene percentage (%)', fontsize=14)
                plt.ylabel('Number of cells', fontsize=14)
                plt.title(f'Distribution of Mitochondrial Gene Percentage', fontsize=16)
                
                plt.xlim(0, min(adata.obs['pct_counts_mt'].max(), 100))
                
                plt.legend(fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, f"pct_counts_mt_distribution.png"))
                plt.close()

                # Filtering
                initial_cell_count = adata.n_obs
                cells_below_200_genes = np.sum(adata.obs['n_genes_by_counts'] < 200)
                cells_above_5000_genes = np.sum(adata.obs['n_genes_by_counts'] > 5000)
                cells_high_mt = np.sum(adata.obs['pct_counts_mt'] > 10)
                
                adata = adata[
                    (adata.obs['n_genes_by_counts'] > 200) &
                    (adata.obs['n_genes_by_counts'] < 5000) &
                    (adata.obs['pct_counts_mt'] < 10),
                    :
                ]
                sc.pp.filter_genes(adata, min_cells=3)    
                sce.pp.scrublet(adata, batch_key='sample') # doublet detection
                qc_adata_list.append(adata)
                filtered_cell_count = adata.n_obs
                total_removed_cells = initial_cell_count -  filtered_cell_count

                self.logger.info(f"Cells before QC: {initial_cell_count}")
                self.logger.info(f"Cells with genes < 200: {cells_below_200_genes}")
                self.logger.info(f"Cells with genes > 5000: {cells_above_5000_genes}")
                self.logger.info(f"Cells with MT > 10%: {cells_high_mt}")
                self.logger.info(f"Total cells removed by QC {total_removed_cells}")
                self.logger.info(f"Cells after QC: {filtered_cell_count}")
                self.logger.info(f"Quality control completed for batch {batch_name}")
            if len(self.adata_list) == 1:
                self.adata = adata
            self.adata_list = qc_adata_list
        except Exception as e:
            self.logger.error('Error during quality control: {}'.format(e))
            sys.exit(1)
    
    def integrate_data(self):
        """_summary_
            Integrate the QCed datasets.
        """
        self.logger.info('Starting data integration')
        try:
            all_var_names = [set(adata.var_names) for adata in self.adata_list]
            total_vars = set.union(*all_var_names)
            total_var_count = len(total_vars)
            self.logger.info(f'Total number of unique genes across all datasets: {total_var_count}')
            
            common_var_names = set.intersection(*all_var_names)
            common_var_count = len(common_var_names)
            self.logger.info(f'Number of common genes: {common_var_count}')
            
            excluded_var_count = total_var_count - common_var_count
            self.logger.info(f'Total number of excluded genes due to mismatched names: {excluded_var_count}')
            
            for i, adata in enumerate(self.adata_list):
                original_var_count = adata.n_vars
                adata = adata[:, list(common_var_names)]
                new_var_count = adata.n_vars
                excluded_vars = original_var_count - new_var_count
                batch_name = adata.obs['batch'][0] if 'batch' in adata.obs.columns else f'Dataset_{i}'
                self.logger.info(f'Dataset {i} ({batch_name}): Excluded {excluded_vars} genes due to mismatched names')
                self.adata_list[i] = adata
            self.adata = self.adata_list[0].concatenate(*self.adata_list[1:], batch_key='batch')
            self.logger.info('Data integration completed')
            output_file = os.path.join(self.output_dir, 'integrated_data.h5ad')
            self.adata.write(output_file)
            self.logger.info(f'Integrated data saved to {output_file}')
        except Exception as e:
            self.logger.error('Error during data integration: {}'.format(e))
            sys.exit(1)
        
    def preprocess_data(self):
        """_summary_
            Normalization and scaling of the data
        """
        self.logger.info('Starting data preprocessing')
        try:
            processed_adata_list = []
            for adata in self.adata_list:
                adata = adata.copy()
                # Normalization
                sc.pp.normalize_total(adata, target_sum=1e4)
                # Logarithmize the data
                sc.pp.log1p(adata)

                # find highly variable genes
                sc.pp.highly_variable_genes(adata, n_top_genes=2000)
                adata.raw = adata
                processed_adata_list.append(adata)
            self.adata_list = processed_adata_list
            if len(self.adata_list) == 1:
                self.adata = adata
        except Exception as e:
            self.logger.error('Error during data preprocessing: {}'.format(e))
            sys.exit(1)
    
    # def annotate_cell_types(self):
    #     """_summary_
    #         Annotate cell types based on Celltypist
    #     """
    #     self.logger.info('Starting cell type annotation')
    #     try:
    #         celltypist.models.download_models()
    #         # TODO: add option to choose model
    #         model = celltypist.models.Model.load(model='Immune_All_High.pkl')
    #         predictions = celltypist.annotate(
    #             self.adata,
    #             model=model,
    #             majority_voting=True,
    #         )
    #         self.adata.obs['cell_type'] = predictions.predicted_labels.predicted_labels
    #     except Exception as e:
    #         self.logger.error('Error during cell type annotation: {}'.format(e))
    #         sys.exit(1)
    
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

    def clustering(self):
        """_summary_
            Perform finding neighbors and clustering on the data
        """
        self.logger.info('Starting clustering')
        try:
            sc.pp.neighbors(self.adata, n_pcs=self.n_dims)
            sc.tl.leiden(self.adata, resolution=0.5, key_added='leiden')
        except Exception as e:
            self.logger.error('Error during clustering: {}'.format(e))
            sys.exit(1)
        
    def annotate_cell_types(
        self,
        marker_genes: Dict[str, List[str]],
        score_perceontile: float = 50,
        min_markers: int = 3,
        min_delta: float = 0.04,
        obs_key: str = 'leiden',
        out_obs_key : str = 'cell_type',
        out_uns_key: str = 'cluster_debug_info'
        
    ):
        """_summary_
            Assign cell types to clusters based on their average expression and the provided marker_genes dictionary.

            Args:
                param marker_genes: Dict[str, List[str]] mapping each cell type to its list of marker genes
                param score_percentile: Percentile at which to set the best_score threshold (default: 50)
                param min_markers: Minimum number of matched marker genes required for assignment
                param min_delta: Minimum gap required between the top and second-best score for confident assignment
       """
        self.logger.info('Starting cluster-based marker annotation')
        
        if self.adata.raw is None:
            raise ValueError("adata.raw is Empty")
        raw_df = self.adata.raw.to_adata().to_df()
        
        clusters = self.adata.obs[obs_key]
        cluster_avg = raw_df.groupby(clusters).mean()
        
        initial_debug = {}
        for cl in cluster_avg.index:
            scores = {}
            details = {}
            for ct, genes in marker_genes.items():
                matched = [g for g in genes if g in cluster_avg.columns]
                ignored = [g for g in genes if g not in cluster_avg.columns]
                score = cluster_avg.loc[cl, matched].mean() if matched else 0
                scores[ct] = score
                details[ct] = (matched, ignored)
            best_ct = max(scores, key=scores.get)
            initial_debug[cl] = {
                'best_cell_type': best_ct,
                'best_score': scores[best_ct],
                'all_scores': scores,
                'all_details': details,
                'matched': details[best_ct][0],
                'ignored': details[best_ct][1]
            }
            
            all_best_scores = [v['best_score'] for v in initial_debug.values()]
            min_score = np.percentile(all_best_scores, score_perceontile)
            
            final_assignment = {}
            final_debug = {}
            for cl, info in initial_debug.items():
                best = info['best_score']
                sorted_scores = sorted(info['all_scores'].values(), reverse=True)
                second = sorted_scores[1] if len(sorted_scores) > 1 else 0
                gap = best - second
                n_matched = len(info['matched'])
                
                if best < min_score or n_matched < min_markers or gap < min_delta:
                    assigned = 'unknown'
                else:
                    assigned = info['best_cell_type']
                final_assignment[cl] = assigned
                final_debug[cl] = {
                    **info,
                    'assigned': assigned,
                    'score_gap': gap,
                    'n_matched': n_matched
                }
                self.logger.info(
                    "[Cluster %s] assigned -> %s (best=%.3f, second=%.3f, gap=%.3f), matched=%d/%d",
                    cl,
                    assigned,
                    best,
                    second,
                    gap,
                    n_matched,
                    len(marker_genes.get(assigned, []))
                )
                
            self.adata.obs[out_obs_key] = clusters.map(final_assignment)
            self.adata.uns[out_uns_key] = final_debug
            
            self.logger.info("Completed cell type annotation")
            
    # TODO: modify
    def run_tsne(
        self,
        color_by: str,
        cmap_name: str = "Reds"
    ):
        """
        
        """
        try:
            self.logger.info(f'Starting tSNE colored by {color_by}')
            if 'X_tsne' not in self.adata.obsm:
                self.logger.info("Computing t-SNE embedding")
                sc.tl.tsne(self.adata, n_pcs=self.n_dims)
            else:
                self.logger.info("t-SNE already computed, skipping")

            base = plt.get_cmap(cmap_name)
            colors = ["lightgrey", base(0.4), base(1.0)]
            custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                f"grey{cmap_name}", colors, N=256
            )

            title = f"t-SNE plot colored by {color_by}"
            sc.pl.tsne(
                self.adata,
                color=color_by,
                cmap=custom_cmap,
                size=5,
                title=title,
                legend_loc='right margin',
                save=f"_{color_by}_tSNE.png"
            )
            self.logger.info("t-SNE plot complete")
       
        except Exception as e:
            self.logger.error('Error during tSNE: {}'.format(e))
            sys.exit(1)  
   
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
    
    # TODO: need to review
    def run_pipeline(self):
        """_summary_
            Run the entire pipeline based on the steps provided by the user
        """
        self.load_data()

        if 'quality_control' in self.steps:
            self.quality_control()
        else:
            self.logger.info('Skipping quality control step')
        
        if 'preprocessing' in self.steps:
            self.preprocess_data()
        else:
            self.logger.info('Skipping preprocessing step')
            
        if 'pca' in self.steps:
            self.run_pca()
        else:
            self.logger.info('Skipping PCA step')
        
        if 'clustering' in self.steps:
            self.clustering()
        else:
            self.logger.info('Skipping clustering step')
            
        if 'annotation' in self.steps:
            self.annotate_cell_types()
        else:
            self.logger.info('Skipping cell type annotation step')     
            
        self.save_results()
        self.logger.info('Pipeline execution completed successfully')

#TODO: delete
# def init_process(log_queue):
#     logger = logging.getLogger('ScRNAseqPipeline')
#     logger.handlers = []
#     logger.setLevel(logging.INFO)
#     queue_handler = QueueHandler(log_queue)
#     logger.addHandler(queue_handler)  
    
def parse_arguments():
    parser = argparse.ArgumentParser(description="Single cell RNA-seq pipeline")
    parser.add_argument('--data_paths', type=str, nargs='+', required=True, help='Path to Cell Ranger output directories (filtered_feature_bc_matrix) or h5ad data formats')
    parser.add_argument('--genes_file', type=str, default=None, help='Path to genes file (required for mtx.gz)')
    parser.add_argument('--meta_file', type=str, default=None, help='Path to metadata file (required for mtx.gz)')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save outputs')
    parser.add_argument('--gene', type=str, help='Gene of interest to compute correlation')
    parser.add_argument('--target_gene', type=str, help='Target gene of interest')
    parser.add_argument('--reference_gene', type=str, help='Reference gene to compare with')
    parser.add_argument('--cell_types', type=str, nargs='*', help='List of cell types to analyze (e.g., T_cells B_cells)')
    parser.add_argument('--n_dims', type=int, default=4, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    
    parser.add_argument('--only_highly_variable_genes', action='store_true', help='Using only highly variable genes in all analysis steps')
    parser.add_argument('--steps', type=str, nargs='*', help='List of steps to run in the pipeline')
    
    parser.add_argument('--skip_correlation', action='store_true', help='Skip the correlation computation step')
    args = parser.parse_args()
    for path in args.data_paths:
        if path.endswith('.mtx.gz'):
            if args.genes_file is None or args.meta_file is None:
                parser.error('--genes_file and --meta_file are required when using .mtx.gz format.')
        
    return args
    
def main():
    args = parse_arguments()
       
    pipeline = ScRNAseqPipeline(
        data_paths=args.data_paths,
        output_dir=args.output_dir,
        gene=args.gene,
        target_gene=args.target_gene,
        cell_types=args.cell_types,
        n_dims=args.n_dims,
        random_state=args.random_state,
        only_highly_variable_genes=args.only_highly_variable_genes,
        steps=args.steps,
    )
    pipeline.run_pipeline()
    
if __name__ == "__main__":
    main()
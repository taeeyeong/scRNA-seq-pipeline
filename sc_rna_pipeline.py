# TODO: sphinx setting
import os
import sys
import argparse
import logging
import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from typing import Dict, List

import celltypist


class ScRNAseqPipeline:
    def __init__(
        self, data_paths, output_dir, genes_file=None, meta_file=None,
        n_dims=4, random_state=42, only_highly_variable_genes=False, steps=[], marker_genes={}
        ):
        self.data_paths = data_paths
        self.genes_file = genes_file
        self.meta_file = meta_file
        self.output_dir = output_dir
        self.n_dims = n_dims
        self.random_state = random_state
        self.adata = None
        
        self.only_highly_variable_genes = only_highly_variable_genes
        self.steps = steps
        self.marker_genes = marker_genes
        
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
                elif data_path.endswith('.mtx'):
                    file_type = 'mtx'
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
                        
                elif file_type == 'mtx':
                    if self.genes_file is None or self.meta_file is None:
                        self.logger.error('genes_file and meta_file must be specified for mtx format.')
                        continue

                    self.logger.info(f'Reading mtx file from {data_path}')

                    try:
                        from scipy import io
                        import anndata as ad
                        with open(data_path, 'rb') as f:
                            X = io.mmread(f).tocsr()

                        genes_df = pd.read_csv(self.genes_file)
                        cells_df = pd.read_csv(self.meta_file)
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
                        self.logger.error(f'Exception reading mtx files: {inner_e}')
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
            for idx, adata in enumerate(self.adata_list):
                if 'sample' not in adata.obs.columns:
                    self.logger.warning(
                        f"Adata at index {idx} has no 'sample' column in .obs; "
                    )
                    adata.obs['sample'] = f"dataset_{idx}"
                adata = adata.copy()
                adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
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
                plt.savefig(os.path.join(self.output_dir, f"{adata.obs['sample'][0]}_n_genes_by_counts_distribution.png"))
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
                initial_gene_count = adata.n_vars
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
                sce.pp.scrublet(adata, batch_key='sample')  # doublet detection

                if 'predicted_doublet' in adata.obs.columns:
<<<<<<< codex/quality_control-함수-doublet-finding-필터링-확인
                    doublet_cells = adata.obs.index[adata.obs['predicted_doublet']].tolist()
                    doublet_count = len(doublet_cells)
                    if doublet_count:
                        self.logger.info(
                            f"Removed {doublet_count} predicted doublets: {', '.join(doublet_cells)}"
                        )
                    adata = adata[~adata.obs['predicted_doublet']].copy()
                else:
                    self.logger.warning(
                        'Scrublet did not add predicted_doublet column; no doublet filtering applied'
                    )
=======
                    doublet_count = int(adata.obs['predicted_doublet'].sum())
                    adata = adata[~adata.obs['predicted_doublet']].copy()
                    self.logger.info(f"Removed {doublet_count} predicted doublets")
                else:
                    self.logger.warning('Scrublet did not add predicted_doublet column; no doublet filtering applied')
>>>>>>> main

                qc_adata_list.append(adata)
                filtered_cell_count = adata.n_obs
                filtered_gene_count = adata.n_vars
                total_removed_cells = initial_cell_count - filtered_cell_count

                self.logger.info(f"Cells before QC: {initial_cell_count}")
                self.logger.info(f"Genes before QC: {initial_gene_count}")
                self.logger.info(f"Cells with genes < 200: {cells_below_200_genes}")
                self.logger.info(f"Cells with genes > 5000: {cells_above_5000_genes}")
                self.logger.info(f"Cells with MT > 10%: {cells_high_mt}")
                self.logger.info(f"Total cells removed by QC {total_removed_cells}")
                self.logger.info(f"Cells after QC: {filtered_cell_count}")
                self.logger.info(f"Genes after QC: {filtered_gene_count}")
                self.logger.info(f"Quality control completed for sample {adata.obs['sample'][0]}")
            if len(self.adata_list) == 1:
                self.adata = adata
            self.adata_list = qc_adata_list
        except Exception as e:
            self.logger.error('Error during quality control: {}'.format(e))
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
                
                adata.raw = adata
                processed_adata_list.append(adata)
            self.adata_list = processed_adata_list
            if len(self.adata_list) == 1:
                self.adata = adata
        except Exception as e:
            self.logger.error('Error during data preprocessing: {}'.format(e))
            sys.exit(1)
            
    def integrate_data(
        self,
        batch_key: str = 'batch',
        join: str = 'inner'   
    ):
        """_summary_
            Integrate the QCed datasets.
            
            Args:
                batch_key: Name of the obs-column that will record each dataset's batch label
                join: How to join variables across datasets ('inner' to keep intersection, 'outer' to keep union)
        """
        self.logger.info('Starting data integration')
        try:
            for idx, ad in enumerate(self.adata_list):
                if 'sample' not in ad.obs.columns:
                    self.logger.warning(f"Adata {idx} has no 'sample' column in .obs; ")
                    ad.obs['sample'] = f"dataset_{idx}"
                if ad.raw is not None:
                    self.logger.info(f"[Dataset {idx}] clearing raw layer to avoid duplicate-label errors")
                    ad.raw = None

            keys = [
                ad.obs['sample'].unique()[0] for ad in self.adata_list
            ]
            
            
            
            combined = anndata.concat(
                self.adata_list,
                join=join,
                label=batch_key,
                keys=keys,
            )
            self.adata = combined
            self.adata.raw = self.adata
            self.logger.info(
                f"Integration finished: "
                f"{combined.n_obs} cells, {combined.n_vars} genes"
            )
        except Exception as e:
            self.logger.error('Error during data integration: {}'.format(e))
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
            # find highly variable genes
            sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
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
                param score_percentile: Percentile at which to set the best_score threshold (default: 50)
                param min_markers: Minimum number of matched marker genes required for assignment
                param min_delta: Minimum gap required between the top and second-best score for confident assignment
       """
        self.logger.info('Starting cluster-based marker annotation')
        
        if self.adata.raw is None:
            raise ValueError("adata.raw is Empty")
        raw_df = self.adata.raw.to_adata().to_df()
        
        marker_genes = self.marker_genes
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
    
    # TODO: error 수정 
    # @staticmethod
    # def _sanitize_uns(adata):
    #     """
    #     Recursively convert unsupported types in adata.uns to supported ones for h5ad.
    #     """
    #     def sanitize(obj):
    #         if isinstance(obj, dict):
    #             return {str(k): sanitize(v) for k, v in obj.items()}
    #         elif isinstance(obj, (tuple, set)):
    #             return [sanitize(v) for v in obj]
    #         elif isinstance(obj, list):
    #             return [sanitize(v) for v in obj]
    #         elif isinstance(obj, (np.integer, np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64)):
    #             return int(obj)
    #         elif isinstance(obj, (np.floating, np.float_, np.float16, np.float32, np.float64)):
    #             return float(obj)
    #         elif isinstance(obj, np.ndarray):
    #             return obj.tolist()
    #         elif isinstance(obj, pd.DataFrame):
    #             return obj.to_dict(orient='list')
    #         elif isinstance(obj, pd.Series):
    #             return obj.to_list()
    #         elif isinstance(obj, (str, int, float, bool, type(None))):
    #             return obj
    #         else:
    #             return str(obj)  # fallback: convert to string

    #     adata.uns = sanitize(adata.uns)
    #     return adata

    # def save_h5ad(self, file_name: str = None):
    #     """
    #     """
    #     fname = file_name if file_name else "final_adata.h5ad"
    #     out_path = os.path.join(self.output_dir, fname)
    #     try:
    #         self.logger.info(f'Sanitizing .uns before saving to {out_path}')
    #         self.adata = self._sanitize_uns(self.adata)
    #         self.adata.write(out_path)
    #         self.logger.info(f"Saved adata to {out_path}")
    #     except Exception as e:
    #         self.logger.error(f"Error saving adata to h5ad: {e}", exc_info=True)
    #         sys.exit(1)
            
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
        
        if 'integrate_data' in self.steps:
            self.integrate_data()
        else:
            self.logger.info('Skipping integration step')
            
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
        
        if 'save_h5ad' in self.steps:
            self.save_h5ad()
        else:
            self.logger.info('Skipping save h5ad step')
            
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
    parser.add_argument('--n_dims', type=int, default=4, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    
    parser.add_argument('--only_highly_variable_genes', action='store_true', help='Using only highly variable genes in all analysis steps')
    parser.add_argument('--steps', type=str, nargs='*', help='List of steps to run in the pipeline')
    parser.add_argument('--marker_genes', type=Dict, default=None, nargs='*', help=' Dict[str, List[str]] mapping each cell type to its list of marker genes')
    
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
        n_dims=args.n_dims,
        random_state=args.random_state,
        only_highly_variable_genes=args.only_highly_variable_genes,
        steps=args.steps,
        marker_genes = args.marker_genes
    )
    pipeline.run_pipeline()
    
if __name__ == "__main__":
    main()
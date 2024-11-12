# TODO: sphinx setting
import os
import sys
import argparse
import logging
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import celltypist
from multiprocessing import current_process
from scipy.stats import pearsonr
from logging.handlers import QueueHandler

class ScRNAseqPipeline:
    def __init__(
        self, data_path, output_dir, gene, cell_types, target_gene=None, reference_gene=None, 
        n_dims=4, random_state=42, skip_qc=False,only_highly_variable_genes=False, 
        skip_preprocessing=False, skip_batch_correction=False,
        skip_pca=False, skip_clustering=False, skip_visualization=False, skip_annotation=False, 
        skip_correlation=False
        ):
        self.data_path = data_path
        self.output_dir = output_dir
        self.gene = gene
        self.target_gene = target_gene
        self.reference_gene = reference_gene
        self.cell_types = cell_types
        self.n_dims = n_dims
        self.random_state = random_state
        self.adata = None
        
        self.only_highly_variable_genes = only_highly_variable_genes
        self.skip_qc = skip_qc
        self.skip_preprocessing = skip_preprocessing
        self.skip_batch_correction = skip_batch_correction
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
            self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-') | self.adata.var_names.str.startswith('mt-')
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
            

    def batch_effect_correction_combat(self):
        """_summary_
            Perform batch effect correction using ComBat.
        """
        self.logger.info('Starting batch effect correction using ComBat')
        try:
            if 'batch' not in self.adata.obs.columns:
                self.logger.error('Batch information not found in adata.obs')
                sys.exit(1)
            
            # apply ComBat to the log1p transformed data
            sc.pp.combat(self.adata, key='batch')
            self.logger.info('Combat batch effect correction completed')
            
        except Exception as e:
            self.logger.error(f'Error during batch effect correction: {e}')
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
        self.logger.info(f'Starting tSNE colored by {color_by}')
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
                save=f'_{color_by}_tSNE.png' 
            )
            self.logger.info(f'tSNE plot saved as _{color_by}_tSNE.png')
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
    def compute_correlation(self, top_n=30, cell_type=None):
        """
        Compute Pearson correlation between the target gene and all other genes.

        Args:
            top_n (int): Number of top correlated genes to return.
            cell_type (str, optional): Specific cell type to compute correlation within. If None, use all cells.
        Returns:
            List[Tuple[str, float]]: List of tuples containing gene names and their correlation coefficients.
        """
        try:
            if cell_type:
                self.logger.info(f'Computing correlation for target gene {self.target_gene} within cell type {cell_type}')
                adata = self.adata[self.adata.obs['predicted_celltype'] == cell_type]
            else:
                self.logger.info(f'Computing correlation for target gene {self.target_gene} across all cells')
                adata = self.adata  # 전체 셀 사용
            if self.target_gene not in adata.var_names:
                self.logger.error(f"Target gene {self.target_gene} not found in the dataset.")
                return []
            if adata.n_obs == 0:
                self.logger.warning('No cells found for the specified cell type.')
                return []
            target_exp = adata[:, self.target_gene].X.toarray().flatten()
            correlations = []
            for gene in adata.var_names:
                gene_exp = adata[:, gene].X.toarray().flatten()
                if np.std(gene_exp) == 0 or np.std(target_exp) == 0:
                    continue
                corr, p_val = pearsonr(target_exp, gene_exp)
                if not np.isnan(corr) and p_val < 0.05:
                    correlations.append((gene, corr))
            if not correlations:
                self.logger.warning('No significant correlated genes found.')
                return []
            # 절대 상관계수 기준으로 정렬 후 상위 top_n 유전자 선택
            top_genes = sorted(correlations, key=lambda x: -abs(x[1]))[:top_n]
            self.logger.info(f'Top {top_n} correlated genes: {top_genes}')
            return top_genes
        except Exception as e:
            self.logger.error(f'Error computing correlation: {e}')
            sys.exit(1)
            
    def compute_correlation_per_cell_type(self, top_n=30):
        """Compute Pearson correlation between the target gene and all other genes within each specified cell type.

        Args:
            top_n (int): Number of top correlated genes to return per cell type.

        Returns:
            Dict[str, List[Tuple[str, float]]]: Dictionary mapping cell types to list of tuples containing gene names and their correlation coefficients.
        """
        try:
            self.logger.info(f'{self.target_gene}와의 상관관계를 각 셀 타입별로 계산합니다.')

            adata = self.adata  # 고변이 유전자 포함 여부에 따른 데이터 사용

            if self.target_gene not in adata.var_names:
                self.logger.error(f"타겟 유전자 {self.target_gene}가 데이터셋에 없습니다.")
                sys.exit(1)

            # 분석할 셀 타입 목록 설정
            if self.cell_types:
                selected_cell_types = [ct for ct in self.cell_types if ct in adata.obs['predicted_celltype'].unique()]
                missing_cell_types = set(self.cell_types) - set(selected_cell_types)
                if missing_cell_types:
                    self.logger.warning(f'다음 셀 타입은 데이터에 없으므로 제외됩니다: {missing_cell_types}')
            else:
                selected_cell_types = adata.obs['cell_type'].unique().tolist()

            self.logger.info(f'분석할 셀 타입: {selected_cell_types}')

            cell_type_gene_corr = {}  # 각 셀 타입별로 유전자와 상관계수 저장

            for cell_type in selected_cell_types:
                self.logger.info(f'셀 타입 처리 중: {cell_type}')
                # 셀 타입 필터링
                cell_indices = adata.obs['predicted_celltype'] == cell_type
                subset = adata[cell_indices]

                if subset.n_obs == 0:
                    self.logger.warning(f'셀 타입 {cell_type}에 해당하는 세포가 없습니다.')
                    continue

                target_exp = subset[:, self.target_gene].X.toarray().flatten()

                correlations = []
                for gene in subset.var_names:
                    if gene == self.target_gene:
                        continue
                    gene_exp = subset[:, gene].X.toarray().flatten()
                    if np.std(gene_exp) == 0 or np.std(target_exp) == 0:
                        continue
                    # 발현 비율 계산
                    expr_ratio = np.mean(gene_exp > 0)
                    if expr_ratio < 0.1:
                        continue  # 발현 비율이 10% 미만인 유전자는 제외
                    corr, p_val = pearsonr(target_exp, gene_exp)
                    if not np.isnan(corr) and p_val < 0.05:
                        correlations.append((gene, corr))

                if not correlations:
                    self.logger.warning(f'셀 타입 {cell_type}에서 유의미한 상관된 유전자가 없습니다.')
                    continue

                # 절대 상관계수 기준으로 정렬 후 상위 top_n 유전자 선택
                top_genes = sorted(correlations, key=lambda x: -abs(x[1]))[:top_n]
                self.logger.info(f'셀 타입 {cell_type}의 상위 {top_n}개 상관된 유전자: {top_genes}')

                cell_type_gene_corr[cell_type] = top_genes

            return cell_type_gene_corr

        except Exception as e:
            self.logger.error(f'셀 타입별 상관관계 계산 중 오류 발생: {e}')
            sys.exit(1)
        # def compute_correlation(self, top_n=30):
        #     """_summary_
    #         Compute correlation between gene of interest and all other genes

    #         Args:
    #             top_n (int): Number of top correlated genes to return
            
    #         Returns:
    #             List[Tuple[str, float, float]]: List of tuples 
    #                 containing gene name, correlation coefficients, and p-value
    #     """
    #     try: 
    #         self.logger.info(f'Process {current_process().name}: Computing correlation between {self.gene} and all other genes')
    #         if self.gene not in self.adata.var_names:
    #             self.logger.error(f"Gene {self.gene} not found in the dataset.")
    #             sys.exit(1)
                
    #         gene_exp = self.adata[:, self.gene].X.toarray().flatten()
    #         all_genes = self.adata.var_names
            
    #         if self.cell_types:
    #             selected_cell_types = [ct for ct in self.cell_types if ct in self.adata.obs['predicted_celltype'].unique()]
    #             missing_cel_types = set(self.cell_types) - set(selected_cell_types)
    #             if missing_cel_types:
    #                 self.logger.warning(f'The following specified cell types were not found in the data and will be ignored: {missing_cel_types}')
    #             top_genes_set = set()
    #             for cell_type in selected_cell_types:
    #                 self.logger.info(f'Processing cell type: {cell_type}')
    #                 cell_indices = self.adata.obs['predicted_celltype'] == cell_type
    #                 subset = self.adata[cell_indices]
    #                 if subset.n_obs == 0:
    #                     self.logger.warning(f'No cells found for cell type: {cell_type}')
    #                     continue
    #                 subset_target_exp = subset[:, self.target_gene].X.toarray().flatten()
    #                 correlations = []
    #                 for gene in subset.var_names:
    #                     gene_exp = subset[:, gene].X.toarray().flatten()
    #                     if np.std(gene_exp) == 0 or np.std(subset_target_exp) == 0:
    #                         continue
    #                     corr, p_val = pearsonr(subset_target_exp, gene_exp)
    #                     correlations.append((gene, corr, p_val))
    #                 significant_genes = [item for item in correlations if item[2] < 0.05]
    #                 if not significant_genes:
    #                     self.logger.warning(f'No significant correlated genes found for cell type: {cell_type}')
    #                     continue
    #                 top_genes = sorted(significant_genes, key=lambda x: -abs(x[1]))[:top_n]
    #                 self.logger.info(f'Top {top_n} correlated genes for cell type {cell_type}: {[gene for gene, _, _ in top_genes]}')
    #                 for gene, _, _  in top_genes:
    #                     top_genes_set.add(gene)
    #                 if not top_genes_set:
    #                     self.logger.warning('No significant correlated genes found across all specified cell types.')
    #                     return []
    #                 self.logger.info(f'Final list of top correlated genes across all cell types: {sorted(top_genes_set)}')
    #                 return sorted(top_genes_set)
    #         else:
    #             correlations = []
    #             for gene in all_genes:
    #                 gene_data = self.adata[:, gene].X.toarray().flatten()
    #                 corr, p_value = pearsonr(gene_exp, gene_data)
    #                 correlations.append((gene, corr, p_value))
    #                 print(f'{gene}: corr={corr}, p-value={p_value}')
    #             significant_genes = [item for item in correlations if item[2] < 0.05]
    #             top_genes = sorted(significant_genes, key=lambda x: -abs(x[1]))[:top_n]
    #             self.logger.info(f'Top {top_n} correlated genes computed.')
    #             return top_genes
        
    #     except Exception as e:
    #         self.logger.error(f'Error computing correlation: {e}')
    #         sys.exit(1)
                
    def _plot_dotplot(self, percentage_expr, average_expr, genes, cell_types):
        """_summary_
            Plot a dot plot of percentage expression and average expression

        Args:
            percentage_expr (pd.DataFrame): DataFrame of percentage expression.
            average_expr (pd.DataFrame): DataFrame of average expression.
            genes (List[str]): List of gene names.
            cell_types (List[str]): List of cell types. 
        """
        try:
            sns.set(style='whitegrid')
            # data transformation
            percentage_expr = percentage_expr.astype(float)
            average_expr = average_expr.astype(float)
            # set graph size
            fig, ax = plt.subplots(figsize=(len(cell_types) * 1.2, len(genes) * 0.8))
            # color
            norm = plt.Normalize(vmin=average_expr.values.min(), vmax=average_expr.values.max())

            min_size = 300
            max_size = 3000

            pct_values = percentage_expr.values.flatten()
            pct_min = pct_values.min()
            pct_max = pct_values.max()

            if pct_max == pct_min:
                pct_max += 1
            x_positions = np.arange(len(cell_types))
            y_positions = np.arange(len(genes))

            for i, gene in enumerate(genes):
                for j, cell_type in enumerate(cell_types):
                    pct = percentage_expr.loc[gene, cell_type]
                    size = min_size + ((pct - pct_min) / (pct_max - pct_min)) * (max_size - min_size)
                    color = average_expr.loc[gene, cell_type]
                    ax.scatter(x_positions[j], y_positions[i], s=size, c=[color], cmap='Reds', norm=norm, edgecolors='gray')
            ax.set_xticks(x_positions)
            ax.set_xticklabels(cell_types, rotation=90, fontsize=18)
            ax.set_yticks(y_positions)
            ax.set_yticklabels(genes, fontsize=18)

            sm = plt.cm.ScalarMappable(cmap='Reds', norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label('Average Expression', fontsize=20)
            ax.set_xlabel('Cell Type', fontsize=18)
            ax.set_ylabel('Gene', fontsize=18)

            legend_percentages = np.linspace(pct_min, pct_max, num=5)
            legend_handles = []
            for pct in legend_percentages:
                size = min_size + ((pct - pct_min) / (pct_max - pct_min)) * (max_size - min_size)
                handle = mlines.Line2D([], [], color='gray', marker='o', linestyle='None',
                                       markersize=np.sqrt(size / np.pi), label=f'{pct:.1f}%')
                legend_handles.append(handle)
            ax.legend(handles=legend_handles, title='% of Cells Expressing', loc='upper right', bbox_to_anchor=(1.2, 1))

            plt.tight_layout()
            # plot save
            plot_filename = os.path.join(self.output_dir, f'{self.gene}_correlated_genes_dotplot.png')
            plt.savefig(plot_filename)
            plt.close(fig)
            self.logger.info(f'Correlation dot plot saved as {plot_filename}')
        except Exception as e:
            self.logger.error(f'Error generating dot plot: {e}')
            sys.exit(1)
        
    def visualize_correlated_genes(self, top_n=30):
        """Visualize the top correlated genes with the target gene across specified cell types.

        Args:
            top_n (int): Number of top correlated genes to visualize.
        """
        try:
            self.logger.info(f'Process {current_process().name}: Visualizing the top {top_n} highly correlated genes with {self.target_gene}')
            
            top_genes = self.compute_correlation(top_n=top_n)
            if not top_genes:
                self.logger.warning('No significant correlated genes found.')
                return
            
            # 유전자 이름과 상관계수를 분리
            top_gene_names = [gene for gene, _ in top_genes]
            gene_correlations = {gene: corr for gene, corr in top_genes}
            
            # 분석할 셀타입 목록을 설정
            if self.cell_types:
                selected_cell_types = [ct for ct in self.cell_types if ct in self.adata.obs['predicted_celltype'].unique()]
                missing_cell_types = set(self.cell_types) - set(selected_cell_types)
                if missing_cell_types:
                    self.logger.warning(f'The following specified cell types were not found in the data and will be ignored: {missing_cell_types}')
            else:
                selected_cell_types = self.adata.obs['predicted_celltype'].unique().tolist()
            
            self.logger.info(f'Analyzing the following cell types: {selected_cell_types}')
            
            # 선택한 셀타입의 데이터를 필터링합니다.
            adata_filtered = self.adata[self.adata.obs['predicted_celltype'].isin(selected_cell_types)]
            
            # 도트 플롯 저장 경로 설정
            plot_filename = f'{self.target_gene}_correlated_genes_dotplot.png'
            sc.settings.figdir = self.output_dir  # 결과 저장 디렉토리 설정
            
            # 도트 플롯 생성 및 저장
            sc.pl.dotplot(
                adata_filtered,
                var_names=top_gene_names,
                groupby='predicted_celltype',
                standard_scale='var',  # 각 유전자를 0에서 1로 정규화
                var_group_labels=[''] * len(top_gene_names),  # 그룹 레이블 제거
                var_group_positions=[(i, i) for i in range(len(top_gene_names))],  # 각 유전자를 개별 그룹으로 처리
                swap_axes=True,  # 축 교환 (유전자가 y축)
                dendrogram=False,
                save=f'_{plot_filename}',
                show=False
            )
            self.logger.info(f'Correlation dot plot saved as {plot_filename}')
        
        except Exception as e:
            self.logger.error(f'Error visualizing correlated genes: {e}')
            sys.exit(1)
            
    def visualize_correlated_genes_combined(self, top_n=30):
  
        """Visualize the top correlated genes with the target gene for each cell type,
        including all selected cell types in the plot.

        Args:
            top_n (int): Number of top correlated genes to visualize per cell type.
        """
        try:
            self.logger.info(f'Visualizing the top {top_n} correlated genes with {self.target_gene} for each cell type, including all selected cell types')

            cell_type_gene_corr = self.compute_correlation_per_cell_type(top_n=top_n)

            if not cell_type_gene_corr:
                self.logger.warning('No significant correlated genes found for any cell type.')
                return

            # 선택한 셀 타입의 데이터만 필터링
            if self.cell_types:
                adata_filtered = self.adata[self.adata.obs['predicted_celltype'].isin(self.cell_types)]
            else:
                adata_filtered = self.adata

            for cell_type, gene_corr_list in cell_type_gene_corr.items():
                if not gene_corr_list:
                    continue

                # 유전자 이름 추출
                top_gene_names = [gene for gene, _ in gene_corr_list]

                # 유전자 목록을 30개씩 나누기
                gene_chunks = [top_gene_names[i:i + 30] for i in range(0, len(top_gene_names), 30)]

                for idx, gene_chunk in enumerate(gene_chunks):
                    # 도트 플롯 저장 경로 설정
                    plot_filename = f'{self.target_gene}_correlated_genes_dotplot_{cell_type}_part{idx+1}.png'
                    sc.settings.figdir = self.output_dir  # 결과 저장 디렉토리 설정

                    # 도트 플롯 생성 및 저장
                    sc.pl.dotplot(
                        adata_filtered,
                        var_names=gene_chunk,
                        groupby='predicted_celltype',
                        standard_scale='var',  # 각 유전자를 0에서 1로 정규화
                        swap_axes=True,
                        dendrogram=False,
                        show=False,
                        save=f'_{plot_filename}'
                    )
                    self.logger.info(f'Correlation dot plot for genes from cell type {cell_type} part {idx+1} saved as {plot_filename}')

        except Exception as e:
            self.logger.error(f'Error visualizing correlated genes per cell type combined: {e}')
            sys.exit(1)
    # def visualize_correlated_genes(self, top_n=30):
    #     try:
    #         self.logger.info(f'Process {current_process().name}: Visualize the top {top_n} highly correlated genes with {self.gene}')
            
    #         top_genes = self.compute_correlation(top_n=top_n)
    #         if not top_genes:
    #             self.logger.warning('No significant correlated genes found.')
    #             return
            
    #         if self.cell_types:
    #             top_gene_names = top_genes
    #         else:
    #             top_gene_names = [item[0] for item in top_genes]
            
    #         if 'cell_type' in self.adata.obs.columns:
    #             cell_type_col = 'cell_type'
    #         elif 'predicted_celltype' in self.adata.obs.columns:
    #             cell_type_col = 'predicted_celltype'
    #         else:
    #             self.logger.error('Cell type annotations not found in adata.obs')
    #             sys.exit(1)

    #         if self.cell_types:
    #             selected_cell_types = [ct for ct in self.cell_types if ct in self.adata.obs[cell_type_col].unique()]
    #             missing_cell_types = set(self.cell_types) - set(selected_cell_types)
    #             if missing_cell_types:
    #                 self.logger.warning(f'The following specified cell types were not found in the data and will be ignored: {missing_cell_types}')
    #         else:
    #             selected_cell_types = self.adata.obs[cell_type_col].unique()
    #         self.logger.info(f'Analyzing the following cell types: {selected_cell_types}')
            
    #         percentage_expr = pd.DataFrame(index=top_gene_names, columns=selected_cell_types)
    #         average_expr = pd.DataFrame(index=top_gene_names, columns=selected_cell_types)
            
    #         for gene in top_gene_names:
    #             for cell_type in selected_cell_types:
    #                 cell_indices = self.adata.obs[cell_type_col] == cell_type
    #                 expr_values = self.adata[cell_indices, gene].X.toarray().flatten()
    #                 # 발현된 세포의 비율 계산
    #                 pct = np.sum(expr_values > 0) / len(expr_values) * 100 if len(expr_values) > 0 else 0
    #                 percentage_expr.loc[gene, cell_type] = pct
    #                 # 발현된 세포들 중 평균 발현량 계산
    #                 if np.sum(expr_values > 0) > 0:
    #                     avg = np.mean(expr_values[expr_values > 0])
    #                 else:
    #                     avg = 0
    #                 average_expr.loc[gene, cell_type] = avg
    #         self.logger.info('Expression percentage and average expression calculated for dot plot.')
    #         # 도트 플롯 생성
    #         self._plot_dotplot(percentage_expr, average_expr, top_gene_names, selected_cell_types)


    #     except Exception as e:
    #         self.logger.error(f'Error visualizing correlated genes: {e}')  
             
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
        
        if not self.skip_batch_correction:
            self.batch_effect_correction_combat()
        else:
            self.logger.info('Skipping batch correction step')
            
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
        self.run_tsne(color_by='METTL3_status') # TODO: 수정
        # self.run_tsne(color_by='AR_METTL3_combined_status') # TODO: 수정
        # self.run_tsne(color_by='predicted_celltype')
        # self.run_tsne(color_by='batch')
        # self.run_tsne(color_by='total_counts')

        
        if not self.skip_visualizaton:
            self.visualize_clusters()
        else:
            self.logger.info('Skipping visualization step')
        
        if not self.skip_correlation:
            self.compute_correlation()
        else:
            self.logger.info('Skipping correlation computation step')
        if not self.cell_types:
            self.visualize_correlated_genes()
        else:
            self.visualize_correlated_genes_combined()
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
    parser.add_argument('--target_gene', type=str, help='Target gene of interest')
    parser.add_argument('--reference_gene', type=str, help='Reference gene to compare with')
    parser.add_argument('--cell_types', type=str, nargs='*', help='List of cell types to analyze (e.g., T_cells B_cells)')
    parser.add_argument('--n_dims', type=int, default=4, help='Number of principal components to use')
    parser.add_argument('--random_state', type=int, default=42, help='Random state for reproducibility')
    
    parser.add_argument('--only_highly_variable_genes', action='store_true', help='Using only highly variable genes in all analysis steps')
    parser.add_argument('--skip_qc', action='store_true', help='Skip the quality control step')
    parser.add_argument('--skip_preprocessing', action='store_true', help='Skip the preprocessing step')
    parser.add_argument('--skip_batch_correction', action='store_true', help='Skip the batch correction step')
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
        target_gene=args.target_gene,
        cell_types=args.cell_types,
        n_dims=args.n_dims,
        random_state=args.random_state,
        only_highly_variable_genes=args.only_highly_variable_genes,
        skip_qc=args.skip_qc,
        skip_preprocessing=args.skip_preprocessing,
        skip_batch_correction=args.skip_batch_correction,
        skip_pca=args.skip_pca,
        skip_clustering=args.skip_clustering,
        skip_visualization=args.skip_visualization,
        skip_annotation=args.skip_annotation,
        skip_correlation=args.skip_correlation
    )
    pipeline.run_pipeline()
    
if __name__ == "__main__":
    main()
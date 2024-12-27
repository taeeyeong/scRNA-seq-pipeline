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
        self, data_paths, output_dir, gene, cell_types, target_gene=None, reference_gene=None, 
        n_dims=4, random_state=42, only_highly_variable_genes=False, steps=[]
        ):
        self.data_paths = data_paths
        self.output_dir = output_dir
        self.gene = gene
        self.target_gene = target_gene
        self.reference_gene = reference_gene
        self.cell_types = cell_types
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
                else:
                    self.logger.error(f'Unsupported data format for {data_path}')
                    sys.exit(1)
                if adata is None:
                    self.logger.error(f'Failed to load data from {data_path}')
                    continue
                if adata.n_obs == 0:
                    self.logger.warning(f'Dataset loaded from {data_path} is empty (no cells)')
                    continue
                
                if not 'batch' in adata.obs.columns:
                    batch_name = str(idx)
                    adata.obs['batch'] = batch_name
                adata_list.append(adata)
                self.logger.info(f'Dataset {idx} loaded with batch name') # TODO: Batch name 추가
            if not adata_list:
                self.logger.error('No datasets loaded successfully. Exiting pipeline.')
                sys.exit(1)
            self.logger.info('All datasets loaded successfully')
            self.adata_list = adata_list
            
            # if self.only_highly_variable_genes:
            #     tg_mask = self.adata.var_names.isin([self.gene])
            #     hvg_mask = self.adata.var.highly_variable
            #     combined_mask = tg_mask | hvg_mask 
            #     self.adata = self.adata[:, combined_mask].copy()
                
            for idx, adata in enumerate(self.adata_list):
                self.logger.info(f'Dataset {idx} contains {adata.n_obs} cells and {adata.n_vars} genes')
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
                # visualize quality control metrics
                # sc.pl.violin(
                #     adata, 
                #     ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                #     jitter=0.4, 
                #     multi_panel=True, 
                #     save=f'{batch_name}_qc_violin.png',
                #     show=False
                # )
                # visualize n_genes_by_counts distribution
                sns.histplot(adata.obs['n_genes_by_counts'], bins=50)
                plt.xlabel('n_genes_by_counts')
                plt.ylabel('Number of cells')
                plt.title('Distribution of n_genes_by_counts')
                plt.savefig(os.path.join(self.output_dir, f"{batch_name}_n_genes_by_counts_distribution.png"))
                plt.close()

                # lowest threshold for number of genes and counts
                X = np.percentile(adata.obs['n_genes_by_counts'], 5) # lower threshold: 5%
                Y = np.percentile(adata.obs['total_counts'], 99) # upper threshold: 1%

                self.logger.info(f"Calculated lower threshold (X): {X}")
                self.logger.info(f"Calculated upper threshold (Y): {Y}")

                # visualize pct_counts_mt distribution
                sns.histplot(adata.obs['pct_counts_mt'], bins=50)
                plt.xlabel('pct_counts_mt')
                plt.ylabel('Number of cells')
                plt.title('Distribution of pct_counts_mt')
                plt.savefig(os.path.join(self.output_dir, f"{batch_name}_pct_counts_mt_distribution.png"))
                plt.close()

                Z = np.percentile(adata.obs['pct_counts_mt'], 90) # upper threshold: 90%
                self.logger.info(f"Calculated upper threshold for pct_counts_mt (Z): {Z}")

                # Filtering
                initial_cell_count = adata.n_obs
                adata = adata[adata.obs.n_genes_by_counts > X, :]
                adata = adata[adata.obs.total_counts < Y, :]
                adata = adata[adata.obs.pct_counts_mt < Z, :]
                qc_adata_list.append(adata)
                filtered_cell_count = adata.n_obs

                self.logger.info(f"Cells before QC: {initial_cell_count}")
                self.logger.info(f"Cells after QC: {filtered_cell_count}")
                self.logger.info(f"Quality control completed for batch {batch_name}")
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
            # Normalization
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            # Logarithmize the data
            sc.pp.log1p(self.adata)
            
            # find highly variable genes
            sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
            # self.adata = self.adata[:, self.adata.var['highly_variable']].copy()
            
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
                adata = self.adata[self.adata.obs['cell_type'] == cell_type]
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
                selected_cell_types = [ct for ct in self.cell_types if ct in adata.obs['cell_type'].unique()]
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
                cell_indices = adata.obs['cell_type'] == cell_type
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
                selected_cell_types = [ct for ct in self.cell_types if ct in self.adata.obs['cell_type'].unique()]
                missing_cell_types = set(self.cell_types) - set(selected_cell_types)
                if missing_cell_types:
                    self.logger.warning(f'The following specified cell types were not found in the data and will be ignored: {missing_cell_types}')
            else:
                selected_cell_types = self.adata.obs['cell_type'].unique().tolist()
            
            self.logger.info(f'Analyzing the following cell types: {selected_cell_types}')
            
            # 선택한 셀타입의 데이터를 필터링합니다.
            adata_filtered = self.adata[self.adata.obs['cell_type'].isin(selected_cell_types)]
            
            # 도트 플롯 저장 경로 설정
            plot_filename = f'{self.target_gene}_correlated_genes_dotplot.png'
            sc.settings.figdir = self.output_dir  # 결과 저장 디렉토리 설정
            
            # 도트 플롯 생성 및 저장
            sc.pl.dotplot(
                adata_filtered,
                var_names=top_gene_names,
                groupby='cell_type',
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
                adata_filtered = self.adata[self.adata.obs['cell_type'].isin(self.cell_types)]
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
                        groupby='cell_type',
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

    def _load_gene_list(self, gene_list_file):
        """_summary_
            Helper method to load gene lists from a file

        Args:
            gene_list_file (str): Path to the file containing the gene list
                Supported formats:
                - CSV file (.csv):
                    gene_name,group
                    METTL3,writer
                    METTL14,writer
                    YTHDF1,reader
                    ...
                - Text file (.txt):
                    # Writers
                    METTL3
                    METTL14
                    ...
                    # Readers
                    YTHDF1
                    ...
        Returns:
            dict: Dictionary mapping gene groups to lists of genes
        """
        try:
            genes_dict = {}
            
            if gene_list_file.endswith('.csv'):
                df = pd.read_csv(gene_list_file)
                for group in df['group'].unique():
                    genes_dict[group] = df[df['group'] == group]['gene_name'].tolist()
            
            elif gene_list_file.endswith('.txt'):
                current_group = None
                with open(gene_list_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('#'):
                            current_group = line[1:].strip().lower()
                            genes_dict[current_group] = []
                        elif line and current_group:
                            genes_dict[current_group].append(line)
            else:
                self.logger.error('Unsupported file format. Please use .csv or .txt files.')
                return None
            
            return genes_dict
            
        except Exception as e:
            self.logger.error(f'Error loading gene list: {e}')
            return None
        
    def perform_differential_expression(self, group_key='Group', genes=None, gene_file=None, 
                                 conditions=None, min_cells=3, pval_threshold=0.05):
        """Perform differential expression analysis for specified genes between groups.

        Args:
            group_key (str): Column name in adata.obs for group comparison (default: 'Group')
            genes (List[str], optional): List of genes to analyze. If None, uses default RNA modification genes
            gene_file (str, optional): Path to file containing gene list
            conditions (List[str]): Two conditions to compare [reference, comparison]
            min_cells (int): Minimum number of cells required per group
            pval_threshold (float): P-value threshold for significance

        Returns:
            dict: Dictionary containing DE results for each cell type
        """
        try:
            self.logger.info('Starting differential expression analysis')
            
            # Default RNA modification-related genes
            default_gene_list = [
                'ALKBH8', 'BCDIN3D', 'BMT2', 'BUD23', 'CBLL1', 'CDK5RAP1', 'CDKAL1', 'CEBPZ', 'CMTR1', 'CMTR2', 
                'DIMT1', 'EMG1', 'FBL', 'FBLL1', 'FDXACB1', 'FMR1', 'FTSJ1', 'FTSJ3', 'HENMT1', 'HSD17B10', 'LARP7', 
                'LCMT2', 'MEPCE', 'METTL1', 'METTL14', 'METTL15', 'METTL16', 'METTL2A', 'METTL2B', 'METTL3', 'METTL4', 
                'METTL5', 'METTL6', 'METTL7A', 'METTL7B', 'METTL8', 'MRM1', 'MRM2', 'MRM3', 'MTERF4', 'NOP2', 'NSUN2', 
                'NSUN3', 'NSUN4', 'NSUN5', 'NSUN6', 'NSUN7', 'PCIF1', 'PRORP', 'RAMAC', 'RBM15', 'RBM15B', 'RNGTT', 'RNMT', 
                'RRNAD1', 'RSAD1', 'SPOUT1', 'TARBP1', 'TFB1M', 'TFB2M', 'TGS1', 'THADA', 'THUMPD2', 'THUMPD3', 'TRDMT1',
                'TRIT1', 'TRMO', 'TRMT1', 'TRMT10A', 'TRMT10B', 'TRMT10C', 'TRMT11', 'TRMT112', 'TRMT12', 'TRMT13',
                'TRMT1L', 'TRMT2A', 'TRMT2B', 'TRMT44', 'TRMT5', 'TRMT6', 'TRMT61A', 'TRMT61B', 'TRMT9B', 'TRMU', 'TYW3', 
                'VIRMA', 'WDR4', 'WDR6', 'WTAP', 'ZC3H13', 'ZCCHC4'
            ]
            
            # Check if group column exists
            if group_key not in self.adata.obs.columns:
                self.logger.error(f'Group column "{group_key}" not found in data')
                return None
            
            # Determine genes to analyze
            analysis_genes = genes if genes else default_gene_list
            if gene_file:
                file_genes = self._load_genes_from_file(gene_file)
                if file_genes is None:
                    return None
                analysis_genes = list(set(analysis_genes + file_genes))
            
            # Filter genes that exist in the dataset
            available_genes = [gene for gene in analysis_genes if gene in self.adata.var_names]
            missing_genes = set(analysis_genes) - set(available_genes)
            
            if missing_genes:
                self.logger.warning(
                    f'The following genes were not found in the dataset: {", ".join(missing_genes)}'
                )
            
            if not available_genes:
                self.logger.error('None of the specified genes were found in the dataset.')
                return None
                
            self.logger.info(f'Analyzing {len(available_genes)} genes: {", ".join(available_genes)}')
            
            # Subset AnnData to include only the genes of interest
            adata_subset = self.adata[:, available_genes].copy()
            
            de_results = {}
            
            # Perform DE analysis for each cell type
            for cell_type in adata_subset.obs['cell_type'].unique():
                self.logger.info(f'Analyzing cell type: {cell_type}')
                
                # Subset data for current cell type
                cell_mask = adata_subset.obs['cell_type'] == cell_type
                cell_adata = adata_subset[cell_mask]
                
                # Check minimum cell number requirement for each group
                group_counts = cell_adata.obs[group_key].value_counts()
                if any(group_counts[cond] < min_cells for cond in conditions):
                    self.logger.warning(
                        f'Skipping {cell_type} - insufficient cells in groups: '
                        f'{dict(group_counts[conditions])}'
                    )
                    continue
                
                # Perform DE analysis
                sc.tl.rank_genes_groups(
                    cell_adata, 
                    groupby=group_key,
                    reference=conditions[0],
                    method='wilcoxon',
                    key_added='de_results'
                )
                
                # Get results for all analyzed genes
                de_stats = sc.get.rank_genes_groups_df(
                    cell_adata,
                    group=conditions[1],
                    key='de_results'
                )
                
                # Filter significant results
                significant_de = de_stats[de_stats['pvals_adj'] < pval_threshold]
                
                if not significant_de.empty:
                    de_results[cell_type] = significant_de
                    
                    # Visualize results
                    self._plot_de_results(
                        de_stats,  # Pass all results for visualization
                        cell_type,
                        conditions,
                        f'de_analysis_{cell_type}_{conditions[1]}_vs_{conditions[0]}.png'
                    )
                    
                    # Save results
                    output_file = os.path.join(
                        self.output_dir, 
                        f'de_results_{cell_type}_{conditions[1]}_vs_{conditions[0]}.csv'
                    )
                    de_stats.to_csv(output_file, index=False)
                    
                    self.logger.info(
                        f'Found {len(significant_de)} significantly differentially expressed genes '
                        f'out of {len(available_genes)} analyzed genes in {cell_type}'
                    )
                else:
                    self.logger.info(
                        f'No significantly differentially expressed genes found in {cell_type} '
                        f'among {len(available_genes)} analyzed genes'
                    )
            
            return de_results


        except Exception as e:
            self.logger.error(f'Error performing differential expression: {e}')
            return None
        
    def _plot_de_results(self, de_stats, cell_type, conditions, filename):
        """Create volcano plot for DE results.
        
        Args:
            de_stats (pd.DataFrame): DE analysis results
            cell_type (str): Cell type being analyzed
            conditions (List[str]): Conditions being compared [reference, comparison]
            filename (str): Output filename
        """
        try:
            # Calculate -log10(pvals_adj)
            de_stats['neg_log10_pval'] = -np.log10(de_stats['pvals_adj'])
            
            # Add significance categories
            de_stats['significance'] = 'Not Significant'
            sig_mask = (de_stats['pvals_adj'] < 0.05) & (abs(de_stats['logfoldchanges']) > 1)
            de_stats.loc[sig_mask, 'significance'] = 'Significant'
            
            # Create plot
            plt.figure(figsize=(10, 8))
            
            # Create scatter plot
            sns.scatterplot(
                data=de_stats,
                x='logfoldchanges',
                y='neg_log10_pval',
                hue='significance',
                palette={'Significant': 'red', 'Not Significant': 'grey'},
                alpha=0.6
            )
            
            # Add threshold lines
            plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
            plt.axvline(x=-1, color='gray', linestyle='--', alpha=0.5)
            plt.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
            
            # Label top genes
            top_genes = de_stats[sig_mask].nlargest(10, 'neg_log10_pval')
            for _, gene in top_genes.iterrows():
                plt.annotate(
                    gene['names'],
                    xy=(gene['logfoldchanges'], gene['neg_log10_pval']),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=8,
                    alpha=0.7
                )
            
            # Customize plot
            plt.title(f'Differential Expression: {cell_type}\n{conditions[1]} vs {conditions[0]}')
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(Adjusted p-value)')
            
            # Add legend
            plt.legend(title='Significance', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Adjust layout and save
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, filename), dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f'Volcano plot saved as {filename}')
            
        except Exception as e:
            self.logger.error(f'Error plotting DE results: {e}')
        
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
            Run the entire pipeline based on the steps provided by the user
        """
        self.load_data()
        if 'quality_control' in self.steps:
            self.quality_control()
        else:
            self.logger.info('Skipping quality control step')
        if len(self.adata_list) > 1:
            self.integrate_data()
        else:
            self.logger.info('Only one dataset found after QC; skipping data integration step')
            self.adata = self.adata_list[0]
        if 'preprocessing' in self.steps:
            self.preprocess_data()
        else:
            self.logger.info('Skipping preprocessing step')
        if 'batch_correction' in self.steps:
            self.batch_effect_correction_combat()
        else:
            self.logger.info('Skipping batch correction step')
            
        if 'annotation' in self.steps:
            self.annotate_cell_types()
        else:
            self.logger.info('Skipping cell type annotation step')     
        if 'pca' in self.steps:
            self.run_pca()
        else:
            self.logger.info('Skipping PCA step')
        if 'clustering' in self.steps:
            self.clustering()
        else:
            self.logger.info('Skipping clustering step')
            
        # self.add_gene_expression_to_obs(['AR', 'METTL3']) # TODO: 수정
        #TODO: add skip step 
        # self.run_tsne(color_by='Group_Biopsy') # TODO: 수정
        # self.run_tsne(color_by='AR_METTL3_combined_status') # TODO: 수정
        # self.run_tsne(color_by='cell_type')
        # self.run_tsne(color_by='batch')
        # self.run_tsne(color_by='total_counts')
          
        if 'visualization' in self.steps:
            self.visualize_clusters()
        else:
            self.logger.info('Skipping visualization step')
        
        if 'correlation' in self.steps:
            self.compute_correlation()
        else:
            self.logger.info('Skipping correlation computation step')
        # if not self.cell_types:
        #     self.visualize_correlated_genes()
        # else:
        #     self.visualize_correlated_genes_combined()
        if 'differential_expression' in self.steps and hasattr(self, 'conditions'):
            self.perform_differential_expression()
        else:
            self.logger.info('Skipping differential expression analysis step')
        self.save_results()
        self.logger.info('Pipeline execution completed successfully')

    def analyze_gene_correlations(self, target_gene='AR', cell_type='Epithelial', genes=None, gene_file=None, 
                            min_expr=0.1, corr_threshold=0.3, pval_threshold=0.05):
        """Analyze correlations between target gene and gene list in specific cell type.
        
        Args:
            target_gene (str): Target gene to correlate against
            cell_type (str): Specific cell type to analyze
            genes (List[str], optional): List of genes to analyze correlations with
            gene_file (str, optional): Path to file containing gene list
            min_expr (float): Minimum expression threshold (fraction of cells)
            corr_threshold (float): Minimum correlation coefficient threshold
            pval_threshold (float): Maximum p-value threshold for significance
            
        Returns:
            pd.DataFrame: Sorted correlation results
        """
        try:
            self.logger.info(f'Starting correlation analysis for {target_gene} in {cell_type}')
            
            # Validate target gene
            if target_gene not in self.adata.var_names:
                self.logger.error(f'Target gene {target_gene} not found in dataset')
                return None
            
            # Validate cell type
            if cell_type not in self.adata.obs['cell_type'].unique():
                self.logger.error(f'Cell type {cell_type} not found in dataset')
                return None
            
            # Get default gene list if none provided
            if genes is None and gene_file is None:
                genes = [
                'ALKBH8', 'BCDIN3D', 'BMT2', 'BUD23', 'CBLL1', 'CDK5RAP1', 'CDKAL1', 'CEBPZ', 'CMTR1', 'CMTR2', 
                'DIMT1', 'EMG1', 'FBL', 'FBLL1', 'FDXACB1', 'FMR1', 'FTSJ1', 'FTSJ3', 'HENMT1', 'HSD17B10', 'LARP7', 
                'LCMT2', 'MEPCE', 'METTL1', 'METTL14', 'METTL15', 'METTL16', 'METTL2A', 'METTL2B', 'METTL3', 'METTL4', 
                'METTL5', 'METTL6', 'METTL7A', 'METTL7B', 'METTL8', 'MRM1', 'MRM2', 'MRM3', 'MTERF4', 'NOP2', 'NSUN2', 
                'NSUN3', 'NSUN4', 'NSUN5', 'NSUN6', 'NSUN7', 'PCIF1', 'PRORP', 'RAMAC', 'RBM15', 'RBM15B', 'RNGTT', 'RNMT', 
                'RRNAD1', 'RSAD1', 'SPOUT1', 'TARBP1', 'TFB1M', 'TFB2M', 'TGS1', 'THADA', 'THUMPD2', 'THUMPD3', 'TRDMT1',
                'TRIT1', 'TRMO', 'TRMT1', 'TRMT10A', 'TRMT10B', 'TRMT10C', 'TRMT11', 'TRMT112', 'TRMT12', 'TRMT13',
                'TRMT1L', 'TRMT2A', 'TRMT2B', 'TRMT44', 'TRMT5', 'TRMT6', 'TRMT61A', 'TRMT61B', 'TRMT9B', 'TRMU', 'TYW3', 
                'VIRMA', 'WDR4', 'WDR6', 'WTAP', 'ZC3H13', 'ZCCHC4'
            ]
            
            # Load genes from file if provided
            if gene_file:
                file_genes = self._load_genes_from_file(gene_file)
                if file_genes is None:
                    return None
                genes = list(set(genes + file_genes)) if genes else file_genes
            
            # Remove target gene from analysis genes if present
            if target_gene in genes:
                genes.remove(target_gene)
            
            # Filter for available genes
            available_genes = [gene for gene in genes if gene in self.adata.var_names]
            missing_genes = set(genes) - set(available_genes)
            
            if missing_genes:
                self.logger.warning(f'Genes not found in dataset: {", ".join(missing_genes)}')
            
            if not available_genes:
                self.logger.error('No valid genes found for correlation analysis')
                return None
            
            # Subset data for cell type
            cell_mask = self.adata.obs['cell_type'] == cell_type
            adata_subset = self.adata[cell_mask]
            
            # Calculate expression frequency
            n_cells = adata_subset.n_obs
            expr_freq = (adata_subset.X > 0).sum(axis=0) / n_cells
            
            # Filter genes by expression threshold
            expr_mask = expr_freq > min_expr
            genes_to_analyze = [gene for gene, mask in zip(available_genes, expr_mask) if mask]
            
            if not genes_to_analyze:
                self.logger.error(f'No genes pass expression threshold of {min_expr}')
                return None
            
            # Get expression data
            target_expr = adata_subset[:, target_gene].X.toarray().flatten()
            
            # Initialize results dictionary
            results = []
            
            # Calculate correlations
            for gene in genes_to_analyze:
                gene_expr = adata_subset[:, gene].X.toarray().flatten()
                corr_coef, p_value = pearsonr(target_expr, gene_expr)
                
                results.append({
                    'gene': gene,
                    'correlation': corr_coef,
                    'p_value': p_value,
                    'expression_frequency': expr_freq[available_genes.index(gene)]
                })
                
            # Convert to DataFrame and sort
            results_df = pd.DataFrame(results)
            
            # Filter significant correlations
            significant_df = results_df[
                (abs(results_df['correlation']) >= corr_threshold) & 
                (results_df['p_value'] <= pval_threshold)
            ]
            
            if significant_df.empty:
                self.logger.info('No significant correlations found')
                return results_df.sort_values('correlation', ascending=False)
            
            # Sort by absolute correlation value
            significant_df = significant_df.sort_values('correlation', ascending=False)
            
            # Save results
            output_file = os.path.join(
                self.output_dir,
                f'correlation_results_{cell_type}_{target_gene}.csv'
            )
            significant_df.to_csv(output_file, index=False)
            
            # Create visualization
            self._plot_correlation_results(
                significant_df,
                target_gene,
                cell_type,
                f'correlation_plot_{cell_type}_{target_gene}.png'
            )
            
            self.logger.info(
                f'Found {len(significant_df)} significant correlations out of {len(genes_to_analyze)} analyzed genes'
            )
            
            return significant_df
            
        except Exception as e:
            self.logger.error(f'Error in correlation analysis: {e}')
            return None

    def _plot_correlation_results(self, results_df, target_gene, cell_type, filename):
        """Create visualization for correlation results.
        
        Args:
            results_df (pd.DataFrame): Correlation analysis results
            target_gene (str): Target gene name
            cell_type (str): Cell type analyzed
            filename (str): Output filename
        """
        try:
            plt.figure(figsize=(12, 6))
            
            # Create bar plot
            sns.barplot(
                data=results_df.head(20),  # Top 20 correlations
                x='gene',
                y='correlation',
                palette='coolwarm'
            )
            
            # Customize plot
            plt.xticks(rotation=45, ha='right')
            plt.title(f'Top Correlations with {target_gene} in {cell_type}')
            plt.xlabel('Gene')
            plt.ylabel('Correlation Coefficient')
            
            # Add significance markers
            for i, p_val in enumerate(results_df.head(20)['p_value']):
                if p_val <= 0.001:
                    plt.text(i, 0, '***', ha='center', va='bottom')
                elif p_val <= 0.01:
                    plt.text(i, 0, '**', ha='center', va='bottom')
                elif p_val <= 0.05:
                    plt.text(i, 0, '*', ha='center', va='bottom')
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, filename), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            self.logger.error(f'Error plotting correlation results: {e}')

#TODO: 삭제
def init_process(log_queue):
    logger = logging.getLogger('ScRNAseqPipeline')
    logger.handlers = []
    logger.setLevel(logging.INFO)
    queue_handler = QueueHandler(log_queue)
    logger.addHandler(queue_handler)  
    
def parse_arguments():
    parser = argparse.ArgumentParser(description="Single cell RNA-seq pipeline")
    parser.add_argument('--data_paths', type=str, nargs='+', required=True, help='Path to Cell Ranger output directories (filtered_feature_bc_matrix) or h5ad data formats')
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
    
    return parser.parse_args()
    
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
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
from sklearn.neighbors import NearestNeighbors
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor, as_completed
from joblib import Parallel, delayed
from tqdm import tqdm
 
# Define Min-Max Scaling Function
def min_max_scale(data, feature_range=(0, 1)):
    if issparse(data):
        data = data.toarray()
    min_val = np.min(data, axis=0)
    max_val = np.max(data, axis=0)
    diff = max_val - min_val  
    scale = np.where(diff != 0, (feature_range[1] - feature_range[0]) / diff, 0)
    scaled_data = scale * (data - min_val) + feature_range[0]
    scaled_data[:, diff == 0] = feature_range[0]
    if isinstance(data, csr_matrix):
        scaled_data = csr_matrix(scaled_data)
    return scaled_data

# Preprocessing Function

def preprocess_adata(adata, min_cells=20):
    # 确保基因名称唯一
    adata.var_names_make_unique()

    sc.pp.normalize_total(adata)
    # 对数转换
    sc.pp.log1p(adata)
    # 使用min-max缩放 (可以考虑 sc.pp.scale 进行Z-score标准化)
    adata.X = min_max_scale(adata.X)

    if isinstance(adata.X, np.ndarray):
        expressed_in_spots = np.sum(adata.X > 0, axis=0)
    else:
    # 对于稀疏矩阵，使用这种方法
        expressed_in_spots = np.sum(adata.X > 0, axis=0).A1  # 使用.A1将numpy矩阵转换为1维数组

# 筛选出至少在20个空间点中有表达的基因
    adata = adata[:, expressed_in_spots >=  min_cells]

    return adata

# Compute Neighbors Function
def compute_neighbors(adata, k_neighbors=6):
    spatial_coords = adata.obsm['spatial']
    nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1).fit(spatial_coords)
    distances, indices = nbrs.kneighbors(spatial_coords)
    neighbors = {i: indices[i][1:].tolist() for i in range(spatial_coords.shape[0])}
    return neighbors

# Merge Results and Compute Max R Values
def merge_files_and_compute_max(results_dict):
    all_results = []
    for exp_thre, results in results_dict.items():
        for result in results:
            result['exp_thre'] = exp_thre
            all_results.append(result)

    results_df = pd.DataFrame(all_results)
    results_df = results_df[['Gene', 'R_g']]
    
    merged_df = results_df.groupby('Gene').max('R_g').reset_index()
    merged_df.sort_values(by='R_g', ascending=False, inplace=True)
    merged_df.index = range(merged_df.shape[0])

    return merged_df


# ------------------ Compute Gene Groups ----------------------
def compute_gene_groups_for_gene(gene_index, gene_name, expression_matrix, exp_thre):
    gene_expression = expression_matrix[:, gene_index].flatten()
    quantile_thre = np.percentile(gene_expression, exp_thre)
    L_g = np.where(gene_expression <= quantile_thre)[0]
    H_g = np.where(gene_expression > quantile_thre)[0]
    return gene_name, L_g, H_g

def compute_gene_groups(adata, exp_thre, num_workers=1):
    expression_matrix = adata.X.toarray() if issparse(adata.X) else adata.X
    results = {}
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(compute_gene_groups_for_gene, i, gene_name, expression_matrix, exp_thre)
                   for i, gene_name in enumerate(adata.var_names)]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Computing gene groups"):
            gene_name, L_g, H_g = future.result()
            results[gene_name] = {'L_g': L_g, 'H_g': H_g}
    adata.var['L_g'] = adata.var.index.map(lambda x: results[x]['L_g'] if x in results else [])
    adata.var['H_g'] = adata.var.index.map(lambda x: results[x]['H_g'] if x in results else [])
    return results

# ------------------ Compute Weight Matrices ------------------
def compute_weight_for_gene(adata, gene_name, neighbors):
    H_g = set(adata.var.loc[gene_name, 'H_g'])
    L_g = set(adata.var.loc[gene_name, 'L_g'])
    W_gs = {}
    for s, n_s in neighbors.items():
        n_s_set = set(n_s)
        if s in H_g:
            intersection_size = len(n_s_set & H_g)
            W_gs[s] = intersection_size / len(n_s_set) if len(n_s_set) > 0 else 0
        elif s in L_g:
            intersection_size = len(n_s_set & L_g)
            W_gs[s] = intersection_size / len(n_s_set) if len(n_s_set) > 0 else 0
    return gene_name, W_gs

def compute_weight_matrices(adata, neighbors, num_workers):
    W_gs_all_genes = {}
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(compute_weight_for_gene, adata, gene_name, neighbors)
                   for gene_name in adata.var_names]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Computing weight matrices"):
            gene_name, W_gs = future.result()
            W_gs_all_genes[gene_name] = W_gs
    return W_gs_all_genes


#  ------------------ Compute R Values ------------------
def compute_r_value_for_gene(gene_name, weights, expression_df, adata):
    H_g = set(adata.var.loc[gene_name, 'H_g'])
    L_g = set(adata.var.loc[gene_name, 'L_g'])
    gene_exp_list = list(expression_df[gene_name])
    mean_expression = np.mean(gene_exp_list)
    var_expression = np.var(gene_exp_list)
    C_g = var_expression / mean_expression if mean_expression != 0 and var_expression != 0 else 1
    H_g_big = [s for s in H_g if weights[s] > 0]
    L_g_big = [s for s in L_g if weights[s] > 0]

    avg_H_g_big_L_g_zero = (sum( gene_exp_list[s] for s in H_g_big) / len(H_g_big)) if len(H_g_big)>0 else 0
    avg_L_g_big_H_g_zero = (sum(gene_exp_list[s] for s in L_g_big) / len(L_g_big)) if len(L_g_big)>0 else 0

    B_g = avg_H_g_big_L_g_zero - avg_L_g_big_H_g_zero
    
    A_g = sum(weights[s] for s in H_g) / len(H_g) if H_g else 0
    
    

    return gene_name, A_g, B_g, C_g

def compute_r_values(adata, W_gs_all_genes, num_workers):
    expression_df = adata.to_df()
    results = []
    R_g_values = {}
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(compute_r_value_for_gene, gene_name, W_gs, expression_df, adata)
                   for gene_name, W_gs in W_gs_all_genes.items()]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Computing R values"):
            gene_name, A_g, B_g, C_g = future.result()
            #R_g = A_g  
            R_g = A_g * B_g * C_g
            results.append({'Gene': gene_name, 'A_g': A_g, 'B_g': B_g, 'C_g': C_g, 'R_g': R_g})
            R_g_values[gene_name] = R_g

    print(f' R_g =A_g * B_g * C_g ')        
    return results, R_g_values


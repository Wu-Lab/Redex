import os
import ot
import torch
import random
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from torch.backends import cudnn
from scipy.sparse import issparse
from scipy.sparse.csc import csc_matrix
from scipy.sparse.csr import csr_matrix
from sklearn.neighbors import NearestNeighbors 
import pandas as pd





def construct_interaction(adata, n_neighbors=3):
    """Constructing spot-to-spot interactive gaph"""
    position = adata.obsm['spatial']
    
    # calculate distance matrix
    distance_matrix = ot.dist(position, position, metric='euclidean')
    n_spot = distance_matrix.shape[0]
    
    adata.obsm['distance_matrix'] = distance_matrix
    
    # find k-nearest neighbors
    interaction = np.zeros([n_spot, n_spot])  
    for i in range(n_spot):
        vec = distance_matrix[i, :]
        distance = vec.argsort()
        for t in range(1, n_neighbors + 1):
            y = distance[t]
            interaction[i, y] = 1
         
    adata.obsm['graph_neigh'] = interaction
    
    #transform adj to symmetrical adj
    adj = interaction
    adj = adj + adj.T
    adj = np.where(adj>1, 1, adj)
    
    adata.obsm['adj'] = adj
    #print('n_neighbors',n_neighbors)
    
def construct_interaction_KNN(adata, n_neighbors=3):
    position = adata.obsm['spatial']
    n_spot = position.shape[0]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1).fit(position)  
    _ , indices = nbrs.kneighbors(position)
    x = indices[:, 0].repeat(n_neighbors)
    y = indices[:, 1:].flatten()
    interaction = np.zeros([n_spot, n_spot])
    interaction[x, y] = 1
    
    adata.obsm['graph_neigh'] = interaction
    
    #transform adj to symmetrical adj
    adj = interaction
    adj = adj + adj.T
    adj = np.where(adj>1, 1, adj)
    
    adata.obsm['adj'] = adj
    print('Graph constructed!')   




def min_max_scale(data, feature_range=(0, 1)):
    if issparse(data):
        data = data.toarray()
    min_val = np.min(data, axis=0)
    max_val = np.max(data, axis=0)
    diff = max_val - min_val  
    scale = np.zeros_like(diff)
    scale[diff != 0] = (feature_range[1] - feature_range[0]) / diff[diff != 0]
    
    scaled_data = scale * (data - min_val) + feature_range[0]
    scaled_data[:, diff == 0] = feature_range[0]
    
    if isinstance(data, csr_matrix):
        scaled_data = csr_matrix(scaled_data)
        
    return scaled_data

# def preprocess(adata,dataset_name,gene_path):
#     sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)
#     sc.pp.scale(adata, zero_center=False, max_value=10)
#     #adata.X = min_max_scale(adata.X)
#     print('preprocess  gene use seuratv3-3000')

def preprocess(adata, dataset_name,gene_path):
    # 根据 dataset_name 动态构建 CSV 文件的路径
    #file_path = f'D:/Users/jtyu/新建文件夹/NEW MODEL WITH NEW GENE COMBINE/SVG/results/{dataset_name}/three_top_rg-results_df.csv'
    #file_path= os.path.join(gene_path, 'SPARX.csv')
    #file_path= os.path.join(gene_path, 'PI.csv')
   # file_path= os.path.join(gene_path, 'SPA.csv')
    file_path= os.path.join(gene_path, 'three_top_rg-results_df.csv')
    # 读取 CSV 文件中的基因列表
    gene_list_df = pd.read_csv(file_path)
    # 假设 CSV 文件中包含一个名为 'Gene' 的列，包含了基因名
    gene_list = gene_list_df['Gene'].values[:3000]  # 获取前3000个基因名
    #gene_list = gene_list_df['Gene'].values[:100]  # 获取前3000个基因名
   # gene_list = gene_list_df['gene_name'].values[:500]  # 获取前3000个基因名
    #gene_list = gene_list_df['gene'].values[:500]  # 获取前3000个基因名
    # 设置这些基因为高变异基因
    adata.var['highly_variable'] = adata.var_names.isin(gene_list)
    # 标准的预处理步骤
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False, max_value=10)
    adata.X = min_max_scale(adata.X)
    print(f'processed {dataset_name} our gene 3000')







def get_feature(adata, deconvolution=False):
    if deconvolution:
       adata_Vars = adata
    else:   
       adata_Vars =  adata[:, adata.var['highly_variable']]
       
    if isinstance(adata_Vars.X, csc_matrix) or isinstance(adata_Vars.X, csr_matrix):
       feat = adata_Vars.X.toarray()[:, ]
    else:
       feat = adata_Vars.X[:, ] 
    
    adata.obsm['feat'] = feat
      
    

    
def normalize_adj(adj):
    """Symmetrically normalize adjacency matrix."""
    adj = sp.coo_matrix(adj)
    rowsum = np.array(adj.sum(1))
    d_inv_sqrt = np.power(rowsum, -0.5).flatten()
    d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.
    d_mat_inv_sqrt = sp.diags(d_inv_sqrt)
    adj = adj.dot(d_mat_inv_sqrt).transpose().dot(d_mat_inv_sqrt)
    return adj.toarray()

def preprocess_adj(adj):
    """Preprocessing of adjacency matrix for simple GCN model and conversion to tuple representation."""
    adj_normalized = normalize_adj(adj)+np.eye(adj.shape[0])
    return adj_normalized 



def fix_seed(seed):
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    cudnn.deterministic = True
    cudnn.benchmark = False
    
    os.environ['PYTHONHASHSEED'] = str(seed)
    os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' 
    
    

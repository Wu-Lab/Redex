import os
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.metrics import normalized_mutual_info_score
from sklearn.decomposition import PCA
import scanpy as sc
import ot
import matplotlib.pyplot as plt


# ============================================================
#  Clustering
# ============================================================

def mclust_R(adata, num_cluster, modelNames='EEE', used_obsm='emb_pca', random_seed=2024):
    """Clustering using the mclust algorithm (requires R and rpy2)."""
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[used_obsm]), num_cluster, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs['mclust'] = mclust_res
    adata.obs['mclust'] = adata.obs['mclust'].astype('int').astype('category')
    return adata


def clustering(adata, n_clusters=7, radius=50, key='emb', method='mclust',
               start=0.1, end=2.0, increment=0.05, refinement=False):
    """
    Spatial clustering based on the learned representation.

    Parameters
    ----------
    adata : AnnData
    n_clusters : int
    radius : int
        Number of neighbors for spatial refinement.
    method : str
        'mclust', 'leiden', or 'louvain'.
    refinement : bool
        Whether to apply spatial refinement.
    """
    pca = PCA(n_components=20, random_state=2024)
    #pca = PCA(n_components=10, random_state=2024)
    embedding = pca.fit_transform(adata.obsm['reconsfeat'].copy())
    adata.obsm['emb_pca'] = embedding

    if method == 'mclust':
        adata = mclust_R(adata, used_obsm='emb_pca', num_cluster=n_clusters)
        adata.obs['domain'] = adata.obs['mclust']
    elif method == 'leiden':
        res = search_res(adata, n_clusters, use_rep='emb_pca', method=method,
                         start=start, end=end, increment=increment)
        sc.tl.leiden(adata, random_state=0, resolution=res)
        adata.obs['domain'] = adata.obs['leiden']
    elif method == 'louvain':
        res = search_res(adata, n_clusters, use_rep='emb_pca', method=method,
                         start=start, end=end, increment=increment)
        sc.tl.louvain(adata, random_state=0, resolution=res)
        adata.obs['domain'] = adata.obs['louvain']

    if refinement:
        new_type = refine_label(adata, radius, key='domain')
        adata.obs['domain'] = new_type


def refine_label(adata, radius=50, key='label'):
    """Refine cluster labels using spatial neighbors."""
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values
    position = adata.obsm['spatial']
    distance = ot.dist(position, position, metric='euclidean')
    n_cell = distance.shape[0]

    for i in range(n_cell):
        vec = distance[i, :]
        index = vec.argsort()
        neigh_type = []
        for j in range(1, n_neigh + 1):
            neigh_type.append(old_type[index[j]])
        max_type = max(neigh_type, key=neigh_type.count)
        new_type.append(max_type)

    new_type = [str(i) for i in list(new_type)]
    return new_type


def search_res(adata, n_clusters, method='leiden', use_rep='emb',
               start=0.1, end=3.0, increment=0.01):
    """Search for the resolution that yields the target number of clusters."""
    print('Searching resolution...')
    label = 0
    sc.pp.neighbors(adata, n_neighbors=30, use_rep=use_rep)
    for res in sorted(list(np.arange(start, end, increment)), reverse=True):
        if method == 'leiden':
            sc.tl.leiden(adata, random_state=0, resolution=res)
            count_unique = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        elif method == 'louvain':
            sc.tl.louvain(adata, random_state=0, resolution=res)
            count_unique = len(pd.DataFrame(adata.obs['louvain']).louvain.unique())
        print(f'resolution={res}, cluster number={count_unique}')
        if count_unique == n_clusters:
            label = 1
            break
    assert label == 1, "Resolution not found. Please try a bigger range or smaller step."
    return res


# ============================================================
#  Evaluation
# ============================================================

def evaluate(adata, pred_key='domain', truth_key='ground_truth'):
    """
    Compute ARI and NMI between predicted and ground-truth labels.
    Automatically filters out spots without ground-truth annotation.

    Parameters
    ----------
    adata : AnnData
    pred_key : str
        Key in adata.obs for predicted cluster labels.
    truth_key : str
        Key in adata.obs for ground-truth labels.

    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData (spots with ground truth only).
    ari : float
    nmi : float
    """
    adata_filtered = adata[~pd.isnull(adata.obs[truth_key])].copy()
    ari = metrics.adjusted_rand_score(adata_filtered.obs[pred_key], adata_filtered.obs[truth_key])
    nmi = normalized_mutual_info_score(adata_filtered.obs[pred_key], adata_filtered.obs[truth_key])
    adata_filtered.uns['ARI'] = ari
    adata_filtered.uns['NMI'] = nmi
    print(f'ARI = {ari:.4f}')
    print(f'NMI = {nmi:.4f}')
    return adata_filtered, ari, nmi


# ============================================================
#  Visualization
# ============================================================

def plot_spatial(adata, color_key, title=None, figsize=(3, 3), spot_size=100, save_path=None):
    """Plot spatial scatter colored by a given obs key."""
    print(f"[Spatial plot] figsize={figsize}, spot_size={spot_size}")
    plt.rcParams["figure.figsize"] = figsize
    sc.pl.spatial(adata, color=color_key, title=title,
                  spot_size=spot_size, show=False)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    plt.show()
    plt.close()


def plot_umap(adata, color_key, title=None, save_path=None, use_rep='emb_pca'):
    """Compute (if needed) and plot UMAP colored by a given obs key."""
    if 'X_umap' not in adata.obsm:
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.umap(adata)

    plt.figure(figsize=(3, 3))
    sc.pl.umap(adata, color=color_key, title=title,
               frameon=False, size=50, show=False)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    plt.show()
    plt.close()


def plot_paga(adata, color_key='domain', title='PAGA', save_path=None):
    """Compute and plot PAGA graph."""
    sc.tl.paga(adata, groups=color_key)
    plt.figure(figsize=(3, 3))
    sc.pl.paga(adata, color=color_key, title=title,
               frameon=False, show=False,
               node_size_scale=6, edge_width_scale=2,
               fontsize=12, fontoutline=1)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    plt.show()
    plt.close()

def find_marker_genes(adata, output_path, cluster_key='mclust', quant_in=30, quant_out=70, eps=1e-9):
    """
    针对每个 domain，计算每个基因的 marker 指标并保存 CSV。
    指标 = domain 内 quant_in 分位值 / (domain 外 quant_out 分位值 + eps)
    """
    os.makedirs(output_path, exist_ok=True)
    expr_df = adata.to_df()
    domains = adata.obs[cluster_key].unique()
    genes = adata.var_names
    marker_results = {}

    for domain in domains:
        in_mask = adata.obs[cluster_key] == domain
        expr_in = expr_df.loc[in_mask]
        expr_out = expr_df.loc[~in_mask]

        ratios = {}
        for gene in genes:
            q_in = np.percentile(expr_in[gene].values, quant_in)
            q_out = np.percentile(expr_out[gene].values, quant_out) if len(expr_out) > 0 else np.nan
            ratios[gene] = q_in / (q_out + eps)

        marker_df = pd.DataFrame(list(ratios.items()), columns=['gene', 'marker_score'])
        marker_df.sort_values('marker_score', ascending=False, inplace=True)
        marker_df.to_csv(os.path.join(output_path, f"domain_{int(domain):02d}_marker_genes.csv"), index=False)
        marker_results[domain] = marker_df

    print(f"Marker genes saved to {output_path}")
    return marker_results
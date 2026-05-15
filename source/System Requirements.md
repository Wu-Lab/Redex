# System Requirements

`Redex` can run on Linux, macOS, and Windows when the required scientific Python packages are available. GPU acceleration is recommended for graph neural network training, but CPU execution is supported for smaller datasets.

Recommended environment:

- Python 3.8 or 3.9
- CUDA-compatible PyTorch installation for GPU usage
- R and the R package `mclust` if using `method="mclust"` in the clustering helper

Core Python dependencies are listed in `requirements.txt`:

```text
numpy==1.22.3
pandas==1.4.2
scipy==1.8.1
scikit-learn==1.1.1
scikit-misc==0.1.4
matplotlib==3.7.5
scanpy==1.9.8
anndata==0.9.2
tqdm==4.64.0
joblib==1.4.2
POT==0.9.4
igraph==0.11.5
rpy2==3.4.1
```

Optional dependencies:

- `rpy2` for calling R `mclust` from Python
- R package `mclust` for Gaussian mixture clustering
- Jupyter or JupyterLab for running the tutorial notebooks

If `mclust` is not available, Redex utilities also support Leiden and Louvain clustering.

# Data preparation

Redex expects an `AnnData` object with:

- feature values in `adata.X`
- spatial coordinates in `adata.obsm["spatial"]`
- unique feature names in `adata.var_names`

All tutorial data can be downloaded from Quark Cloud:

[https://pan.quark.cn/s/2f695d3f9e90](https://pan.quark.cn/s/2f695d3f9e90)

After downloading, place the datasets under the repository root following the layouts used in the tutorial notebooks.

For 10x Visium datasets, a typical loader is:

```python
import scanpy as sc

adata = sc.read_visium("./DLPFC/151509")
adata.var_names_make_unique()
```

For h5ad datasets:

```python
import scanpy as sc

adata = sc.read_h5ad("./data/example.h5ad")
adata.var_names_make_unique()
```

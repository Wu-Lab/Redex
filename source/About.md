# About

`Redex` is a Python package for spatial omics analysis. It identifies spatially variable features based on regional differential expression and uses those features for spatial domain segmentation and marker feature discovery.

Redex contains three connected modules:

1. SVF selection ranks genes, peaks, or derived features by combining spatial coherence, regional expression contrast, and global expression heterogeneity.
2. Spatial clustering learns spot embeddings with a graph attention network in a self-supervised contrastive learning framework.
3. Marker feature selection ranks region-specific features by comparing within-domain and out-of-domain expression.

Redex can be applied to spatial transcriptomics, spatial ATAC-seq, and spatial gene-pair co-expression feature matrices.

```{image} images/Fig1.jpg
:alt: Redex workflow overview
:width: 95%
```

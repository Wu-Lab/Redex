# About

`Redex` is a Python package for spatial omics analysis. It identifies spatially variable features based on regional differential expression and uses those features for spatial domain segmentation and marker feature discovery.

Redex supports a workflow for spatial omics analysis:

1. **Spatially variable feature prioritization**  
   Redex ranks genes, peaks, or derived spatial features using the Redex score.

2. **Spatial domain segmentation**  
   Top-ranked features are used to learn spatially informed embeddings for identifying spatial domains.

3. **Marker feature detection**  
   Region-specific markers are detected by comparing within-domain and out-of-domain expression patterns.

4. **Downstream spatial interpretation**  
   Redex supports spatially variable gene detection, spatially variable peak identification, gene-gene relationship discovery, trajectory analysis, spatial domain segmentation, and marker feature identification.

Redex can be applied to spatial transcriptomics, spatial ATAC-seq, and derived spatial feature matrices such as local gene-pair co-expression profiles.
```{image} images/Fig1.jpg
:alt: Redex workflow overview
:width: 95%
```

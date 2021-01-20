

# Dimensional Reduction, Clustering, and Differential Gene expression testing

1) In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above.

---



<h3>  How to store and interact with dimensional reduction information (such as the output from RunPCA)  </h3>

https://satijalab.org/seurat/v3.0/dim_reduction_vignette.html

is this pre or post PCA?

load data

```R
library(Seurat)
pbmc <- readRDS(file = "../data/pbmc3k_final.rds")
```

<h4> cell.embeddings </h4> stores the coordinates for each cell in low-dimensional space.

<h4> feature.loadings</h4> stores the weight for each feature along each dimension of the embedding

<h4>feature.loadings.projected</h4> Seurat typically calculate the dimensional reduction on a subset of genes (for example, high-variance genes), and then project that structure onto the entire dataset (all genes). The results of that projection (calculated with ProjectDim) are stored in this slot. Note that the cell loadings will remain unchanged after projection but there are now feature loadings for all feature

<h4>stdev</h4> The standard deviations of each dimension. Most often used with PCA (storing the square roots of the eigenvalues of the covariance matrix) and can be useful when looking at the drop off in the amount of variance that is explained by each successive dimension.

<h4>key</h4> Sets the column names for the cell.embeddings and feature.loadings matrices. For example, for PCA, the column names are PC1, PC2, etc., so the key is “PC”.

<h4>jackstraw</h4> Stores the results of the jackstraw procedure run using this dimensional reduction technique. Currently supported only for PCA.

<h4>misc</h4> Bonus slot to store any other information you might want


---


2) The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component.

3) Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default.

4) To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space using the top 15 PCs.

5) In order to cluster the cells based on similarity of expression, we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.

6) For the myeloid, vascular, and macroglial cells, we performed similar analyses as described above, with a few modifications. In order to identify reproducible sub-clusters of cells, we performed the same graph-based clustering through a range of PCs, “k.neighbor” and “resolution” parameters and inspected cluster memberships for stable configurations. For the myeloid, vascular, and macroglia, we took the top 12, 11, and 8 PCs and set resolutions to 0.5, 0.3, and 0.45 respectively.

7) To identify marker genes for each cluster, we used the FindAllMarkers() function using default parameters, which implements a Wilcoxon Rank Sum test comparing gene expression of cells within a given cluster versus all other cells.

8) We repeated this analysis to identify marker genes distinguishing subsets within a cell-type.

9) To infer the functional relevance of sub-clusters, we performed gene ontology enrichment analyses on the top 50 differentially expressed genes using Fisher’s Exact test as implemented in the topGO R package. For the enrichment analyses of the gene expression changes in astrocytes, our initial analysis revealed very few differentially expressed genes between the uninjured and 1dpi astrocytes, which we attributed to the low numbers of uninjured astrocytes captured. Therefore, we supplemented our uninjured astrocyte dataset with ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset50. We also supplemented our uninjured OPC dataset in order to validate that our uninjured cells were more transcriptional similar to the external reference cells than to our injured cells.

10) To account for differences in sequencing depth between our dataset and the external dataset, we performed differential expression tests using MAST as implemented in Seurat. We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis.

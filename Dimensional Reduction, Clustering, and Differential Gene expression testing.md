

# Dimensional Reduction, Clustering, and Differential Gene expression testing

---

1) In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above.

---



<h3>  How to store and interact with dimensional reduction information (such as the output from RunPCA)  </h3>

"We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results)"

"By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA."

"Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA." (this step already done?)

https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

https://satijalab.org/seurat/v3.0/dim_reduction_vignette.html



---


2) The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component.


---

he first thing to look at is the PCA scree-plot, showing the proportion of variance explained by each component. We are looking for a “knee” in the plot, where additional PCs do not bring much more new information.

For this purpose, Seurat provides the function PCElbowPlot, that displays the standard-deviation of each PC.

`PCElbowPlot(sobj, num.pc = 40)`


https://gtpb.github.io/ADER18S/pages/tutorial-seurat-mca

--- 

3) Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default.

4) To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space using the top 15 PCs.

5) In order to cluster the cells based on similarity of expression, we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.

6) For the myeloid, vascular, and macroglial cells, we performed similar analyses as described above, with a few modifications. In order to identify reproducible sub-clusters of cells, we performed the same graph-based clustering through a range of PCs, “k.neighbor” and “resolution” parameters and inspected cluster memberships for stable configurations. For the myeloid, vascular, and macroglia, we took the top 12, 11, and 8 PCs and set resolutions to 0.5, 0.3, and 0.45 respectively.

7) To identify marker genes for each cluster, we used the FindAllMarkers() function using default parameters, which implements a Wilcoxon Rank Sum test comparing gene expression of cells within a given cluster versus all other cells.

8) We repeated this analysis to identify marker genes distinguishing subsets within a cell-type.

9) To infer the functional relevance of sub-clusters, we performed gene ontology enrichment analyses on the top 50 differentially expressed genes using Fisher’s Exact test as implemented in the topGO R package. For the enrichment analyses of the gene expression changes in astrocytes, our initial analysis revealed very few differentially expressed genes between the uninjured and 1dpi astrocytes, which we attributed to the low numbers of uninjured astrocytes captured. Therefore, we supplemented our uninjured astrocyte dataset with ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset50. We also supplemented our uninjured OPC dataset in order to validate that our uninjured cells were more transcriptional similar to the external reference cells than to our injured cells.

10) To account for differences in sequencing depth between our dataset and the external dataset, we performed differential expression tests using MAST as implemented in Seurat. We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis.

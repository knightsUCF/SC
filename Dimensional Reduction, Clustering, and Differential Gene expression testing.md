

# Dimensional Reduction, Clustering, and Differential Gene expression testing

---

<h3> 1) In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above. </h3>

---



<h3>  How to store and interact with dimensional reduction information (such as the output from RunPCA)  </h3>

"We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results)"

"By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA."

"Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA."

https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

https://satijalab.org/seurat/v3.0/dim_reduction_vignette.html



---


2) The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component.


---

The first thing to look at is the PCA scree-plot, showing the proportion of variance explained by each component. We are looking for a “knee” in the plot, where additional PCs do not bring much more new information.

For this purpose, Seurat provides the function PCElbowPlot, that displays the standard-deviation of each PC.

`PCElbowPlot(sobj, num.pc = 40)`


https://gtpb.github.io/ADER18S/pages/tutorial-seurat-mca

--- 

3) Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default.

---

Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We first determine the k-nearest neighbors of each cell. We use this knn graph to construct the SNN graph by calculating the neighborhood overlap (Jaccard index) between every cell and its k.param nearest neighbors.

https://rdrr.io/github/satijalab/seurat/man/FindNeighbors.html


---

4) To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space using the top 15 PCs.

---

```R
# Plot UMAP, coloring cells by cell type (currently stored in object@ident)
DimPlot(pbmc, reduction = "umap")
```

https://satijalab.org/seurat/v3.0/interaction_vignette.html

---

5) In order to cluster the cells based on similarity of expression, we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.

---

Cluster Determination

Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B.


```R
FindClusters(object, genes.use = NULL, pc.use = NULL, k.param = 30,
  k.scale = 25, plot.SNN = FALSE, prune.SNN = 1/15, save.SNN = FALSE,
  reuse.SNN = FALSE, do.sparse = FALSE, modularity.fxn = 1,
  resolution = 0.8, algorithm = 1, n.start = 100, n.iter = 10,
  random.seed = 0, print.output = TRUE)
```

https://www.rdocumentation.org/packages/Seurat/versions/1.4.0/topics/FindClusters

---

6) For the myeloid, vascular, and macroglial cells, we performed similar analyses as described above, with a few modifications. In order to identify reproducible sub-clusters of cells, we performed the same graph-based clustering through a range of PCs, “k.neighbor” and “resolution” parameters and inspected cluster memberships for stable configurations. For the myeloid, vascular, and macroglia, we took the top 12, 11, and 8 PCs and set resolutions to 0.5, 0.3, and 0.45 respectively.

---

from example above:

```resolution = 0.8, algorithm = 1, n.start = 100, n.iter = 10,```

---

7) To identify marker genes for each cluster, we used the FindAllMarkers() function using default parameters, which implements a Wilcoxon Rank Sum test comparing gene expression of cells within a given cluster versus all other cells.

---


```R
FindAllMarkers(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 0.01,
  ...
)
```


https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindAllMarkers

---

8) We repeated this analysis to identify marker genes distinguishing subsets within a cell-type.

---

```R
min.cells.feature = 3,
min.cells.group = 3,
```

---

9) To infer the functional relevance of sub-clusters, we performed gene ontology enrichment analyses on the top 50 differentially expressed genes using Fisher’s Exact test as implemented in the topGO R package. For the enrichment analyses of the gene expression changes in astrocytes, our initial analysis revealed very few differentially expressed genes between the uninjured and 1dpi astrocytes, which we attributed to the low numbers of uninjured astrocytes captured. Therefore, we supplemented our uninjured astrocyte dataset with ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset. We also supplemented our uninjured OPC dataset in order to validate that our uninjured cells were more transcriptional similar to the external reference cells than to our injured cells.

---

 "Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster"
 
 (Bioconductor package)
 
 ```R
 # install org.Mm.eg.db if not already installed (for mouse only)
# install.packages("org.Mm.eg.db") 
cluster0 <- SubsetData(experiment.aggregate, ident.use = 0)
expr <- cluster0@data
# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.75)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
    GOdata <- new("topGOdata",
        ontology = "BP", # use biological process ontology
        allGenes = geneList,
        geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
```
 

---

10) To account for differences in sequencing depth between our dataset and the external dataset, we performed differential expression tests using MAST as implemented in Seurat. We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis.

---

```
# Test for DE features using the MAST package
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "MAST"))
```

https://satijalab.org/seurat/v3.0/de_vignette.html



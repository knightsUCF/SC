# Analysis

<br>

<h2> Overview </h2>


Our goal is to replicate the results of the single cell paper using the Seurat library, and R.

<br>

<i>"Single cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord"</i>

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1.full

<i>Seurat</i>

https://satijalab.org/seurat/

<br>

<h2> Outline </h2>

We are using this outline to replicate the results in the study: 

https://github.com/knightsUCF/SC/blob/main/Outline.md

We are starting on the section: "Dimensional Reduction, Clustering, and Differential Gene expression testing"

TODO: genomic data sequencing might require starting on the data preprocessing section. Currently we are using the toy data given by the Seurat tutorial.

<br>

<h2> Libraries </h2>

We will need to import the following libraries:

```R
library(dplyr)
library(Seurat)
library(patchwork)
```

<br>

<h2> Processing Data </h2>

<i> "1) In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above." </i>

```R
variable_genes_count = 10

pbmc.data <- Read10X(data.dir = 'data/filtered_gene_bc_matrices/hg19')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), variable_genes_count)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

<br>

<h2> Selecting Principal Components </h2>

<i> "2) The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component." </i>

```R
ElbowPlot(pbmc, 15) # only returns a graph - https://github.com/satijalab/seurat/blob/b56d194939379460db23380426d3896b54d91ab6/R/visualization.R

# TODO: select only these "top 15 principal components"
```

<br> 

<h2> Finding Neighbors </h2>

<i> "3) Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default." </i>

```R
pbmc <- FindNeighbors(
  pbmc,
  reduction = "pca",
  dims = 1:10,
  assay = NULL,
  features = NULL,
  k.param = 20,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "rann",
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  do.plot = FALSE,
  graph.name = NULL,
)
```




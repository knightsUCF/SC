library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)



# 1) In order to measure the greatest gene expression variation among all the SCI cells,
# we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above.

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



# 2) The top 15 principal components (PCs) were selected based on the “elbow”
# point heuristic in a scree plot which quantifies the contribution of variance
# by each principal component.


PCElbowPlot(pbmc, num.pc = 15)


library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)





load_seurat_object <- function() {
        variable_genes_count = 10 # replace with 2000 later when custom data is processed

        pbmc.data <- Read10X(data.dir = '../Data/Sample/filtered_gene_bc_matrices/hg19')
        pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
        pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^M                  T-")
        pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
        pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
        pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
        top10 <- head(VariableFeatures(pbmc), variable_genes_count)
        all.genes <- rownames(pbmc)
        pbmc <- ScaleData(pbmc, features = all.genes)
        pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
        return (pbmc)
}



analyze <- function(seurat_object) {

        seurat_object <- FindNeighbors(
          seurat_object,
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
          do.plot = FALSE,  graph.name = NULL,
        )

        seurat_object <- FindClusters(seurat_object, resolution = 0.5)

        markers <- FindAllMarkers(
          seurat_object,
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
        )
}


analyze(load_seurat_object())


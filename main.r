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

ElbowPlot(pbmc, 15) # only returns a graph - https://github.com/satijalab/seurat/blob/b56d194939379460db23380426d3896b54d91ab6/R/visualization.R

# TODO: select only these "top 15 principal components"



# 3) Using these components, a nearest-neighbor graph and shared-nearest-neighbor
# graph were generated with “k.neighbors” set to 20 by default.

FindNeighbors(
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



# 4) To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space
# using the top 15 PCs.

# DimPlot(pbmc, reduction = "umap") # Plot UMAP, coloring cells by cell type, reticulate::py_install(packages ='umap-learn'), https://github.com/satijalab/seurat/issues/631


# 5) In order to cluster the cells based on similarity of expression,
# we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.


# FindClusters(pbmc, genes.use = NULL, pc.use = NULL, k.param = 30,
#             k.scale = 25, plot.SNN = FALSE, prune.SNN = 1/15, save.SNN = FALSE,
#             reuse.SNN = FALSE, do.sparse = FALSE, modularity.fxn = 1,
#             resolution = 0.8, algorithm = 1, n.start = 100, n.iter = 10,
#             random.seed = 0, verbose = T)

pbmc <- FindNeighbors(pbmc, dims = 1:10) # need this to run FindClusters() TODO: pick correct parameters for FindNeighbors() above
pbmc <- FindClusters(pbmc, resolution = 0.5)


# 6) For the myeloid, vascular, and macroglial cells,
# we performed similar analyses as described above, with a few modifications.
# In order to identify reproducible sub-clusters of cells,
# we performed the same graph-based clustering through a range of PCs,
# “k.neighbor” and “resolution” parameters and inspected cluster memberships for stable configurations.
# For the myeloid, vascular, and macroglia, we took the top 12, 11,
# and 8 PCs and set resolutions to 0.5, 0.3, and 0.45 respectively.

# TODO: repeat of above steps



# 7) To identify marker genes for each cluster, we used the FindAllMarkers() function
# using default parameters, which implements a Wilcoxon Rank Sum test comparing
# gene expression of cells within a given cluster versus all other cells.

FindAllMarkers(
  pbmc,
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



# 8) We repeated this analysis to identify marker genes distinguishing subsets within a cell-type.

# TODO: research


# 9) To infer the functional relevance of sub-clusters,
# we performed gene ontology enrichment analyses on the top 50 differentially
# expressed genes using Fisher’s Exact test as implemented in the topGO R package.
# For the enrichment analyses of the gene expression changes in astrocytes,
# our initial analysis revealed very few differentially expressed genes
# between the uninjured and 1dpi astrocytes, which we attributed to the
# low numbers of uninjured astrocytes captured.
# Therefore, we supplemented our uninjured astrocyte dataset with 
# ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset.
# We also supplemented our uninjured OPC dataset in order to validate
# that our uninjured cells were more transcriptional similar to the
# external reference cells than to our injured cells.


# TODO: research


# 10) To account for differences in sequencing depth between our dataset and the external dataset,
# we performed differential expression tests using MAST as implemented in Seurat.
# We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis.


# TODO: research


# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Calculating gene variances
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Calculating feature variances of standardized and clipped values
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            Centering and scaling data matrix
#          |====================================================================================================================================| 100%
#          PC_ 1 
#          Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
#          FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
#          PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
#          Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
#          CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
#          MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
#          PC_ 2 
#          Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
#          HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
#          BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
#          Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
#          CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
#          TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
#          PC_ 3 
#          Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
#          HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
#          PLAC8, BLNK, MALAT1, SMIM14, PLD4, LAT2, IGLL5, P2RX5, SWAP70, FCGR2B 
#          Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
#          HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
#          NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
#          PC_ 4 
#          Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HLA-DPB1, HIST1H2AC, PF4, TCL1A 
#          SDPR, HLA-DPA1, HLA-DRB1, HLA-DQA2, HLA-DRA, PPBP, LINC00926, GNG11, HLA-DRB5, SPARC 
#          GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, GZMB, CLU, TUBB1 
#          Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
#          AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
#          LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
#          PC_ 5 
#          Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
#          GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, CCL5 
#          RBP7, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
#          Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
#          PTGES3, LILRB2, MAL, CD27, HN1, CD2, GDI2, ANXA5, CORO1B, TUBA1B 
#          FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB 
#          Computing nearest neighbor graph
#          Computing SNN
#          Computing nearest neighbor graph
#          Computing SNN
#          Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#          
#          Number of nodes: 2638
#          Number of edges: 95965
#          
#          Running Louvain algorithm...
#          0%   10   20   30   40   50   60   70   80   90   100%
#            [----|----|----|----|----|----|----|----|----|----|
#               **************************************************|
#               Maximum modularity in 10 random starts: 0.8723
#             Number of communities: 9
#             Elapsed time: 0 seconds
#             Calculating cluster 0
#             For a more efficient implementation of the Wilcoxon Rank Sum Test,
#             (default method for FindMarkers) please install the limma package
#             --------------------------------------------
#               install.packages('BiocManager')
#             BiocManager::install('limma')
#             --------------------------------------------
#               After installation of limma, Seurat will automatically use the more 
#             efficient implementation (no further action necessary).
#             This message will be shown once per session
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=38s  
#             Calculating cluster 1
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=56s  
#             Calculating cluster 2
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=27s  
#             Calculating cluster 3
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=42s  
#             Calculating cluster 4
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=35s  
#             Calculating cluster 5
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=52s  
#             Calculating cluster 6
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=58s  
#             Calculating cluster 7
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=56s  
#             Calculating cluster 8
#             |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 29s

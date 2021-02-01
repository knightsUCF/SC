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

TODO: raw genomic data might require starting on the data preprocessing section. Currently we are using the toy data given by the Seurat tutorial.

<br>

<h2> Libraries </h2>

We will need to import the following libraries:

```R
library(dplyr)
library(Seurat)
library(patchwork)
library(topGO)
```

<br>

<h2> 1. Processing Data </h2>

<i> "In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above." </i>

```R
variable_genes_count = 10 # replace with 2000 later when custom data is processed

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

<h4>Output</h4>


	Performing log-normalization
	0%   10   20   30   40   50   60   70   80   90   100%
	[----|----|----|----|----|----|----|----|----|----|
	**************************************************|
	Calculating gene variances
	0%   10   20   30   40   50   60   70   80   90   100%
	[----|----|----|----|----|----|----|----|----|----|
	**************************************************|
	Calculating feature variances of standardized and clipped values
	0%   10   20   30   40   50   60   70   80   90   100%
	[----|----|----|----|----|----|----|----|----|----|
	**************************************************|
	Centering and scaling data matrix
	  |==============================================================================================================================================================================================| 100%
	PC_ 1 
	Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
		   FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
		   PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
	Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
		   CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
		   MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
	PC_ 2 
	Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
		   HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
		   BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
	Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
		   CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
		   TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
	PC_ 3 
	Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
		   HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
		   PLAC8, BLNK, MALAT1, SMIM14, PLD4, LAT2, IGLL5, P2RX5, SWAP70, FCGR2B 
	Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
		   HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
		   NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
	PC_ 4 
	Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HLA-DPB1, HIST1H2AC, PF4, TCL1A 
		   SDPR, HLA-DPA1, HLA-DRB1, HLA-DQA2, HLA-DRA, PPBP, LINC00926, GNG11, HLA-DRB5, SPARC 
		   GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, GZMB, CLU, TUBB1 
	Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
		   AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
		   LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
	PC_ 5 
	Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
		   GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, CCL5 
		   RBP7, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
	Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
		   PTGES3, LILRB2, MAL, CD27, HN1, CD2, GDI2, ANXA5, CORO1B, TUBA1B 
		   FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB 
     


<br>

<h2> 2. Selecting Principal Components </h2>

<i> "The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component." </i>

```R
plot(ElbowPlot(pbmc, 15)) # only returns a graph - https://github.com/satijalab/seurat/blob/b56d194939379460db23380426d3896b54d91ab6/R/visualization.R

# TODO: select only these "top 15 principal components"
```

<h3>Output</h3>

![Elbow Plot](https://github.com/knightsUCF/SC/blob/main/charts/Rplot01.png)

<br> 

<h2> 3. Finding Neighbors </h2>

<i> "Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default." </i>

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

<h4>Output</h4>

	Computing nearest neighbor graph
	Computing SNN


<br>

<h2>4. UMAP Plot </h2>

<i>To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space using the top 15 PCs. </i>

```R
# DimPlot(pbmc, reduction = "umap") # Plot UMAP, coloring cells by cell type, reticulate::py_install(packages ='umap-learn'), https://github.com/satijalab/seurat/issues/631
```

TODO: resolve dependencies

<br>

<h2>5. Finding Clusters</h2>

<i>In order to cluster the cells based on similarity of expression, we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.</i>

TODO: verify that we are running the method on the correct "shared-nearest-neighbor graph"

```R
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

<h4>Output</h4>

```
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 2638
Number of edges: 96033

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8720
Number of communities: 9
Elapsed time: 0 seconds
```

<br>

<h2>6. Repeating Process on Different Set</h2>

<i>For the myeloid, vascular, and macroglial cells, we performed similar analyses as described above, with a few modifications. In order to identify reproducible sub-clusters of cells, we performed the same graph-based clustering through a range of PCs, “k.neighbor” and “resolution” parameters and inspected cluster memberships for stable configurations. For the myeloid, vascular, and macroglia, we took the top 12, 11, and 8 PCs and set resolutions to 0.5, 0.3, and 0.45 respectively.</i>

TODO: waiting on 10X genomic data processing pipeline on custom data, and which set of cells this should be run on

<br>

<h2>7. Finding Marker Genes</h2>

<i>To identify marker genes for each cluster, we used the FindAllMarkers() function, using default parameters, which implements a Wilcoxon Rank Sum test comparing, gene expression of cells within a given cluster versus all other cells.</i>

```R
markers <- FindAllMarkers(
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
```

<h4>Output</h4>

```
Calculating cluster 0
For a more efficient implementation of the Wilcoxon Rank Sum Test,
(default method for FindMarkers) please install the limma package
--------------------------------------------
install.packages('BiocManager')
BiocManager::install('limma')
--------------------------------------------
After installation of limma, Seurat will automatically use the more 
efficient implementation (no further action necessary).
This message will be shown once per session
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=36s  
Calculating cluster 1
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=23s  
Calculating cluster 2
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=53s  
Calculating cluster 3
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=35s  
Calculating cluster 4
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=27s  
Calculating cluster 5
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=51s  
Calculating cluster 6
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=53s  
Calculating cluster 7
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=51s  
Calculating cluster 8
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 28s
```

<h4>Output of print(markers) </h4>

```
                  p_val avg_log2FC pct.1 pct.2     p_val_adj cluster      gene
RPS12     2.008629e-140  0.7256738 1.000 0.991 2.754633e-136       0     RPS12
RPS27     2.624075e-140  0.7242847 0.999 0.992 3.598656e-136       0     RPS27
RPS6      1.280169e-138  0.6742630 1.000 0.995 1.755623e-134       0      RPS6
RPL32     4.358823e-135  0.6121027 0.999 0.995 5.977689e-131       0     RPL32
RPS14     3.618793e-128  0.6179756 1.000 0.994 4.962812e-124       0     RPS14
CYBA      1.090337e-124 -1.5739355 0.661 0.914 1.495288e-120       0      CYBA
RPS25     5.612007e-124  0.7506936 0.997 0.975 7.696306e-120       0     RPS25
RPL9      2.513370e-118  0.7518221 0.994 0.972 3.446835e-114       0      RPL9
CD74      5.800646e-117 -2.6837317 0.677 0.904 7.955006e-113       0      CD74
RPL13     3.048839e-116  0.5679197 1.000 0.996 4.181178e-112       0     RPL13
RPL31     6.670396e-114  0.7342977 0.996 0.964 9.147782e-110       0     RPL31
RPS3      1.132020e-113  0.5981477 1.000 0.994 1.552453e-109       0      RPS3
HLA-DRB1  1.716574e-112 -3.0434037 0.102 0.588 2.354109e-108       0  HLA-DRB1
RPL3      3.522082e-112  0.6070706 0.997 0.995 4.830183e-108       0      RPL3
RPL21     1.455167e-109  0.6677666 0.997 0.991 1.995616e-105       0     RPL21
RPS3A     1.075148e-108  0.7754659 0.991 0.977 1.474457e-104       0     RPS3A
RPL30     6.585600e-108  0.6582449 0.997 0.980 9.031491e-104       0     RPL30
LDHB      1.963031e-107  1.0532589 0.901 0.594 2.692101e-103       0      LDHB
HLA-DRA   1.488495e-105 -3.7142006 0.260 0.668 2.041322e-101       0   HLA-DRA
HLA-DPB1  5.609579e-104 -2.8850339 0.195 0.635 7.692977e-100       0  HLA-DPB1
S100A4    9.359936e-104 -1.8001237 0.683 0.859  1.283622e-99       0    S100A4
RPS15A    2.966671e-103  0.6208797 0.997 0.983  4.068492e-99       0    RPS15A
RPLP2     6.976268e-101  0.5870425 1.000 0.990  9.567255e-97       0     RPLP2
HLA-DPA1   1.565361e-99 -2.7757683 0.172 0.609  2.146736e-95       0  HLA-DPA1
RPS27A     4.351761e-98  0.7361393 0.994 0.967  5.968004e-94       0    RPS27A
MALAT1     1.745862e-97  0.6755897 1.000 0.999  2.394275e-93       0    MALAT1
EEF1A1     4.229218e-96  0.5243268 0.994 0.991  5.799949e-92       0    EEF1A1
RPS13      2.677032e-93  0.6792220 0.984 0.962  3.671282e-89       0     RPS13
RPL13A     6.590289e-93  0.4783244 1.000 0.999  9.037923e-89       0    RPL13A
RPL27A     2.994498e-92  0.5368410 0.999 0.988  4.106655e-88       0    RPL27A
OAZ1       5.987877e-92 -1.3195331 0.818 0.946  8.211775e-88       0      OAZ1
LGALS1     7.190789e-92 -2.8670834 0.172 0.586  9.861448e-88       0    LGALS1
RPS18      1.812479e-91  0.5141226 1.000 0.996  2.485633e-87       0     RPS18
CLIC1      8.449825e-89 -1.4229960 0.346 0.753  1.158809e-84       0     CLIC1
RPL11      1.168775e-88  0.4568759 1.000 0.995  1.602857e-84       0     RPL11
RPS29      1.284058e-88  0.7439036 0.990 0.904  1.760958e-84       0     RPS29
HLA-DRB5   3.747079e-88 -2.4167021 0.040 0.463  5.138745e-84       0  HLA-DRB5
RPS28      5.222372e-88  0.6260073 0.991 0.969  7.161960e-84       0     RPS28
RPL35A     1.312573e-87  0.5607366 0.997 0.981  1.800062e-83       0    RPL35A
RPL23A     4.051184e-86  0.5649971 1.000 0.988  5.555793e-82       0    RPL23A
RPS23      1.014288e-84  0.5703270 0.997 0.977  1.390994e-80       0     RPS23
SH3BGRL3   1.874841e-84 -1.1250592 0.736 0.903  2.571158e-80       0  SH3BGRL3
CCR7       1.606796e-82  1.3300401 0.436 0.110  2.203560e-78       0      CCR7
TPT1       1.565057e-80  0.5806329 0.996 0.982  2.146320e-76       0      TPT1
SRGN       1.373584e-79 -1.6235500 0.343 0.691  1.883733e-75       0      SRGN
RPS20      2.705246e-79  0.6656386 0.987 0.943  3.709974e-75       0     RPS20
ACTB       9.343728e-79 -0.8878646 0.990 0.994  1.281399e-74       0      ACTB
CD3D       4.198081e-77  0.9519779 0.838 0.406  5.757249e-73       0      CD3D
TYROBP     5.411737e-77 -3.5452746 0.128 0.495  7.421656e-73       0    TYROBP
RPS4X      1.098372e-76  0.4956895 0.997 0.992  1.506307e-72       0     RPS4X
RPL5       1.941828e-76  0.5932938 0.989 0.969  2.663023e-72       0      RPL5
RPS16      7.993484e-74  0.5290406 0.996 0.982  1.096226e-69       0     RPS16
SAT1       5.721654e-73 -2.0045430 0.353 0.676  7.846677e-69       0      SAT1
FTH1       2.057708e-72 -1.9685786 0.984 0.992  2.821941e-68       0      FTH1
FCER1G     3.494393e-72 -2.8510248 0.095 0.453  4.792210e-68       0    FCER1G
S100A11    1.413537e-71 -1.9170808 0.271 0.618  1.938524e-67       0   S100A11
ARPC3      4.200093e-71 -1.0278629 0.654 0.876  5.760008e-67       0     ARPC3
RPL19      4.624807e-68  0.3826991 0.999 0.993  6.342460e-64       0     RPL19
RPL10      4.099133e-67  0.3462410 1.000 0.997  5.621551e-63       0     RPL10
RPL10A     7.537575e-67  0.5072973 0.991 0.985  1.033703e-62       0    RPL10A
FTL        1.564496e-66 -2.2975545 0.971 0.996  2.145550e-62       0       FTL
ANXA2      3.730564e-66 -1.5285999 0.166 0.535  5.116096e-62       0     ANXA2
RPSA       9.342163e-66  0.5666680 0.989 0.948  1.281184e-61       0      RPSA
TYMP       3.284469e-63 -2.2080492 0.102 0.437  4.504321e-59       0      TYMP
CST3       3.202518e-61 -3.5340070 0.168 0.481  4.391934e-57       0      CST3
RPS5       5.274845e-61  0.5009655 0.990 0.975  7.233923e-57       0      RPS5
RPS10      7.147870e-60  0.5255493 0.987 0.970  9.802589e-56       0     RPS10
RPL14      8.305738e-60  0.4181006 0.994 0.985  1.139049e-55       0     RPL14
MYL6       1.615021e-59 -0.9062210 0.753 0.888  2.214840e-55       0      MYL6
LST1       6.881383e-58 -2.7390873 0.128 0.439  9.437129e-54       0      LST1
S100A6     2.896255e-56 -1.3813321 0.712 0.808  3.971924e-52       0    S100A6
RPS15      5.461292e-56  0.3950233 0.999 0.989  7.489616e-52       0     RPS15
PSAP       8.147086e-56 -1.7114077 0.301 0.584  1.117291e-51       0      PSAP
RPL18      8.318644e-56  0.5063594 0.997 0.979  1.140819e-51       0     RPL18
NPM1       1.858081e-55  0.7281530 0.910 0.801  2.548173e-51       0      NPM1
RPLP0      2.609336e-55  0.5153455 0.984 0.939  3.578443e-51       0     RPLP0
RPL36      3.245240e-55  0.5218402 0.984 0.950  4.450522e-51       0     RPL36
GSTP1      6.054495e-55 -1.7127985 0.230 0.543  8.303134e-51       0     GSTP1
CD3E       2.316054e-54  0.8655242 0.726 0.399  3.176237e-50       0      CD3E
HLA-DQB1   2.905680e-54 -2.0844468 0.024 0.318  3.984850e-50       0  HLA-DQB1
CTSS       1.403283e-53 -1.8754835 0.331 0.591  1.924462e-49       0      CTSS
S100A9     1.846301e-53 -4.3769204 0.141 0.435  2.532018e-49       0    S100A9
HLA-DMA    8.217096e-53 -1.7120408 0.046 0.345  1.126893e-48       0   HLA-DMA
RPL4       1.757749e-51  0.5224684 0.990 0.957  2.410577e-47       0      RPL4
PFN1       1.408071e-50 -0.7061830 0.915 0.956  1.931028e-46       0      PFN1
GAPDH      1.477443e-50 -0.9423717 0.789 0.873  2.026165e-46       0     GAPDH
NOSIP      3.191555e-50  0.9992232 0.628 0.358  4.376899e-46       0     NOSIP
ARPC1B     6.218794e-50 -1.0270372 0.604 0.800  8.528455e-46       0    ARPC1B
COTL1      8.281971e-50 -1.6426851 0.475 0.685  1.135789e-45       0     COTL1
HLA-DQA1   1.816661e-49 -2.1570174 0.020 0.290  2.491369e-45       0  HLA-DQA1
PYCARD     1.869452e-49 -1.4766044 0.149 0.451  2.563766e-45       0    PYCARD
SPI1       2.279180e-49 -1.6238341 0.017 0.284  3.125667e-45       0      SPI1
LEF1       3.324866e-49  1.0526440 0.336 0.104  4.559722e-45       0      LEF1
RPL36A     2.403282e-48  0.5579918 0.956 0.890  3.295861e-44       0    RPL36A
ARPC2      2.663878e-48 -0.7342768 0.679 0.846  3.653242e-44       0     ARPC2
EEF1B2     2.099667e-47  0.5433143 0.940 0.869  2.879483e-43       0    EEF1B2
RPL22      2.861041e-47  0.6468201 0.901 0.780  3.923631e-43       0     RPL22
RPL7       8.797324e-47  0.5357066 0.977 0.968  1.206465e-42       0      RPL7
RPS8       1.804472e-46  0.4303529 0.997 0.984  2.474652e-42       0      RPS8
CD68       2.073468e-46 -1.6285734 0.030 0.292  2.843554e-42       0      CD68
LGALS3     1.177011e-45 -1.6725725 0.050 0.317  1.614153e-41       0    LGALS3
S100A8     8.669629e-45 -3.8982173 0.083 0.346  1.188953e-40       0    S100A8
RPL18A     1.106938e-44  0.3160600 0.999 0.992  1.518055e-40       0    RPL18A
RPS19      1.787996e-44  0.3524235 0.997 0.992  2.452057e-40       0     RPS19
RPS2       2.052882e-44  0.3153629 1.000 0.995  2.815322e-40       0      RPS2
PRKCQ-AS1  2.498572e-44  1.0271608 0.331 0.110  3.426542e-40       0 PRKCQ-AS1
FGR        2.642524e-44 -1.2814713 0.029 0.283  3.623957e-40       0       FGR
FCN1       5.700018e-44 -2.5350172 0.108 0.364  7.817005e-40       0      FCN1
MT-CO1     1.516102e-43 -0.4402246 0.990 0.996  2.079182e-39       0    MT-CO1
PRELID1    2.110668e-43 -1.1677078 0.263 0.545  2.894570e-39       0   PRELID1
CFD        2.987619e-43 -1.8522415 0.046 0.296  4.097221e-39       0       CFD
CEBPD      3.115612e-43 -1.4639507 0.037 0.292  4.272750e-39       0     CEBPD
RPL17      5.045664e-43  0.4359679 0.971 0.959  6.919624e-39       0     RPL17
PIK3IP1    7.614450e-43  0.9378420 0.438 0.185  1.044246e-38       0   PIK3IP1
LGALS2     1.060495e-42 -2.3137677 0.034 0.278  1.454363e-38       0    LGALS2
JUNB       2.031849e-42  0.6844529 0.914 0.904  2.786477e-38       0      JUNB
LYZ        2.307461e-42 -3.7006082 0.476 0.651  3.164452e-38       0       LYZ
LY86       7.859702e-42 -1.3249133 0.019 0.254  1.077880e-37       0      LY86
IFI30      9.640453e-42 -1.3725064 0.001 0.225  1.322092e-37       0     IFI30
GRN        9.942286e-42 -1.4990093 0.057 0.309  1.363485e-37       0       GRN
TMEM66     1.403899e-41  0.7213513 0.841 0.701  1.925307e-37       0    TMEM66
GLTSCR2    5.826768e-40  0.5634697 0.910 0.806  7.990830e-36       0   GLTSCR2
NPC2       8.380077e-40 -1.5214836 0.241 0.482  1.149244e-35       0      NPC2
IFITM3     1.032521e-39 -1.9324130 0.052 0.289  1.415999e-35       0    IFITM3
RPS26      1.048492e-39  0.4896995 0.958 0.900  1.437902e-35       0     RPS26
BTG1       1.700477e-39  0.5691811 0.938 0.860  2.332034e-35       0      BTG1
AP1S2      2.249470e-39 -1.3601929 0.108 0.366  3.084923e-35       0     AP1S2
FHIT       3.318236e-39  0.8686414 0.197 0.041  4.550628e-35       0      FHIT
CD7        7.678552e-39  0.6843976 0.572 0.294  1.053037e-34       0       CD7
SERPINA1   1.149394e-38 -1.5328694 0.069 0.311  1.576278e-34       0  SERPINA1
RPL38      1.167259e-38  0.5126637 0.941 0.876  1.600779e-34       0     RPL38
RAC1       1.783532e-38 -1.0841010 0.244 0.511  2.445935e-34       0      RAC1
RPL35      7.936208e-38  0.4140194 0.987 0.966  1.088372e-33       0     RPL35
SERF2      1.524437e-37 -0.6499521 0.782 0.883  2.090613e-33       0     SERF2
RPL24      7.315735e-37  0.4172226 0.977 0.944  1.003280e-32       0     RPL24
ATP6V0B    3.084259e-36 -1.1281389 0.148 0.406  4.229752e-32       0   ATP6V0B
RPS21      3.622080e-36  0.4985400 0.933 0.841  4.967321e-32       0     RPS21
GABARAP    1.611449e-35 -1.1373864 0.331 0.547  2.209941e-31       0   GABARAP
CTSH       2.711083e-35 -1.0464299 0.023 0.233  3.717980e-31       0      CTSH
CFP        3.621986e-35 -1.3076650 0.046 0.266  4.967191e-31       0       CFP
IL7R       3.894244e-35  0.7256307 0.597 0.333  5.340566e-31       0      IL7R
RPLP1      4.617630e-35  0.2967898 1.000 0.996  6.332618e-31       0     RPLP1
 [ reached 'max' / getOption("max.print") -- omitted 7928 rows ]
 ```
<br>

<h2>8. Repeat Analysis to Distinguish Subsets</h2>
	
<i>We repeated this analysis to identify marker genes distinguishing subsets within a cell-type.</i>

TODO: start on this once the 10X genomics data pipeline processes the custom data

<br>

<h2>9. Gene Ontology Enrichment Analysis</h2>

<i>To infer the functional relevance of sub-clusters, we performed gene ontology enrichment analyses on the top 50 differentially, expressed genes using Fisher’s Exact test as implemented in the topGO R package. For the enrichment analyses of the gene expression changes in astrocytes, our initial analysis revealed very few differentially expressed genes between the uninjured and 1dpi astrocytes, which we attributed to the low numbers of uninjured astrocytes captured. Therefore, we supplemented our uninjured astrocyte dataset with ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset. We also supplemented our uninjured OPC dataset in order to validate that our uninjured cells were more transcriptional similar to the external reference cells than to our injured cells.</i>

TODO: review documentation - https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

# to speed up testing time while experimenting, we can bypass the clustering step by saving the data with: pbmc <- readRDS("pdmc.RDS")



<br>

<h2>10. Differential Expression Tests</h2>
	
<i>To account for differences in sequencing depth between our dataset and the external dataset, we performed differential expression tests using MAST as implemented in Seurat. We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis.</i>





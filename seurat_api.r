library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)




# D A T A  H A N D L I N G  #####################################################################


# for speeding up workflow, we can save and load the Seurat object, without having to rerun clustering

save_data <- function(data) {
        saveRDS(data, file = "SeuratObject.RDS")
}


load_data <- function() {
        data <- readRDS("SeuratObject.RDS")
        return (data)
}


output_seurat_object_data_features <- function(seurat_object) {
        print('general info: ')
        print(seurat_object)
        print('dimensions: ')
        print(dim(x = seurat_object))
        print('row names: ')
        print(head(x = rownames(x = seurat_object)))
        print('column names: ')
        print(head(x = colnames(x = seurat_object)))
        print('a vector of names of associated objects: ')
        print(names(x = seurat_object))
}


get_slot_names <- function(seurat_object) {
	return (slotNames(seurat_object))
}


get_assays <- function(seurat_object) {
	# output
	# "RNA"  "pca"  "tsne"
	# TODO: double check these are all the available assays
	return (names(x = seurat_object))
}


get_rna <- function(seurat_object) {
	# output
	# Assay data with 13714 features for 2638 cells
	# Top 10 variable features:
 	# PPBP, DOK3, NFE2L2, ARVCF, YPEL2, UBE2D4, FAM210B, CTB-113I20.2, GBGT1, GMPPA
	return (seurat_object[['RNA']])
}


get_pca <- function(seurat_object) {
	# output
	# A dimensional reduction object with key tSNE_
 	# Number of dimensions: 2
 	# Projected dimensional reduction calculated: FALSE
 	# Jackstraw run: FALSE
	return (seurat_object[['tsne']])
}


get_assay_data <- function(seurat_object, slot_name, start_row, end_row, start_column, end_column) {
	# slot_name = 'scale.data'
	# output
	#                AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC
	# AL627309.1       -0.06547546    -0.10052277    -0.05804007
	# AP006222.2       -0.02690776    -0.02820169    -0.04508318
	# RP11-206L10.2    -0.03596234    -0.17689415    -0.09997719
	return (GetAssayData(object = seurat_object, slot = slot_name)[start_row:end_row, start_column:end_column])
}


get_available_cell_level_meta_data <- function(seurat_object) {
	# [1] "nGene"        "nUMI"         "orig.ident"   "percent.mito" "res.0.6"
	return (colnames(x = seurat_object[[]]))
}


get_multiple_cell_level_meta_data_values <- function(seurat_object) {
        # returns just the top with head
	# replace with parameters from get_available_cell_level_meta_data(), hardcoded for because we can specify multiple parameters
	# output
	#                nUMI percent.mito
	# AAACATACAACCAC 2421  0.030177759
	# AAACATTGAGCTAC 4903  0.037935958
	# AAACATTGATCAGC 3149  0.008897363
	# AAACCGTGCTTCCG 2639  0.017430845
	# AAACCGTGTATGCG  981  0.012244898
	# AAACGCACTGGTAC 2164  0.016643551
	return (head(x = seurat_object[[c('nUMI', 'percent.mito')]]))
}


get_mean_dispersion_and_dispersion_scaled <- function(seurat_object) {
	# output
	#                      mean dispersion dispersion.scaled
	# AL627309.1    0.013555659   1.432845        -0.6236875
	# AP006222.2    0.004695980   1.458631        -0.5728009
	# RP11-206L10.2 0.005672517   1.325459        -0.8356099
	# RP11-206L10.9 0.002644177   0.859264        -1.7556304
	# LINC00115     0.027437275   1.457477        -0.5750770
	# NOC2L         0.376037723   1.876440        -0.4162432	
	return (head(x = HVFInfo(object = seurat_object)))
}


get_variable_features <- function(seurat_object) {
	# note: VariableFeatures both accesses and sets the vector of variable features
	# output:
	# [1] "PPBP"   "DOK3"   "NFE2L2" "ARVCF"  "YPEL2"  "UBE2D4"
	return head(x = VariableFeatures(object = seurat_object))
}


get_standard_deviations <- function(seurat_object) {
	# output
	# [1] 5.666584 4.326466 3.952192 3.638124 2.191529 1.996551 1.877891 1.798251
        # [9] 1.766873 1.753684 1.731568 1.720525 1.718079 1.715879 1.707009 1.702660
        # [17] 1.697318 1.692549 1.686149 1.683967
	return (Stdev(object = seurat_object, reduction.use = 'pca'))
}




# S E U R A T  A P I  ###########################################################################


load_and_process_seurat_object <- function() {
	variable_genes_count = 10 # replace with 2000 later when custom data is processed
        pbmc.data <- Read10X(data.dir = '../Data/Sample/filtered_gene_bc_matrices/hg19')
	pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^M		  T-")
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

	return (seurat_object)
}





# so = analyze(load_and_process_seurat_object())
# save_data(so)


so = load_data()

# print(so)




output_seurat_object_data_features(so)


print(get_slot_names(so))



print(get_assays(so))


print(get_assay_data(so, 'scale.data', 1, 3, 1, 3))





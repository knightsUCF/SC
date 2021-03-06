library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)




# D A T A  H A N D L I N G  #####################################################################


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


get_tsne <- function(seurat_object) {
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
	# replace with parameters from get_available_cell_level_meta_data(), hardcoded for now because we can specify multiple parameters
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
	return (head(x = VariableFeatures(object = seurat_object)))
}


get_standard_deviations <- function(seurat_object) {
	# output
	# [1] 5.666584 4.326466 3.952192 3.638124 2.191529 1.996551 1.877891 1.798251
        # [9] 1.766873 1.753684 1.731568 1.720525 1.718079 1.715879 1.707009 1.702660
        # [17] 1.697318 1.692549 1.686149 1.683967
	return (Stdev(object = seurat_object, reduction.use = 'pca'))
}



# R N A ########################################################################################


get_rna <- function(seurat_object) {
	# output
	# Assay data with 13714 features for 2638 cells
	# Top 10 variable features:
 	# PPBP, DOK3, NFE2L2, ARVCF, YPEL2, UBE2D4, FAM210B, CTB-113I20.2, GBGT1, GMPPA
	return (seurat_object[['RNA']])
}


get_rna_dimensions <- function(seurat_object) {
    # output
    # [1] 13714  2638
    rna = seurat_object[['RNA']]
    return (dim(x = rna))
}


get_rna_row_data <- function(seurat_object) {
    # output
    # [1] "AL627309.1"    "AP006222.2"    "RP11-206L10.2" "RP11-206L10.9"
    # [5] "LINC00115"     "NOC2L"
    rna = seurat_object[['RNA']]
    return (x = rownames(x = rna))
}


# we are repeating ourselves here, but for now both method names will be useful to show how the data is structured

get_genes <- function(seurat_object) {
    # output
    # [1] "AL627309.1"    "AP006222.2"    "RP11-206L10.2" "RP11-206L10.9"
    # [5] "LINC00115"     "NOC2L"
    rna = seurat_object[['RNA']]
    return (x = rownames(x = rna))
}


get_rna_column_data <- function(seurat_object) {
    # output
    # [1] "AAACATACAACCAC" "AAACATTGAGCTAC" "AAACATTGATCAGC" "AAACCGTGCTTCCG"
    # [5] "AAACCGTGTATGCG" "AAACGCACTGGTAC"
    rna = seurat_object[['RNA']]
    return (x = colnames(x = rna))
}


get_gene_sequence <- function(seurat_object) {
    # output
    # [1] "AAACATACAACCAC" "AAACATTGAGCTAC" "AAACATTGATCAGC" "AAACCGTGCTTCCG"
    # [5] "AAACCGTGTATGCG" "AAACGCACTGGTAC"
    rna = seurat_object[['RNA']]
    return (x = colnames(x = rna))
}


get_rna_data_by_column_and_row <- function(seurat_object, start_row, end_row, start_column, end_column) {
    # output:
    # 3 x 3 sparse Matrix of class "dgCMatrix"
    #          AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC
    # AL627309.1                 .              .              .
    # AP006222.2                 .              .              .
    # RP11-206L10.2   
    rna = seurat_object[['RNA']]
    return (rna[start_row:end_row, start_column:end_column])
}


get_rna_data_from_specific_slot <- function(seurat_object, slot_name, start_row, end_row, start_column, end_column) {
    # slot name = 'scale.data'
    # output:
    #              AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC
    # AL627309.1       -0.06547546    -0.10052277    -0.05804007
    # AP006222.2       -0.02690776    -0.02820169    -0.04508318
    # RP11-206L10.2    -0.03596234    -0.17689415    -0.09997719
    rna = seurat_object[['RNA']]
    return (GetAssayData(object = rna, slot = 'scale.data')[start_row:end_row, start_column:end_column])
}


get_rna_available_column_meta_data <- function(seurat_object) {
    # output
    # [1] "mean"              "dispersion"        "dispersion.scaled"
    rna = seurat_object[['RNA']]
    return (colnames(x = rna[[]]))
}


get_rna_mean_dispersion_and_dispersion_scaled <- function(seurat_object) {
    # output
    #                     mean dispersion dispersion.scaled
    # AL627309.1    0.013555659   1.432845        -0.6236875
    # AP006222.2    0.004695980   1.458631        -0.5728009
    # RP11-206L10.2 0.005672517   1.325459        -0.8356099
    # can pull specific values with: (x = rna[[c('mean', 'dispersion')]]
    # also can turn into a name vector with: head(x = rna[['mean', drop = TRUE]])
    rna = seurat_object[['RNA']]
    return (x = HVFInfo(object = rna))
}


get_rna_variable_features <- function(seurat_object) {
    # output
    # [1] "PPBP"   "DOK3"   "NFE2L2" "ARVCF"  "YPEL2"  "UBE2D4"
    rna = seurat_object[['RNA']]
    return (x = VariableFeatures(object = rna))
}


set_rna_key <- function(seurat_object, key_name) {
    # output
    # "myRNA_"
    rna = seurat_object[['RNA']]
    Key(object = rna) <- key_name
    return (Key(object = rna))
}


access_rna_key <- function(seurat_object) {
    # output
    # "rna_"
    rna = seurat_object[['RNA']]
    return (Key(object = rna))
}


get_rna_feature <- function(seurat_object, feature_name) {
    # output
    # feature_name = 'rna_MS4A1'
    #               rna_MS4A1
    # AAACATACAACCAC  0.000000
    # AAACATTGAGCTAC  2.583047
    # AAACATTGATCAGC  0.000000
    # AAACCGTGCTTCCG  0.000000
    # AAACCGTGTATGCG  0.000000
    # AAACGCACTGGTAC  0.000000
    return (x = FetchData(object = seurat_object, vars.fetch = feature_name))
}







# S E U R A T  A P I  ###########################################################################



show_available_seurat_methods <- function() {
    # output
    # [: access expression data from the data slot
    # [[: access feature-level metadata
    # [[<-: add feature-level metadata
    # colMeans: calculate means across columns (cells) of any expression matrix within the Assay
    # colSums: calculate sums across columns (cells) of any expression matrix within the Assay
    # dimnames: get a list with row (feature) and column (cell) names
    # dim: get the number of features (in data) and cells in the Assay
    # GetAssayData: pull one of the expression matrices within the Assay
    # HVFInfo:
    # Key: get the key assigned to the Assay
    # Key<-: ...
    # merge: ...
    # RenameCells: ...
    # rowMeans: calculate means across rows (features) of any expression matrix within the Assay
    # rowSums: calculate sums across rows (features) of any expression matrix within the Assay
    # SetAssayData: add data to or replace one of the expresion matrices within the Assay
    # SubsetData: ...
    # VariableFeatures: pull the names of features designated as variable
    # VariableFeatures<-: assign a vector of features that are considered variable
    # WhichCells: ...

    print(utils::methods(class = 'Assay'))
}


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





# save once
# so = analyze(load_and_process_seurat_object())
# save_data(so)

# then reload to save time
so = load_data()

# access seurat object data features
# output_seurat_object_data_features(so)
# print(get_slot_names(so))
# print(get_assays(so))
# print(get_assay_data(so, 'scale.data', 1, 3, 1, 3))







print(get_genes(so))

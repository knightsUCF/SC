# Single Cell Analysis

<h3> Contents </h3>

1.0 Overview

1.1 Nature Background Article Notes

2.0 Seurat

2.1 Installing Seurat

2.2 Seurat Documentation

2.3 Seurat Tutorials


<br>

# 1.0 Overview

Goal: replicate results of "Single cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord"

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1

Background article:

https://www.nature.com/articles/s12276-018-0071-8


<h3> 1.1 Nature Background Article Notes </h3>

* "In this review, we will focus on technical challenges in single-cell isolation and library preparation and on computational analysis pipelines available for analyzing scRNA-seq data."


* "Currently, however, the majority of transcriptome analysis experiments continue to be based on the assumption that cells from a given tissue are homogeneous, and thus, these studies are likely to miss important cell-to-cell variability."

* "For example, the ability to find and characterize outlier cells within a population has potential implications for furthering our understanding of drug resistance and relapse in cancer treatment"

https://www.nature.com/articles/s12276-018-0071-8

TODO: review article

# 2.0 Seurat


https://satijalab.org/seurat

https://github.com/satijalab/seurat

<h3> 2.1 Installing Seurat </h3>

Enter these commands in R:

    install.packages("devtools")

    library(devtools)

    install_github("Displayr/flipPlots")

    library(flipPlots)

    remotes::install_github("satijalab/seurat", ref = "release/4.0.0")

    remotes::install_github("jlmelville/uwot")
    
    
 <h3> 2.2 Seurat Documentation </h3>
 
 https://satijalab.org/seurat
 
 <h3> 2.3 Seurat Tutorials </h3>
 
* Guided tutorial --- 2,700 PBMCs

* Multiple Dataset Integration and Label Transfer

* Weighted Nearest Neighbor Analysis

* Multimodal Reference Mapping

* Analysis of spatial datasets

* Mouse Cell Atlas, 250K cells

* Multimodal analysis

* Stimulated vs Control PBMCs

* SCTransform

* scATAC-seq + scRNA-seq integration

* Mixscape

https://satijalab.org/seurat/vignettes.html


<h3> 2.3 Guided Clustering Seurat Tutorial </h3>

https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

to get this example working, at least on a mac we have to change the line to:

    pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
    
when data folder is on Desktop:

    setwd("/Users/x/Desktop")
    
then we can get the Seurat object:


    library(dplyr)
    library(Seurat)
    library(patchwork)

    pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")


    # Initialize the Seurat object with the raw (non-normalized data).
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

    print(pbmc)

    # An object of class Seurat 
    # 13714 features across 2700 samples within 1 assay 
    # Active assay: RNA (13714 features, 0 variable features)






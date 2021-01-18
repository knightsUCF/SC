# Single Cell Analysis

<h3> Contents </h3>

<h4> 1.0 Overview </h4>

1.1 Goal

1.2 Background Article Notes

<h4> 2.0 Seurat </h4>

2.1 Installing Seurat

2.2 Seurat Documentation

2.3 Guided Clustering Seurat Tutorial (Processing a Seurat Object)

<h4> 3.0 Replicating Results </h4>

3.1 Steps

<h4> 4.0 Data </h4>

4.1 Importing Data into Seurat

4.2 Feature Counts

4.3 Installing Bioconductor for Feature Counts

4.4 Using Feature Counts

4.5 Samtools

4.6 Big Data

<br>

# 1.0 Overview

<h3> 1.1 Goal </h3>

Goal: replicate results of "Single cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord"

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1

Background article:

https://www.nature.com/articles/s12276-018-0071-8


<h3> 1.2 Background Article Notes </h3>

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


<h3> 2.3 Guided Clustering Seurat Tutorial (Processing a Seurat Object) </h3>

https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

to get this example working, at least on a Mac we have to change the line to:

    pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
    
when data folder is on Desktop:

    setwd("/Users/x/Desktop")
    
then we can get the Seurat object:

```R
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
```

# 3.0 Replicating Results

Here is a place to start replicating the results of the single cell analysis.

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1.full

"To generate the full SCI dataset, all samples were processed and combined using Seurat v353. After filtering each sample count matrix for genes that were expressed in at least 10 cells, each dataset was independently normalized and scaled using the SCTransform function. SCTransform performs a variance-stabilizing transformation in which genes are grouped according to mean expression in order to smoothen model parameters for negative binomial regression54. To remove cell-cycle genes as a confounding source of variation, mitochondrial percentage and cell cycle scores based on the expression of canonical G2M and S phase markers were computed for each cell. Cell cycle genes were provided through the Seurat tutorial. These score values were then used as input for the “vars.to.regress” argument in the SCTransfrom() function. This operation generates a “corrected” expression matrix by building a regression model on these variables for each gene.

To identify shared and unique molecular cell-types across datasets and time-points, sample expression matrices were batch-corrected using Seurat’s Data Integration workflow55. In brief, data integration uses canonical correlation analysis between two samples to identify mutual nearest neighbors.

These “anchors” are then used to generate a corrected gene expression matrix based on the consistency of the anchors between cells, effectively performing a batch correction.

For the full SCI dataset, the 2000 most variables genes were used as input for the “anchor. features” argument of the FindIntegrationAnchors() function, where the variance of a gene was measured as the residual divided by the expected variance under the SCTransfrom() model. This resulted in a single, batch-corrected expression matrix for containing all cells.

Dimensional Reduction, Clustering, and Differential Gene expression testing

In order to measure the greatest gene expression variation among all the SCI cells, we first performed PCA on the batch-corrected expression matrix for the top 2000 variable genes taken from above. The top 15 principal components (PCs) were selected based on the “elbow” point heuristic in a scree plot which quantifies the contribution of variance by each principal component. Using these components, a nearest-neighbor graph and shared-nearest-neighbor graph were generated with “k.neighbors” set to 20 by default. To visualize the cells, we generate a UMAP plot with default Seurat parameters using cell coordinates in PCA-space using the top 15 PCs. In order to cluster the cells based on similarity of expression, we ran the FindClusters() function on the shared-nearest-neighbor graph with default parameters.

To infer the functional relevance of sub-clusters, we performed gene ontology enrichment analyses on the top 50 differentially expressed genes using Fisher’s Exact test as implemented in the topGO R package48. For the enrichment analyses of the gene expression changes in astrocytes, our initial analysis revealed very few differentially expressed genes between the uninjured and 1dpi astrocytes, which we attributed to the low numbers of uninjured astrocytes captured. Therefore, we supplemented our uninjured astrocyte dataset with ACNT1 and ACNT2 astrocyte data from the previously published mouse CNS single-cell atlas dataset50. We also supplemented our uninjured OPC dataset in order to validate that our uninjured cells were more transcriptional similar to the external reference cells than to our injured cells. To account for differences in sequencing depth between our dataset and the external dataset, we performed differential expression tests using MAST as implemented in Seurat56. We used all differentially expressed genes (p_val_adj < 0.001) as input for gene ontology analysis."


<h3> 3.1 Steps </h3>

Let's start with these steps

steps to "generate the full SCI dataset"

process all samples and combine with Seurat

filter each sample count matrix for genes that were expressed in at least 10 cells

normalize and scale each data set with the SCTransform function ("performs variance stabiling transformation in which genes are grouped according to mean expression in order to smoothen model parameters for negative binomial regression"

removing cell-cycle genes as a confounding source of variation: "mitochondrial percentage and cell cycle scores based on the expression of canonical G2M and S phase markers were computed for each cell (Cell cycle genes were provided through the Seurat tutorial)"

These score values were then used as input for the “vars.to.regress” argument in the SCTransfrom() function. (This operation generates a “corrected” expression matrix by building a regression model on these variables for each gene.)

To identify shared and unique molecular cell-types across datasets and time-points, sample expression matrices were batch-corrected using Seurat’s Data Integration workflow (In brief, data integration uses canonical correlation analysis between two samples to identify mutual nearest neighbors.)

These “anchors” are then used to generate a corrected gene expression matrix based on the consistency of the anchors between cells, effectively performing a batch correction.


# 4.0 Data

"In this study, we used sc-RNAseq to generate a single cell transcriptomics dataset of virtually all cell types that comprise the uninjured and injured spinal cord at 1, 3, and 7dpi. From this dataset, we were able to obtain unique molecular signatures of multiple cell types as well as their subpopulations present throughout the acute injury phase."

"These cell types were grouped into three categories, namely myeloid, vascular, and macroglia, for further analysis as described below."

http://mousebrain.org/

http://mousebrain.org/downloads.html

https://www.ncbi.nlm.nih.gov/biosample/SAMN17339713/

SRA: SRS8040434


"To assess the cellular heterogeneity among all cell populations at the injury site, we obtained a total of 51,843 cells from uninjured and 1, 3, and 7dpi tissue, which resulted in a total of 15 distinct clusters when visualized on a UMAP plot"

https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6854062

(click on run, and then data access)

<h3> 4.1 Importing Data into Seurat </h3>

"Seurat objects and their inputs don't actually contain the nucleotide sequence information you are looking for. If you are looking for the full fasta sequence for a specific gene or transcript isoform you can find easily find those from ensembl or NCBI browsers. If you are looking for the sequences for the reads that match those transcripts in your samples then you'll need to dig into the BAM file (assuming cell ranger output). You visualize those reads and where they map using interactive tools like IGV."

BAM files not supported: https://github.com/satijalab/seurat/issues/2488

"This is not supported in Seurat, you will need to create a count matrix before loading data into Seurat."

"how to convert my barcoded bam file to a count matrix."

https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html


<h3> 4.2 Feature Counts </h3>

"Output of counting = A count matrix, with genes as rows and samples are columns"

"we will be using the featureCounts tool to get the gene counts. We picked this tool because it is accurate, fast and is relatively easy to use."

"can be used to count both RNA-seq and genomic DNA-seq reads."

"featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations"

http://bioinf.wehi.edu.au/featureCounts/

<h3> 4.3 Installing Bioconductor for Feature Counts </h3>

http://bioconductor.org/install/

In R type:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")


BiocManager::install(c("Rsubread"))
```

<h3> 4.4 Using Feature Counts </h3>

Input for counting = multiple BAM files + 1 GTF file

Simply speaking, the genomic coordinates of where the read is mapped (BAM) are cross-referenced with the genomic coordinates of whichever feature you are interested in counting expression of (GTF), it can be exons, genes or transcripts.

(srun is used to submit a job for execution in real time)

Start an interactive session with 4 cores:

    $ srun --pty -p short -t 0-12:00 -c 4 --mem 8G --reservation=HBC /bin/bash
    
TODO: 8 GB Mac machine crashing, need to reformat 32 GB ram Windows machine to Linux

<h3> 4.5 Samtools </h3>

TODO: install Samtools

<h3> 4.6 Big Data Todos </h3>

Research:

- HPC cluster

- Google life sciences api (bioinformatics tools built with docker containers, unlimited computation)

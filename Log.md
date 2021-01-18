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



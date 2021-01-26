# Analysis



<h2> Introduction </h2>


Our goal is to replicate the results of the single cell paper using the Seurat library, and R.

<br>

<i>"Single cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord"</i>

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1.full

<i>Seurat</i>

https://satijalab.org/seurat/


<h2> Outline </h2>

We are using this outline to replicate the results in the study: 

https://github.com/knightsUCF/SC/blob/main/Outline.md

We are starting on the section: "Dimensional Reduction, Clustering, and Differential Gene expression testing"

TODO: genomic data sequencing might require starting on the data preprocessing section. Currently we are using the toy data given by the Seurat tutorial.


<h2> Libraries </h2>

We will need to import the following libraries:

```R
library(dplyr)
library(Seurat)
library(patchwork)
```

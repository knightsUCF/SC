# Analysis

<h3> Introduction </h3>

Our goal is to replicate the results of the single cell paper using the Seurat library, and R.

"Single cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord": https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1.full

Seurat: https://satijalab.org/seurat/

---

<h3> Outline </h3>

We are using this outline to replicate the results in the study: 

https://github.com/knightsUCF/SC/blob/main/Outline.md

We are starting on the section: "Dimensional Reduction, Clustering, and Differential Gene expression testing"

TODO: genomic data sequencing might require starting on the data preprocessing section. Currently we are using the toy data given by the Seurat tutorial.

---

<h3> Libraries </h3>

We will need to import the following libraries:

```R
library(dplyr)
library(Seurat)
library(patchwork)
```

# TopGO Gene List Format

https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/



# Clusters Overview

https://github.com/ucdavis-bioinformatics-training/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/blob/master/day3/scRNA_Workshop-PART6.Rmd


# Markers / Genes

"A genetic marker is a gene or DNA sequence with a known location on a chromosome that can be used to identify individuals or species. ... A genetic marker may be a short DNA sequence, such as a sequence surrounding a single base-pair change (single nucleotide polymorphism, SNP), or a long one, like minisatellites.

Genetic marker - Wikipedia"

# Seurat / TopGO Workflow

https://rpubs.com/kshekhar/349874

https://bioinformatics.stackexchange.com/questions/5225/script-to-allow-gene-set-enrichment-analysis-of-10x-genomics-data-in-r

"Seruat will give you a list of genes which it thinks are upregulated in a particular cluster. Look at the functions that talk about marker genes - these functions basically do a DE analysis of the genes in one cluster compared to the others.

Then take that list and feed it to any standard GO analysis tool. Have a look at the topGO topKEGG and geneSetTest functions in the limma package, the GOStats package and the gsea. All should be suitable. The GOSeq package is designed to compensate for gene length bias in RNA-seq experiments. As 10X only samples 3'-tags, there shouldn't be any gene-length bias in the data, so this shouldn't be an issue."


# TopGO with Seurat

https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day3/scRNA_Workshop-PART6.html


# Custom Python Environment

https://help.dreamhost.com/hc/en-us/articles/115000695551-Installing-and-using-virtualenv-with-Python-3

https://support.rstudio.com/hc/en-us/articles/360023654474-Installing-and-Configuring-Python-with-RStudio


# Selecting Mirror

When a package, such as UMAP is not working, we can try installing from another mirror:

chooseCRANmirror()

# Seurat Documentation

https://cran.r-project.org/web/packages/Seurat/Seurat.pdf

# Check Dots

"Seurat.checkdotsFor functions that have ...  as a parameter, this controls the behavior when anitem isn’t used. Can be one of warn, stop, or silent"

# Plan

- while waiting for the dependency issue to be resolved, continue

- we have to continue on the Mac, so check one more time if somehow we can resolve this, if not we just have to keep moving, while not getting our feet stuck in time sinks, just yet

- after we exhaust the other steps we can determine what is left in terms of fixing dependencies

- upgrading to latest Ubuntu, not enough space: https://www.zdnet.com/article/how-to-upgrade-from-ubuntu-linux-16-04-to-18-04/

- entering these two things into the R console, will show where the error is coming from:

options(error=recover) 

options(show.error.locations=TRUE)

Error in CheckDots(...) : argument is missing, with no default

Enter a frame number, or 0 to exit   

1: source("~/Documents/SC/sc.r")

2: withVisible(eval(ei, envir))

3: eval(ei, envir)

4: eval(ei, envir)

5: sc.r#37: FindNeighbors(pbmc, reduction = "pca", dims = 1:10, assay = NULL, features = NULL, k.param = 20, compute.SNN = TRUE, prune.SNN = 1/15, nn.me

6: FindNeighbors.Seurat(pbmc, reduction = "pca", dims = 1:10, assay = NULL, features = NULL, k.param = 20, compute.SNN = TRUE, prune.SNN = 1/15, nn.meth

7: CheckDots(...)

Selection: 7

Called from: eval(substitute(browser(skipCalls = skip), list(skip = 7 - which)), 
    envir = sys.frame(which))
    

     CheckDots(...)
      if (!is.null(x = dims)) {
        assay <- DefaultAssay(object = object[[reduction]])
        data.use <- Embeddings(object = object[[reduction]])
        if (max(dims) > ncol(x = data.use)) {
          stop("More dimensions specified in dims than have been computed")
        }
        data.use <- data.use[, dims]
        neighbor.graphs <- FindNeighbors(object = data.use, 
          k.param = k.param, compute.SNN = compute.SNN, prune.SNN = prune.SNN, 
          nn.method = nn.method, annoy.metric = annoy.metric, 
          nn.eps = nn.eps, verbose = verbose, force.recalc = force.recalc, 
          ...)
      }



# Set Working Directory

When opening project set the working directory. #TODO later find out how we can automatically set up this working directory.

# DPLYR

Dependency not installing on Linux, works on Mac.

https://cran.r-project.org/web/packages/dplyr/index.html


# Cran Servers

https://cran.r-project.org/mirrors.html


# General Plan

- try to complete the rest of part 1

- investigate the genomics data pipeline


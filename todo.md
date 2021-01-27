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


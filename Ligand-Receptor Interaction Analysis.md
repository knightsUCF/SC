# Ligand-Receptor Interaction Analysis

1. To infer potential ligand-receptor interactions between two cell-types, we adapted the method used in CellPhoneDB15. We first pulled a reference list of human ligand-receptor pairs published previously51 and converted the genes into mouse orthologs using the Ensembl biomaRt package49.

2. We defined the ligand-receptor score as the mean of the average log-normalized expression of the receptor gene in one cell-type and the average log-normalized expression of the ligand gene in a second cell-type.

3. To identify enriched ligand-receptor interactions, we applied a permutation test to identify interactions scores that are enriched in a specific <ligand cell A, receptor cell B, time-point> combination.

4. For each of 1000 permutations, we randomly shuffled the cell-type and time-point labels and calculated an interaction scores for all possible <ligand cell, receptor cell, time-point> combinations. Repeating this 1000 times generated a null distribution of interaction scores for each ligand-receptor pair.

5. We compared the interaction scores of the actual (ligand cell A, receptor cell B, time-point) labels to the null distribution and calculated p-values as the proportion of null scores which are equal to or greater than the actual interaction score.

<br>

<h2> 1. Getting Mouse Orthologs for Human Ligand-Receptor Pairs </h2>

https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/

Reference paper they used: https://www.nature.com/articles/s41586-018-0698-6

We want to get the gene ID.

<i>"Linking Ensembl and Uniprot identification

We assigned to the custom-curated interaction list all the Ensembl gene identifications by matching information from Uniprot and Ensembl by the gene name."</i>


<h3> Ensembl Orthologous </h3>

![Ensembl Orthologous](https://github.com/knightsUCF/SC/blob/main/charts/Ensembl%20Orthologous.png)

<br>

<h3> Ensembl Export Settings </h3>

![Ensembl Export Settings](https://github.com/knightsUCF/SC/blob/main/charts/Ensembl%20Export%20Settings.png)


TODO: review export settings from Ensembl

<br>

<h3> Mouse Orthologous Data </h3>

While we are reviewing the paper for the proper gene IDs, we pulled a random mouse orthologous gene to get the data process started: https://github.com/knightsUCF/SC/blob/main/mouse%20orthologous.txt

<br>

<h2>2. Defining the Ligand-Receptor Score</h2

<i>We defined the ligand-receptor score as the mean of the average log-normalized expression of the receptor gene in one cell-type and the average log-normalized expression of the ligand gene in a second cell-type.</i>

```R
print(head(AverageExpression(object = pbmc)))
```

<h3>Output</h3>

TODO: verify - https://github.com/satijalab/seurat/issues/193

```
$RNA
                               0            1            2            3
AL627309.1          6.128664e-03 5.927264e-03 4.854338e-02 0.000000e+00
AP006222.2          0.000000e+00 8.206078e-03 1.088471e-02 0.000000e+00
RP11-206L10.2       7.453092e-03 0.000000e+00 0.000000e+00 2.065031e-02
RP11-206L10.9       0.000000e+00 0.000000e+00 1.050116e-02 0.000000e+00
LINC00115           1.911893e-02 2.469048e-02 3.753737e-02 3.888541e-02
NOC2L               4.974632e-01 3.598115e-01 2.725375e-01 5.865349e-01
...

                               4            5            6            7
AL627309.1          2.054586e-02 0.000000e+00   0.00000000   0.00000000
AP006222.2          1.191488e-02 0.000000e+00   0.00000000   0.00000000
RP11-206L10.2       0.000000e+00 0.000000e+00   0.00000000   0.08462847
RP11-206L10.9       0.000000e+00 1.200008e-02   0.00000000   0.00000000
LINC00115           1.948277e-02 1.469374e-02   0.05855154   0.00000000


                             8
AL627309.1            0.0000000
AP006222.2            0.0000000
RP11-206L10.2         0.0000000
RP11-206L10.9         0.0000000
LINC00115             0.0000000
NOC2L                 0.0000000
KLHL17                0.0804829
PLEKHN1               0.0000000
RP11-54O7.17          0.0000000
HES4                  0.1609658
...

```

<br>

<h2>3. Identifying enriched ligand-receptor interactions</h2>

<i>To identify enriched ligand-receptor interactions, we applied a permutation test to identify interactions scores that are enriched in a specific <ligand cell A, receptor cell B, time-point> combination.</i>

TODO: get Ubuntu docker

https://rdrr.io/github/kendomaniac/rCASC/man/seuratPermutation.html

```R
 system("wget http://130.192.119.59/public/section4.1_examples.zip")
 unzip("section4.1_examples.zip")
 setwd("section4.1_examples")
 system("wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz")
 system("gzip -d Homo_sapiens.GRCh38.94.gtf.gz")
 system("mv Homo_sapiens.GRCh38.94.gtf genome.gtf")
 scannobyGtf(group="docker", file=paste(getwd(),"bmsnkn_5x100cells.txt",sep="/"),
             gtf.name="genome.gtf", biotype="protein_coding", 
             mt=TRUE, ribo.proteins=TRUE,umiXgene=3)
 
 seuratBootstrap(group="docker",scratch.folder="/data/scratch/",
      file=paste(getwd(), "annotated_bmsnkn_5x100cells.txt", sep="/"), 
      nPerm=160, permAtTime=8, percent=10, separator="\t",
      logTen=0, pcaDimensions=6, seed=111)
```
<br>

<h2>4. Calculating interaction scores</h2>

<i>For each of 1000 permutations, we randomly shuffled the cell-type and time-point labels and calculated an interaction scores for all possible <ligand cell, receptor cell, time-point> combinations. Repeating this 1000 times generated a null distribution of interaction scores for each ligand-receptor pair.</i>
  
 TODO
 
 <br>
 
 <h2>5. Comparing interaction scores</h2>
 
 <i>We compared the interaction scores of the actual (ligand cell A, receptor cell B, time-point) labels to the null distribution and calculated p-values as the proportion of null scores which are equal to or greater than the actual interaction score.</i>
 
 TODO

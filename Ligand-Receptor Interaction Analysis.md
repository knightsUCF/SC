# Ligand-Receptor Interaction Analysis

1. To infer potential ligand-receptor interactions between two cell-types, we adapted the method used in CellPhoneDB15. We first pulled a reference list of human ligand-receptor pairs published previously51 and converted the genes into mouse orthologs using the Ensembl biomaRt package49.

2. We defined the ligand-receptor score as the mean of the average log-normalized expression of the receptor gene in one cell-type and the average log-normalized expression of the ligand gene in a second cell-type.

3. To identify enriched ligand-receptor interactions, we applied a permutation test to identify interactions scores that are enriched in a specific <ligand cell A, receptor cell B, time-point> combination.

4. For each of 1000 permutations, we randomly shuffled the cell-type and time-point labels and calculated an interaction scores for all possible <ligand cell, receptor cell, time-point> combinations. Repeating this 1000 times generated a null distribution of interaction scores for each ligand-receptor pair.

5. We compared the interaction scores of the actual (ligand cell A, receptor cell B, time-point) labels to the null distribution and calculated p-values as the proportion of null scores which are equal to or greater than the actual interaction score.



<h2> 1. Getting Mouse Orthologs for Human Ligand-Receptor Pairs </h2>

https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/

Reference paper they used: https://www.nature.com/articles/s41586-018-0698-6

We want to get the gene ID.

<i>"Linking Ensembl and Uniprot identification

We assigned to the custom-curated interaction list all the Ensembl gene identifications by matching information from Uniprot and Ensembl by the gene name."</i>

While we are reviewing the paper for the proper gene IDs, we pulled a random mouse orthologous gene to get the data process started: https://github.com/knightsUCF/SC/blob/main/mouse%20orthologous.txt

TODO: review these export settings from Ensembl:



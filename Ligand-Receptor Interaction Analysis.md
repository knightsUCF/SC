# Ligand-Receptor Interaction Analysis

To infer potential ligand-receptor interactions between two cell-types, we adapted the method used in CellPhoneDB15. We first pulled a reference list of human ligand-receptor pairs published previously51 and converted the genes into mouse orthologs using the Ensembl biomaRt package49.

We defined the ligand-receptor score as the mean of the average log-normalized expression of the receptor gene in one cell-type and the average log-normalized expression of the ligand gene in a second cell-type.

To identify enriched ligand-receptor interactions, we applied a permutation test to identify interactions scores that are enriched in a specific <ligand cell A, receptor cell B, time-point> combination.

For each of 1000 permutations, we randomly shuffled the cell-type and time-point labels and calculated an interaction scores for all possible <ligand cell, receptor cell, time-point> combinations. Repeating this 1000 times generated a null distribution of interaction scores for each ligand-receptor pair.

We compared the interaction scores of the actual (ligand cell A, receptor cell B, time-point) labels to the null distribution and calculated p-values as the proportion of null scores which are equal to or greater than the actual interaction score.

# Outline

Numbered outline of research steps.

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1.full



<h3> Abstract section </h3>

One way to investigate cellular responses to injury is through transcriptomic profiling of cell types that comprise the injury site.


1) Used sc-RNAseq generate a single cell transcriptomics dataset of virtually all cell types that comprise the uninjured and injured spinal cord at 1, 3, and 7dpi

2) From this dataset, we were able to obtain unique molecular signatures of multiple cell types as well as their subpopulations present throughout the acute injury phase. 

3) By assessing expression of ligandreceptor pairs on different cell types, we were able to gain insight into potential signaling relationships that mediat angiogenesis, gliosis, and fibrosis.


<h3> Results section </h3>

1) To assess the cellular heterogeneity among all cell populations at the injury site, we obtained a total of 51,843 cells from uninjured and 1, 3, and 7dpi tissue, which resulted in a total of 15 distinct clusters when visualized on a UMAP plot (Fig. 1a, b).

2) These 15 clusters represented all major cell types that are known to comprise the SCI site including microglia, monocytes, macrophages, neutrophils, dendritic cells, astrocytes, oligodendrocytes, OPCs, neurons, fibroblasts, pericytes, ependymal cells, and endothelial cells. These cell types were grouped into three categories, namely myeloid, vascular, and macroglia, for further analysis as described below.

3) Neurons were excluded from the analysis because they were not expected to survive our dissociation protocol and thus subject to a selection bias.

4) Lymphocytes were also excluded because they are primarily involved in autoimmunity after SCI10, which is beyond the scope of this study. Cell types pertaining to each cluster were identified using annotated lineage markers (Extended Data Fig. 1).

5) The highest differentially expressed genes (DEGs) provided a unique molecular signature for each cell type (Fig. 1c, d), which in most cases were different from canonical markers used in the literature.

6) For example, the highest DEGs in OPCs were non-canonical genes such as tnr (Tenascin-R) and lhfpl3 (lipoma HMGIG fusion partner), which displayed better specificity than canonical OPC markers such as pdgfra and cspg4 that were expressed in multiple cell types (Extended Data Fig. 1).

7)  Interestingly, while certain marker genes were expressed both before and after injury, others such as postn in fibroblasts, changed expression in an injury-dependent manner. 

8) While postn was not expressed in uninjured fibroblasts, there was a graded increase as injury progressed (Fig. 1c). We validated this by genetic lineage tracing in PostnEYFP mice, which showed EYFP cells present in the fibrotic scar and overlying meninges but absent in surrounding spinal cord tissue. (Extended Data Fig. 2). 

9) Taken together, our analysis of DEGs between major cell types uncovers highly specific molecular identifiers, many of which are non-canonical and display temporal specificity.


<h3> Myeloid analysis reveals temporal changes in macrophage and microglial subtypes </h3>

1) To determine the heterogeneity within the myeloid population, clustering analysis was performed on myeloid cells and visualized on a separate UMAP (Fig. 2a, b), which revealed two large clusters corresponding to microglia and peripherally-derived myeloid cells as identified by annotated markers (Extended Data Fig. 3). 

2) We identified 6 microglial subtypes.

3) Homeostatic microglia were identified based on its expression of several annotated markers of steady-state microglia, such as p2ry12, siglech, and tmem119 (Fig. 2).

4) We identified four DAM subtypes, which were labeled DAM-A to D. DAM-A was identified by high expression of the low density lipoprotein receptor msr1 and low expression of the purinergic receptor p2ry12, and comprised 100% of microglia present at 1dpi. DAM-B and C expressed moderate level of p2ry12 and low level of msr1, and their expression profile was similar to homeostatic microglia, suggesting an intermediate state as some DAM-A revert back to a more homeostatic state. DAM-D had low levels of both msr1 and p2ry12, and was best identified by high expression of the growth factor igf1 and the cholesterol-binding protein apoe. 

5) The last microglia subtype was identified as dividing microglia due to their high expression of cell cycle-related genes (mki67 and top2a). Dividing microglia shared several DEGs with DAM-D such as igf1, spp1, and fabp5 (Fig. 2e), and the fact that expansion of DAM-D from 3 to 7 dpi coincides with reduction in dividing microglia (Fig. 2c) suggests that dividing microglia may give rise to DAM-D after injury. A putative model depicting the relationship between DAM subtypes is illustrated in Extended Data Fig. 4a.

6) We used flow cytometry to validate the presence of DAMs in vivo (gating strategy in Extended Data Fig 5). 

7) Microglia were gated on P2ry12 and Msr1 expression based on our sequencing data (Fig. 2d). As expected, over 90% of microglia present in the uninjured spinal cord were in the homeostatic state (P2ry12hi/Msr1lo), and this decreased to 20% at 1dpi (Fig. 2f, g).

8) At 1dpi, there was a large increase in Msr1hi microglia, consistent with the appearance of DAM-A. However, the majority of Msr1hi microglia were also P2ry12hi, representing a transition state between homeostatic and DAM-A, which are expected to be P2ry12lo/Msr1hi based on our sequencing results. 

9)  At 7dpi, we observed a significant decrease in Msr1hi microglia and a partial return of homeostatic microglia, which this was not observed in our sequencing data.

10) Taken together, our flow cytometry data support the appearance of DAM-A microglia subtype after SCI in vivo, but the temporal effects are more graded than those predicted from our sequencing data perhaps due to a delay in manifestation of gene expression changes at the protein level.

11) The peripherally-derived myeloid cluster revealed two monocyte and two macrophage subtypes in addition to border-associated macrophages and dendritic cells. The two monocyte subtypes corresponded to the classical Ly6Chi and the non-classical Ly6Clo monocytes (Fig. 3c) with the latter representing the largest population at 1dpi. 3 and 7dpi were dominated by the two macrophage subtypes, Macrophage-A and Macrophage-B respectively. 

12) Both subtypes expressed the lysosomal gene cd63, but were distinguished by preferential expression of heme oxygenase hmox in Macrophage-A and apoe in Macrophage-B. 

13) In conclusion, our data reveals the presence of multiple cellular states in the monocyte-macrophage lineage that display temporal progression toward a more pro-inflammatory state.

14) To validate the presence of Macrophage-A and Macrophage-B subtypes in vivo, we first isolated macrophages based on CD63hi expression (gating strategy in Extended Data Fig. 6). 

15) Further separation on ApoE and CD11b expression revealed two distinct clusters that were consistent with Macrophage-A (CD63hi/ApoElo/CD11bmed) and Macrophage-B (CD63hi/ApoEhi/CD11bhi)subtypes (Fig. 3d, e).

16) Consistent with our sequencing data, the monocytes were the predominant myeloid populations at 1 dpi, and subsequently decreased at 7dpi, whereas Macrophage-A (ApoElo) and Macrophage-B (ApoEhi) were the most represented macrophage subtypes at 1 and 7dpi, respectively (Fig. 3e).

<h3> Vascular heterogeneity analysis identifies tip cell dynamics </h3>

1. To determine the heterogeneity of vascular cells, clustering analysis was performed only on the vascular clusters (endothelial cells, fibroblasts, and pericytes from Fig. 1) and visualized on a separate UMAP, which revealed endothelial cell subtypes on one side and perivascular mural cells on the other. 

2. We identified each cluster using annotated markers from a previous sc-RNAseq study of the brain vasculature13.

3. We identified fibroblasts, pericytes, and vascular smooth muscle cells (VSMC) as three distinct populations of perivascular mural cells that were identified based on their expression of col1a1, kcnj8, and acta2 respectively (Extended Data Fig. 7b). 

4. We also identified an unknown vascular subtype (U-Vascular) that clustered with mural cells due to its molecular similarity with pericytes, but also expressed endothelial cell markers (Extended Data Fig. 7c).

5. Next, we identified an arterial, a venous, and two capillary subtypes based on annotated markers13. Arterial endothelial cells were identified by expression of gkn3 and stmn2, whereas venous endothelial cells were identified by slc38a5 and icam1. Capillary endothelial cells were identified by the expression of general endothelial cell markers ly6a and cldn5, combined with the lack of selective arterial and venous markers (Fig. 4d, Extended Data Fig. 7c). The fifth endothelial cluster was identified as tip cells based on their expression of the canonical marker apln.

6. Whereas the other endothelial subtypes did not show large temporal changes in proportion, the proportion of tip cells increased significantly at 1dpi, and decreased to basal levels by 7dpi (Fig. 4c).

7. This tip cell temporal profile was validated using in situ hybridization for apln combined with immunostaining for podocalyxin as an endothelial cell marker (Fig 4f).

8. In the uninjured cord, apln transcripts were detected scattered throughout the gray matter, consistent with previous reports of apln expression in neurons14. 

9. However, we did not detect any endothelial cells expressing apln in uninjured tissue sections, suggesting that tip cells detected in uninjured tissue in our sc-RNAseq data was most likely due to tissue processing. 

10. Taken together, our analysis identified all known major vascular cell types at the injury site, including previously undescribed tip cell molecular and temporal profiles.

<h3> Cellular interactions via angiopoietin and VEGF signaling during angiogenesis </h3>

1) To gain insight into mechanisms of cellular interactions during angiogenesis after SCI, we adapted CellPhoneDB15 to calculate “interaction scores” based on the average expression levels of a ligand and its receptor between two cells (Fig. 5a). We focused on angiopoietin (Angpt) and vascular endothelial growth factor (Vegf) signaling due to their well-known roles in angiogenesis. The highest interaction score was detected between tip cells and endothelial cells along the Angpt2-Tie2 pathway at 1dpi (Fig. 5b). Expression analysis showed that Angpt2 is expressed most highly by tip cells at 1dpi followed by a gradual decrease over the next 7 days, and Tie2 (and Tie1) is expressed in endothelial cells at all time points (Fig. 5c). Interestingly, Angpt1 signaling to endothelial/tip cells shifted from VSMC at 1dpi to astrocytes at 3 and 7dpi (Fig. 5b).

2) Vegfa binding to Vegfr1 and Vegfr2 on endothelial cells facilitates the proliferation, survival and directional sprouting of tip cell filipodia during angiogenesis17–19. The highest interaction scores for these ligand-receptor pairs were associated with monocytes/macrophages and astrocytes (Fig. 5d). 

3) Vegfa expression was highest in these two cell types at 1dpi, and Vegfr1 and Vegfr2 receptors were highly expressed by endothelial cells at all time points, suggesting that the major cues for new vessel formation is derived from the infiltrating myeloid cells and astrocytes (Fig. 5e). 

4) Strikingly, the strongest interactions amongst the Vegf family members were associated with placental growth factor (Plgf) binding to Vegfr1 (Fig. 5d).

<h3> Analysis of macroglia heterogeneity reveals astrocyte and OPC-specific roles during gliosis </h3>

1) To assess macroglia heterogeneity, oligodendrocyte lineage cells, astrocytes, and ependymal cells were clustered and visualized on a separate UMAP (Fig. 6a, b), which showed spatial segregation corresponding to their major cell type.

2) We identified three ependymal and an astroependymal subtypes with distinct expression and temporal profiles (Fig. 6c, e, Extended Data Fig 8). 

3) The astropendymal cells were best identified by their expression of the crystallin crym, and shared common markers with both astrocytes (e.g. timp1, gfap) and ependymal cells (e.g. vim, tmsb10).

4) The oligodendrocyte lineage cells segregated into previously described clusters that showed a spatial progression from OPC to dividing OPC and preoligodendrocytes, and finally to oligodendrocytes (Fig. 6a).

5) These populations were identified by prototypical markers, with the exception of OPCs, which were best identified by tenascin tnr (Fig. 6e, Extended Data Fig. 8b). 

6) Whereas mature oligodendrocytes were predominant in the uninjured spinal cord, the proportion of OPCs, dividing OPCs, and preoligodendrocytes gradually increased over the next 7dpi. 

7) To compare reactive OPCs and astrocytes, we performed Gene Ontology (GO) Enrichment Analysis for biological processes associated with DEGs between each time point (Fig. 6f).

8) At 1dpi, top biological processes for both astrocytes and OPCs pertain to translation and biogenesis. By 3dpi, astrocytes are defined by processes related to neurogenesis and gliogenesis, whereas OPCs are defined by mitosis, which reflect the active state of proliferation and differentiation for both cell types during this time.

9) To assess the potential effects of reactive astrocytes and OPCs on axonal growth, we compared the expression levels of axon growth inhibitory molecules using a previously curated list24 (Fig. 6g).

10) Interestingly, inhibitory proteoglycans such as acan, bcan, ncan, and vcan were expressed preferentially by OPCs.

<h3> Macrophage-Mediated Mechanisms of Gliosis and Fibrosis </h3>

1) Since IL6 cytokine family members are the main STAT3 activators, we assessed their expression across all cell types and found oncostatin M (osm) expressed at highest levels in myeloid cells, il6 expressed at highest levels in fibroblasts, and clcf1 expressed at highest levels in astroependymal cells (Fig. 7a).

2) Other IL6 cytokines were either not expressed highly or dropped out of our sequencing analysis.

3) Oncostatin M receptor (osmr) and the signaling coreceptor gp130 (il6st) were expressed highly in fibroblasts and astrocytes, and interaction scores for IL-6 cytokine family members were highest for OSM signaling between Ly6Clo monocytes/Macrophage-B and astrocytes/fibroblasts (Fig. 7b).

4)  To validate the sequencing data, we performed double in situ hybridization for osm and cd11b and found a significant increase in osm mRNA in cd11b+ myeloid cells compared to non-myeloid cells (Fig. 7c, e).

5) To assess OSMR expression, we used immunohistochemistry and found increasing expression of OSMR in both Pdgfr-β+ fibroblasts and GFAP+ astrocytes in the fibrotic and glial scar, respectively (Fig. 7d, f, g).

6)  Taken together, our results strongly suggest that the gp130 signaling pathway induced by Osm is a common mechanism by which astrocytes and fibroblasts are preferentially activated by specific monocyte/macrophage subtypes after SCI.

7) To identify other distinct and common pathways by which macrophage subtypes mediate astrogliosis and fibrosis, we calculated interaction scores for all known ligand-receptor pairs between Macrophage-A/B and astrocytes/fibroblasts at 3 and 7dpi (Fig. 8).

8) While many ligands expressed by macrophages signaled to receptors on both astrocytes and fibroblasts (i.e. common pathways), there were many more ligand-receptor pairs unique to macrophage-fibroblast interactions than macrophage-astrocyte interactions (i.e. distinct pathways). These unique macrophage-fibroblast interactions included signaling related to IL1α/β, Vegfa/b, Pdgfα, and Tgfβ1.

9)  Overall, the highest interactions scores were associated with Spp1 and Apoe signaling, which were common to both astrocytes and fibroblasts.

10) Macrophages also displayed subtype specificity in signaling to astrocytes and fibroblasts. For example, Jag1-Notch signaling was largely specific to Macrophage-A-astrocyte interaction, whereas IL1α/β-IL1r1 signaling was specific to Macrophage-B-fibroblast interaction. In summary, our analysis highlights the utility of our sc-RNAseq dataset in identifying potential signaling mechanisms that mediate astrogliosis and fibrosis after SCI.


<h3> Discussion </h3>

1) We found that myeloid subtypes display distinct temporal regulation of angiogenesis, gliosis, and fibrosis. These results support previous studies showing that myeloid depletion leads to reduced angiogenesis, astrogliosis, and fibrosis after CNS injury3, 27, 28, and provide further insight by identifying potential contributions of specific myeloid subtypes.

<h4> 2) Our finding that Macrophage-A and B subtypes do not correspond to the M1/M2 nomenclature (Extended Data Fig. 5) is consistent with findings from other sc-RNAseq studies29–31. </h4>

3) The vascular analysis revealed novel insight in the contribution of tip cells and astrocytes during angiogenesis after SCI. The data indicate that tip cells are highly dynamic; they appear quickly at 1dpi and are mostly gone by 7pi.

4) Our single cell RNA-seq dataset is the first comprehensive transcriptomic analysis that captures virtually all cells that contribute to the injury site pathology after SCI. This dataset can be used to assess not only heterogeneity of the cells that comprise the injury site, but also to assess signaling mechanisms that underlie cellular interactions at the injury site.


# Methods Section

<h3> Mice & Spinal Cord Injuries </h3
  
Wet lab, TODO: review

<h3> Tissue Dissociation </h3>

Wet lab, TODO: review

<h3> Histology </h3>

Wet lab, TODO: review

<h3> Immunohistochemistry </h3> 

Wet lab, TODO: review

<h3> In situ hybridization </h3>

Wet lab, TODO: review

<h3> Tissue Quantifications </h3>

Wet lab, TODO: review

<h3> Flow Cytometry </h3>

Wet lab, TODO: review

<h3> Single cell RNA-sequencing using 10X Genomics platform </h3>

1) Two biological replicates for each time point for a total of 8 samples were sequenced with a median of ~8000 cells per sample.

2) The median of mean reads per cell across samples was > 55,000 and median sequencing saturation was 81.2%.

3) Libraries for all samples were prepared according to Chromium Single Cell 3’ Library and Gel Bead Kit v2 instructions (10X Genomics, PN-120237) and indexed with Chromium i7 Multiplex Kits (10X Genomics, PN-120262).

4) Primers contain i) a 16nt 10X barcode, ii) a 10nt Unique Molecular Identifier, iii) a poly-dT primer sequence, and iv) an Illumina Read 1 (R1) sequence to produce single-stranded, barcoded complementary DNA (cDNA) from poly-adenylated mRNA. 

5) Reads were then filtered to remove UMIs and barcodes with single base substitution errors and finally used for UMI counting. The output was a count matrix containing all UMI counts for every droplet. CellRanger v2.2.0 was used for the first set of replicates and v3.0.1 for the second set of replicates.


<h3> Preprocessing and Quality Control </h3>

1) First, cell barcodes were ranked according to UMI count and visualized in a log-total UMI count vs log-rank plot. A spline curve was fit to the data to identify knee and inflection points. At least all data points above the knee were considered cell-containing droplets.

2) In order to further distinguish cells from data below the knee, we used the emptyDrops function from the DropletUtils R package52 using the following fixed parameters for all sample: lower = 250; max fit.bounds = 1e06; FDR = 0.001; ignore = 10.

3) Some parameters varied between samples accordingly: “retain” was set to knee point values, and “lower” was set to inflection point values. The result was a filtered count matrix of all putative cell-containing droplets. In order to distinguish low quality cells, we considered cell-level metrics such as library size, percentage of UMIs mapping to mitochondrial genes, and doublet detection algorithm outputs. 

4) We observed that removing cells in the lowest quantiles of total UMI counts preferentially selected against endothelial cells, while removing cells in the highest quantiles of mitochondrial percentage preferentially selected against astrocytes.

5) To avoid bias due to any one metric, we used a multivariate approach to automatic outlier detection as implemented in the Scater R package45, which performs principal component analysis on the quality control metrics to determine outliers.

6) In order to remove potential doublets, we applied the Scrublet Python package47 for each individual sample. Cells from the data that have high local densities of simulated doublets are flagged and removed.

7) We set the expected_doublet_rate for each sample according to the estimated doublet rate per cells sequenced as published by 10X, and default values for all other parameters.

8) During downstream analysis, we observed small, unidentifiable clusters that co-expressed marker genes for two different cell-types, and removed these cells as they were likely to be multiplets not detected by Scrublet. 

9) Additionally, we identified a cluster enriched for hemoglobin genes (Hbb-bs, Hbb-bt, etc.) which was removed from downstream analyses.

<h3> Normalization and Batch Correction </h3>

1) To generate the full SCI dataset, all samples were processed and combined using Seurat v353.

2) After filtering each sample count matrix for genes that were expressed in at least 10 cells, each dataset was independently normalized and scaled using the SCTransform function.

3) To remove cell-cycle genes as a confounding source of variation, mitochondrial percentage and cell cycle scores based on the expression of canonical G2M and S phase markers were computed for each cell. Cell cycle genes were provided through the Seurat tutorial.

4) These score values were then used as input for the “vars.to.regress” argument in the SCTransfrom() function.

5) To identify shared and unique molecular cell-types across datasets and time-points, sample expression matrices were batch-corrected using Seurat’s Data Integration workflow.

6) For the full SCI dataset, the 2000 most variables genes were used as input for the “anchor. features” argument of the FindIntegrationAnchors() function, where the variance of a gene was measured as the residual divided by the expected variance under the SCTransfrom() model. This resulted in a single, batch-corrected expression matrix for containing all cells.

7) For the analysis of the myeloid cells, we tested a similar batch correction as described above. However, we observed that some microglia from the uninjured spinal cord were classified under a peripheral macrophage subset but not under BA-Macrophage (data not shown). 

8) We instead (potential misclassification) combined sample datasets across time-points for each replicate and subsequently performed normalization for each of the two larger datasets as above (i.e. normalize across times and batch-correct across replicates). For myeloid integration, the 2000 most variable genes were used. For the vascular cells, no obvious misclassifications were noted. 

9) We proceeded with normalization and batch correction as above, using the 2000 most variable genes for downstream analysis.

10) For the macroglia cells, we again observe no obvious misclassifications and proceeded with normalization and batch correction. However because of the low numbers of macroglia cells in one of the uninjured sample replicates, we adjusted the “k.filter” parameter to 100 for the FindIntegrationAnchors() function.

<h3> Dimensional Reduction, Clustering, and Differential Gene expression testing </h3>









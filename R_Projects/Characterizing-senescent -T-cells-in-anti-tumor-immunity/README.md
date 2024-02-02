# Characterizing Senescent T cells in Anti-Tumor-Immunity

This project is looking at how T cells found within the B16F10 mouse melanoma model, a melanoma that posses high metastatic capacity, are characterized by a senescent signature.
This genotype and phenotype is different from exhaustion and maybe a subset of interest for reversing T cell dysfunction.

Transcriptomic data (Clariom assay) was performed by a collaborator, pre-processing was completed using thermofisher's proprietary software for the assay. 
I generated differentially expressed genes and visualized them through volcano plots, a heatmap and performed principal component analysis to show distinct populations of T cells given their various stimuli (B16 cells, MC38 cells, or LCMV).
I then used gene ontology to annotate pathways associated with the differentially expressed genes for improved interpretation of transcriptomic changes in the T-cells. 
** This data is within this repository as the manuscript is still being prepared, however once accepted for publication the data will be available.**

To better understand if this T cell phenotype, I used publically deposited single cell RNA-seq data 
([GSE173351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173351) and 
[GSE185206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185206))
from patients who responded and didn't respond to PD-1 checkpoint blockade. Here I hypothesized that those whom didn't respond would be enriched for the same differentially expressed genes identified in the B16 model.
While those who did respond would lack the genes identified from the T cells in the B16 model. 


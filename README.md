# CompBio
Repository of Computational Biology Side projects

## About Me
I am a Ph.D Candidate in Immunology and Cancer Biology. I'm passionate about improving our understanding of onco-immunology through next generation computational methods. 
I am constantly learning and always looking to learn a new coding language and improve. 

## Contact Information
Jordon Inkol, Ontario Canada
- [Email](mailto:jinkol@uoguelph.com?subject=[GitHub]%20Source%20Han%20Sans)
- [LinkedIn](https://www.linkedin.com/in/jordon-inkol-145932257)
- [Kaggle](https://www.kaggle.com/jordoninkol)

## Projects completed in Python

### Analysis of Human Head and Neck Squamous Cell Carcinoma Immune cells

 This analysis was performed on data originally from [**Zheng et al.** published in Nature Communications in 2020](https://www.nature.com/articles/s41467-020-20019-0). In this study they examined the immune landscape of 7 patients with head and neck squamous cell carcinoma and compared to matched normal samples from each respective patient. 

  For this project I sought to perform my own scRNA-seq analysis in Python using only scverse libraries to complete the pre-processing and analysis/visualization. These included scanpy and scvi-tools to perform normalization and identification of novel cell type markers. I further expanded on this analyzing potential novel ligand-receptor interactions between Dendritic cells and T cells as well as Macrophages and T cells by employing CellPhoneDB from the [Teichmann lab](http://www.teichlab.org/) at the Sanger Institute. 

  Unfortunately the h5ad files generated during this project were too large to upload, however I am more than happy to provide them upon request via email. 


## Projects completed in R
### [Gasdermin E Expression in Human Tumors](R_Projects/Gasdermin_E_Expresion_in_Human_Tumors)
This is a small R project looking at copy number variation and RNA expression of Gasdermin E in human breast, ovarian and glioblastoma tumors.
Expression of RNA transcripts was also stratified by grade and stage, where applicable, to identify any changes in expression that correlate with increased disease severity.

All data is from The Cancer Genome Atlas

### [Expression of Mouse Colon Adenocarcinoma Neoantigens](R_Projects/Mouse-Tumor-Neoantigen)

I preformed RNA-sequencing on MC-38 mouse colon adenocarcinoma cells with or without expression of neoantigen or model antigens (AdFL, Ad9mer, and gp33). This dataset was used to quantify and identify expression of other neoantigens and tumor associated antigens found within these cells.

I preprocessed the raw data using BWA, Samtools, GATK, and Picard through interfacing with an HPC(Graham Cluster:SHARCNet/Digital Research Alliance of Canada) via bash. Downstream analysis and visualization was completed using EdgeR and ggplot2 in R. This data is currently in a publication under review at Cancer Immunology Research (CIR-23-0639R)

### [Characterizing Senescent T cells in Anti-tumor Immunity](R_Projects/Characterizing-senescent-T-cells-in-anti-tumor-immunity)
This project is looking at how T cells found within the B16F10 mouse melanoma model, a melanoma that posses high metastatic capacity, are characterized by a senescent signature. This genotype and phenotype is different from exhaustion and maybe a subset of interest for reversing T cell dysfunction.

Transcriptomic data (Clariom assay) was performed by a collaborator, pre-processing was completed using thermofisher's proprietary software for the assay. I generated differentially expressed genes and visualized them through volcano plots, a heatmap and performed principal component analysis to show distinct populations of T cells given their various stimuli (B16 cells, MC38 cells, or LCMV). I then used gene ontology to annotate pathways associated with the differentially expressed genes for improved interpretation of transcriptomic changes in the T-cells. ** This data is within this repository as the manuscript is still being prepared, however once accepted for publication the data will be available.**

To better understand if this T cell phenotype, I used publically deposited single cell RNA-seq data (GSE173351 and GSE185206) from patients who responded and didn't respond to PD-1 checkpoint blockade. Here I hypothesized that those whom didn't respond would be enriched for the same differentially expressed genes identified in the B16 model. While those who did respond would lack the genes identified from the T cells in the B16 model.



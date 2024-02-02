#Load in the necessary packages to complete the analysis;
#Seurat is the main scRNA-seq analysis library
#ggplot2, dyplr,and pathwork all help with visualizations

library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(patchwork)

#Define the functions for downstream use
do.QC <- function(obj){
  
  metadata <- obj@meta.data
  metadata$cells <- rownames(metadata)
  metadata <- metadata %>%
    dplyr::rename(sample = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  obj@meta.data <- metadata
  
  plot1 <- ggplot(metadata, aes(x=sample, fill = sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  plot2 <- ggplot(metadata, aes(color=sample, x=nUMI, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  plot3 <- ggplot(metadata, aes(color=sample, x=nGene, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300) 
  
  return(list(plot1,plot2,plot3))
  
}
do.filter <- function(sample){
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
  sample <- subset(sample, subset = nFeature_RNA >200 & nFeature_RNA < 2500 & percent.mt < 20)
  sample <- SCTransform(sample, vst.flavor = "v2", verbose = FALSE)
  return(sample)
}
do.RPCA <- function(sample){
  sample <- NormalizeData(sample, verbose = FALSE)
  sample <- FindVariableFeatures(sample, verbose = FALSE)
}
do.DR <- function(sample){
  sample <- ScaleData(sample)  
  sample <- RunPCA(sample)
  sample <- RunUMAP(sample, dims = 1:30, verbose = FALSE)
  sample <- FindNeighbors(sample, dims = 1:30, verbose = FALSE)
  sample <- FindClusters(sample)
}

#Set working directory to parent directory
setwd("./SCRNA-Seq/non_responder/")


#This for loop loads in files (MD01-LN...) in the current working directory
for (file in c("MD01-LN_1", "MD01-LN_2", "MD01-LN_3")){
  #Here we are loading in the 10X scRNA-seq data post Cell ranger processing
  seurat_data <- Read10X(data.dir = paste0("./", file))
  #Create a seurat object for each min previously loaded, the criteria are each gene must be found in at least 3 cells and each cell must have atleast 200 genes
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 200, project = file)
  #name each seurat object after its respective file name
  assign(file, seurat_obj)
}

nonresponder.combined <- merge(`MD01-LN_1`, y = c(`MD01-LN_2`, `MD01-LN_3`), add.cell.ids = nonresponder, project = nonresponder)

nonresponder.QC <- do.QC(nonresponder.combined)

plot(nonresponder.QC[[1]])
plot(nonresponder.QC[[2]])
plot(nonresponder.QC[[3]])

nonresponder.list <- list(`MD01-LN_1`, `MD01-LN_2`, `MD01-LN_3`)

nonresponder.final <- merge(nonresp.list[[1]], y = c(nonresp.list[[2]],nonresp.list[[3]]), add.cell.ids = nonresponder, project = nonresponder)

nonresponder.final <- do.filter(nonresponder.final)

nonresponder.final <- lapply(nonresponder.final, do.filter)

nonresponder.final <- lapply(nonresponder.final, do.RPCA)

nonresponder.final <- lapply(nonresponder.final, do.DR)


#Set the working directory for the responders
setwd("./SCRNA-Seq/responder/")

#Re-iterate through the responder patients using the same pipeline as with the non-responders
for (file in c("msk1263.gene", "msk1302.gene", "msk1263.tm.r5")){
  #Here we are loading in the 10X scRNA-seq data post Cell ranger processing
  seurat_data <- Read10X(data.dir = paste0("./", file))
  #Create a seurat object for each min previously loaded, the criteria are each gene must be found in at least 3 cells and each cell must have atleast 200 genes
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 200, project = file)
  #name each seurat object after its respective file name
  assign(file, seurat_obj)
}

responder.combined <- merge(`msk1263.gene`, y = c(`msk1302.gene`, `msk1263.tm.r5`), add.cell.ids = responder, project = responder)

responder.QC <- do.QC(responder.combined)

plot(responder.QC[[1]])
plot(responder.QC[[2]])
plot(responder.QC[[3]])

responder.list <- list(`msk1263.gene`, `msk1302.gene`, `msk1263.tm.r5`)

responder.final <- merge(responder.list[[1]], y = c(responder.list[[2]],resp.list[[3]]), add.cell.ids = responder, project = responder)

responder.final <- do.filter(responder.final)

responder.final <- lapply(responder.final, do.filter)

responder.final <- lapply(responder.final, do.RPCA)

responder.final <- lapply(responder.final, do.DR)



#Visualize the clusters via Dim plot functionin Seurat

DimPlot(nonresponder.list[1], label = TRUE)
DimPlot(nonresponder.list[2] = TRUE)
DimPlot(nonresponder.list[3] = TRUE)

DimPlot(responder.list[1], label = TRUE)
DimPlot(responder.list[2], label = TRUE)
DimPlot(responder.list[3], label = TRUE)



#Visualize the T cell senescence genes identified in previous figures/data from T cell exposed to different tumor model stimuli

FeaturePlot(nonresponder.list[1], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))
FeaturePlot(nonresponder.list[2], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))
FeaturePlot(nonresponder.list[3], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))


FeaturePlot(responder.list[1], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))
FeaturePlot(responder.list[2], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))
FeaturePlot(responder.list[3], features = c("CD8A", "CDKN1A","H2AFX", "CD27", "CD28"))







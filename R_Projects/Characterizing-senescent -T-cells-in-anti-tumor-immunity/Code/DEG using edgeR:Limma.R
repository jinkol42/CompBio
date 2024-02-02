#Load in the required libraries
library(edgeR)
library(dplyr)

#Here the data from the output of the Clariom assay is loaded

B16 <- read.delim("~/Desktop/Senescence/Expression Data (ClariomS) natural 2.csv", comment.char="#")

#In this experiment Central Memory T cells were included in the raw data however they are not of interest so they were filtered out. 

deg <- select(B16, -c(TCM1:TCM.M3))
deg <- round(deg)

# Here we use EdgeR to normalize the counts, assign group variables for comparison and then fit a glmQL model. 
d0 <- DGEList(counts = deg, group = group)
d0 <- calcNormFactors(d0)

group <- factor(c(1,1,1,2,2,2,2,3,3,3,3))

design <- model.matrix(~0+group, data = d0$samples)

colnames(design) <- levels(d0$samples$group)

d0 <- estimateDisp(d0, design)

fit <- glmQLFit(d0, design)

#Here we set the comparrisons to find unique transcripts of T cells exposed to different stimuli (MC38 cells, B16 cells, or LCMV )
qlfB16 <- glmQLFTest(fit, contrast = c(-0.5,1,-0.5))
qlfMC38 <- glmQLFTest(fit, contrast = c(-0.5,-0.5,1))
qlfLCMV <- glmQLFTest(fit, contrast = c(1,-0.5,-0.5))

#Save the outputs from EdgeR and export them as csv files.
B16vall <- as.data.frame(qlfB16)
MC38vall <- as.data.frame(qlfMC38)
LCMVvall <- as.data.frame(qlfLCMV)

write.csv(B16vall, "\\B16vall.csv")
write.csv(MC38vall, "\\MC38vall.csv")
write.csv(LCMVvall, "\\LCMVvall.csv")

#Here we filter out any genes that have less than a 1 log-fold change and a p-value less than 0.01
Mc38v1 <- filter(MC38vall, logFC >= abs(1) & PValue < 0.01)
B16v1 <- filter(B16vall, logFC >= abs(1) & PValue < 0.01)
LCMVv1 <- filter(LCMVvall, logFC >= abs(1) & PValue < 0.01)


#Save the filtered datasets as csv files for analysis downstream using Gene Ontology
write.csv(Mc38v1, "\\Mc38GO.csv")
write.csv(B16v1, "\\B16GO.csv")
write.csv(LCMVv1, "\\LCMVGO.csv")


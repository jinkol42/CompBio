#Load in required libraries/packages
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ggtext)
library(ggsci)
library(PCAtools)

#load in the required dataset and prep that data frame for downstream analysis

Expression.Data <- read.delim("~/Desktop/Senescence/Expression Data (ClariomS) natural 2.csv", comment.char="#")

B16 <- select(Expression.Data, -c(2))
B16 <- B16[-1,]
colnames(B16) <- B16[1,]
B16 <- B16[-1,]

#Create a list of stimuli names to create a column
cond <- factor(c("LCMV", "LCMV", "B16", "B16", "B16", "B16", "MC38", "MC38", "MC38", "MC38"))


# Select the columns to generate PCA and calculate loadings
B16.pca <- select(B16, LCMV2:MC38.4)

B16.pca.round <- as.numeric(B16.pca)
B16.pca.round <- round(B16.pca)

B16.pca$LCMV2 = as.numeric(B16.pca$LCMV2)

coldata <- data.frame(row.names = colnames(B16.pca), cond)

dds <- DESeqDataSetFromMatrix(countData = round(B16.pca), colData = coldata, design = ~cond)

vts <- vst(dds, blind = FALSE, fitType = c("mean"))

pcafig <- plotPCA(vts, intgroup = "cond", returnData= TRUE)

datapca <- as.data.frame(pcafig)

datapca$cond <- factor(datapca$cond, levels = c("LCMV", "B16", "MC38"))

#Visualize PCA using ggplot2
pca <- ggplot(datapca, aes(x=PC1, y=PC2, color = cond))+ 
  geom_point(size = 8)+ 
  geom_vline(xintercept = c(16.5, 3.5), linetype="dashed", color ="gray60" )+
  geom_hline(yintercept = c(22, -1), linetype="dashed", color ="gray60")+
  labs(x = "PC1 (66.4%)", y = "PC2 (20.5%)", title = "PCA")+
  theme(panel.border = element_rect(color = "gray20", fill = NA, size = 1),panel.background = element_rect(fill = 'transparent'), 
        plot.background = element_rect(fill = 'transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = 'transparent'), legend.key = element_rect(fill = 'transparent'), legend.title = element_text(size = "30"),
        legend.text = element_markdown(size = "30"), axis.text = element_text(size = "30", color = "gray20"), axis.title = element_text(size = "30"), axis.title.x = element_text(vjust = -1.5),
        plot.margin = margin(25,25,15,25), plot.title = element_text(size = "30"))+
  scale_colour_npg(name="Groups", labels=c("T<sub>CM-LCMV</sub>","T<sub>TP-B16</sub>", "T<sub>TP-MC38</sub>"))





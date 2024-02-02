#Load the required libraries
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

#Load in the raw data and filter based upon at least 1 log fold change and a P value less than 0.01
p <- filter(pheat_data, B16.logFC >= abs(1) & `B.16 PValue` < 0.01 | LCMV.logFC >= abs(1) & LCMV.PValue < 0.01 | MC38.logFC >= abs(1) & MC38.PValue < 0.01)
p.log <- select(p, NAME, B16.logFC, MC38.logFC, LCMV.logFC)
colnames(p.log) <- c("Gene", "B16", "MC38", "LCMV")


#Here we generate a matrix based upon all the genes, the top 300 changed genes, or the top 100 changed genes. 
t300p.log <- head(p.log,200)
t100.log <- rbind(head(p.log,50), tail(p.log,50))

t300.mat <- as.matrix(t300p.log[,-1])
p.mat <- as.matrix(p.log[,-1])
t100.mat <- as.matrix(t100.log[,-1])

#change the row names to based upon the gene name column
rownames(t300.mat) <- t300p.log$Gene
rownames(p.mat) <- p.log$Gene
rownames(t100.mat) <- t100.log$Gene

#Create the matrix required for plotting
base_mean = rowMeans(p.mat)
p.mat_scaled = t(scale(t(p.mat)))
type = gsub("s//d+_", "", colnames(p.mat))

base_mean = rowMeans(t300.mat)
t300.mat_scaled = t(scale(t(t300.mat)))
type = gsub("s//d+_", "", colnames(t300.mat))

base_mean = rowMeans(t100.mat)
t100.mat_scaled = t(scale(t(t100.mat)))
type = gsub("s//d+_", "", colnames(t100.mat))

#Annotate the genes of interest/ particular genes relative to T cell senescence
target_genes <- c("Zap70", "Cd8a", "Cd3e", "Cd27", "Gzmk", "Cd69", "Eomes", "Irf7", "Cdkn1a")
r <- which(p.log$Gene %in% target_genes)

## All genes included

#Create a dummy variable that can store our genes of interest to be called later in the plotting function
r <- which(p.log$Gene %in% target_genes)

#Defined color range for gene expression scale
col_fun = colorRamp2(c(-2,0,2), c("#3C5488B2", "white","#DC0000B2"))

#This annotation is used by pheatmap to define the legend as well as colors to coordinate to the 'stimulus'
ha = HeatmapAnnotation(Samples = type,
                       col = list(Samples = c("LCMV" = "#1b9e77", "MC38" = "#7570b3", "B16" = "#d95f02")),
                       show_legend = FALSE,
                       border = TRUE)

lgd = Legend(title = "Samples",
             labels = gt_render(c("T<sub>TP-B16</sub>", "T<sub>CM-LCMV</sub>", "T<sub>TP-MC38</sub>")),
             legend_gp = gpar(fill = c("#d95f02", "#1b9e77", "#7570b3" )),
             border = TRUE) 

#These are the row positions of the our target genes/ genes of interest that are being stored in the dummy variable 'r'
tg = rowAnnotation(foo = anno_mark(at = c(1,358,405,408,583,599,1039),
                                   labels = target_genes[1:10],
))
#Plot the actual heatmap using pheatmap
allgeneshtmp <- Heatmap(p.mat_scaled, column_dend_side = "top", 
                        clustering_distance_rows = "pearson", 
                        col = col_fun,
                        column_order = c("LCMV", "MC38", "B16"),
                        top_annotation = ha,
                        show_column_names = FALSE,
                        width = unit(6,"cm"),
                        border_gp = gpar(col="black"),
                        rect_gp = gpar(col = "black", lwd=0.03),
                        show_row_names = FALSE,
                        heatmap_legend_param = list(
                          title = "", at = c(-2, -1, 0, 1, 2),
                          legend_height = unit(5, "cm"),
                          border = "transparent"
                        )
)
#This renders the heatmap
allgenesfinal.htmp  <- draw(allgeneshtmp + tg, padding = unit(c(2,2,2,2), "mm"), heatmap_legend_list = list(lgd)) 

##Top 500 genes selected only
col_fun = colorRamp2(c(-2,0,2), c("#3C5488B2", "white","#DC0000B2"))

ha = HeatmapAnnotation(Samples = type,
                       col = list(Samples = c("LCMV" = "#1b9e77", "MC38" = "#7570b3", "B16" = "#d95f02")),
                       show_legend = FALSE,
                       border = TRUE)

lgd = Legend(title = "Samples",
             labels = gt_render(c("T<sub>TP-B16</sub>", "T<sub>CM-LCMV</sub>", "T<sub>TP-MC38</sub>")),
             legend_gp = gpar(fill = c("#d95f02", "#1b9e77", "#7570b3" )),
             border = TRUE) 

tg = rowAnnotation(foo = anno_mark(at = c(1,358,599,405,900,408,583,482, 1039),
                                   labels = target_genes[1:10],
))

t300htmp <- Heatmap(t300.mat_scaled, column_dend_side = "top", 
                    clustering_distance_rows = "pearson", 
                    col = col_fun,
                    column_order = c("LCMV", "MC38", "B16"),
                    top_annotation = ha,
                    show_column_names = FALSE,
                    width = unit(6,"cm"),
                    border_gp = gpar(col="black"),
                    rect_gp = gpar(col = "black", lwd=0.03),
                    row_names_gp = gpar (fontsize = "4"),
                    heatmap_legend_param = list(
                      title = "", at = c(-2, -1, 0, 1, 2),
                      legend_height = unit(5, "cm"),
                      border = "transparent"
                    )
)

t300final.htmp  <- draw(t300htmp, padding = unit(c(2,2,2,2), "mm"), heatmap_legend_list = list(lgd))


#random sample of 100 genes
col_fun = colorRamp2(c(-2,0,2), c("#3C5488B2", "white","#DC0000B2"))

ha = HeatmapAnnotation(Samples = type,
                       col = list(Samples = c("LCMV" = "#1b9e77", "MC38" = "#7570b3", "B16" = "#d95f02")),
                       show_legend = FALSE,
                       border = TRUE)

lgd = Legend(title = "Samples",
             labels = gt_render(c("T<sub>TP-B16</sub>", "T<sub>CM-LCMV</sub>", "T<sub>TP-MC38</sub>")),
             legend_gp = gpar(fill = c("#d95f02", "#1b9e77", "#7570b3" )),
             border = TRUE) 

tg = rowAnnotation(foo = anno_mark(at = c(1,358,599,405,900,408,583,482, 1039),
                                   labels = target_genes[1:10],
))

t100htmp <- Heatmap(t100.mat_scaled, column_dend_side = "top", 
                    clustering_distance_rows = "pearson", 
                    col = col_fun,
                    column_order = c("LCMV", "MC38", "B16"),
                    top_annotation = ha,
                    show_column_names = FALSE,
                    width = unit(6,"cm"),
                    border_gp = gpar(col="black"),
                    rect_gp = gpar(col = "black", lwd=0.03),
                    row_names_gp = gpar (fontsize = "7"),
                    heatmap_legend_param = list(
                      title = "", at = c(-2, -1, 0, 1, 2),
                      legend_height = unit(5, "cm"),
                      border = "transparent"
                    )
)

t100final.htmp  <- draw(t100htmp, padding = unit(c(2,2,2,2), "mm"), heatmap_legend_list = list(lgd))



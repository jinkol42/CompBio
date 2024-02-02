#Load in the required packages and libraries
library(ggplot2)
library(ggtext)
library(dplyr)
library(ggrepel)

#Load in the required data sets

Expression.Data <- read.csv("~/Desktop/Senescence/B16vLCMV.csv")
Expression.Data <- read.csv("~/Desktop/Senescence/MC38vB16.csv")


#Remove na values from B16vLCMV dataframe
B16vLCMV <- na.omit(B16vLCMV)

B16vLCMV %>% # Identify the xlim values
  pull(Log_Change) %>%
  min() %>%
  floor() 

B16vLCMV %>% # Identify the ylim values
  pull(Log_Change) %>%
  max() %>%
  ceiling() 

max(abs(-5), abs(9))

# label genes up, down or ns from DEGS list
B16vLCMV <- B16vLCMV %>%
  mutate(gene_type = case_when(Log_Change >= 1 & `FDR_P-val` <= 0.01 ~ "up",
                               Log_Change <= -1 & `FDR_P-val` <= 0.01 ~ "down",
                               TRUE ~ "ns"))

#Dictate colors, point size and alpha depending on DEG level
cols <- c("up" = "#DC0000B2", "down" = "#3C5488B2", "ns" = "grey")
sizes <- c("up" = 1.5, "down" = 1.5, "ns" = 1)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

sig.genes <- rbind(head(B16vLCMV, n=5), tail(B16vLCMV, n=5))
sig.genes <- as.data.frame(sig.genes)

top.up <- head(sig.genes, n=5)
top.down <- tail(sig.genes, n=5)



#Volcano plot for B16 vs LCMV T cells
vol.plot <-ggplot(B16vLCMV, aes(x = Log_Change,
             y= -log10(`FDR_P-val`),
             fill= gene_type,
             size = gene_type,
             alpha = gene_type))+
  labs(
    x = "Log<sub>2</sub>(Fold Change)",
    y = "-Log<sub>10</sub>(FDR *p* value)",
    title = "T<sub>TP-B16</sub> vs. T<sub>CM-LCMV</sub>"
  )+
  geom_point(colour = "black",
             shape = 21)+
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed")+
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed")+
  geom_text_repel(data = subset(sig_genes_B16vLCMV, Log_Change > 0),
                   aes(label = Gene_Symbol, size = "8",
                       segment.curvature = -0.1,
                       segment.ncp = 5,
                       segment.angle = 30),
                  nudge_x = 7 - subset(sig_genes_B16vLCMV, Log_Change >0)$Log_Change,
                  nudge_y = 0.5,
                  hjust = 0,
                  direction = "y",
                  min.segment.length = 0
                   )+
  geom_text_repel(data = subset(sig_genes_B16vLCMV, Log_Change < 0),
                  aes(label = Gene_Symbol, size = "8",
                      segment.curvature = -0.1,
                      segment.ncp = 5,
                      segment.angle = 30),
                  nudge_x = -5 - subset(sig_genes_B16vLCMV, Log_Change <0)$Log_Change,
                  hjust = 1,
                  direction = "y"
  )+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,12))+
  scale_fill_manual(values = cols)+
  scale_size_manual(values = sizes)+
  scale_alpha_manual(values = alphas)+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(25,25,15,20),
        legend.position = "none",
        plot.title = element_markdown(size = "16", margin = margin(0,0,20,0)),
        plot.title.position = 'plot',
        axis.title.x = element_markdown(margin = margin(10,1,1,1), size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(color = "black", size = "10"))+
  coord_cartesian(expand = FALSE, clip = 'off')


# Rename columns for ease of calling

colnames(MC38vB16)[5] = "Log_Change"
colnames(sig.total.mc38)[5] = "Log_Change"

# label genes up, down or ns from DEGS list, 1 Log2 fold change is the same as 2 fold change
MC38vB16 <- MC38vB16 %>%
  mutate(gene_type = case_when(Log.Change >= 1 & `FDR_P-val` <= 0.01 ~ "up",
                               Log.Change <= -1 & `FDR_P-val` <= 0.01 ~ "down",
                               TRUE ~ "ns"))

#Create a lists of the significantly up and downregulated genes for plotting
sig.up.mc38 <- MC38vB16 %>%
  filter(Gene_Symbol %in% c("Fam132a", "Ermap", "Add2", "Atp1b2", "Tfr2"))

sig.down.mc38 <- MC38vB16 %>%
  filter(Gene_Symbol %in% c("Crmp1", "Ccr4", "Cd5", "Pdcd1", "Izumo1r"))

sig.total.mc38 <- rbind(sig.up.mc38, sig.down.mc38)

# Plot B16 vs MC38 Volcano
vol.mc38 <-ggplot(MC38vB16, aes(x = Log_Change,
                                y= -log10(`FDR_P-val`),
                                fill= gene_type,
                                size = gene_type,
                                alpha = gene_type))+
  labs(
    x = "Log<sub>2</sub>(Fold Change)",
    y = "-Log<sub>10</sub>(FDR *p* value)",
    title = "T<sub>TP-B16</sub> vs. T<sub>TP-MC38</sub>"
  )+
  geom_point(colour = "black",
             shape = 21)+
  geom_hline(yintercept = -log10(0.01),
             linetype = "dashed")+
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed")+
  geom_text_repel(data = subset(sig.total.mc38, Log_Change > 0),
                  aes(label = Gene_Symbol, size = "8",
                      segment.curvature = -0.1,
                      segment.ncp = 5,
                      segment.angle = 30),
                  nudge_x = 7 - subset(sig.total.mc38, Log_Change >0)$Log_Change,
                  nudge_y = 0.5,
                  hjust = 0,
                  direction = "y",
                  min.segment.length = 0
  )+
  geom_text_repel(data = subset(sig.total.mc38, Log_Change < 0),
                  aes(label = Gene_Symbol, size = "8",
                      segment.curvature = -0.1,
                      segment.ncp = 5,
                      segment.angle = 30),
                  nudge_x = -5 - subset(sig.total.mc38, Log_Change <0)$Log_Change,
                  hjust = 1,
                  direction = "y"
  )+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,12))+
  scale_fill_manual(values = cols)+
  scale_size_manual(values = sizes)+
  scale_alpha_manual(values = alphas)+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(25,25,15,20),
        legend.position = "none",
        plot.title = element_markdown(size = "16", margin = margin(0,0,20,0)),
        plot.title.position = 'plot',
        axis.title.x = element_markdown(margin = margin(10,1,1,1), size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.text = element_text(color = "black", size = "10"))+
  coord_cartesian(expand = FALSE, clip = 'off')

                

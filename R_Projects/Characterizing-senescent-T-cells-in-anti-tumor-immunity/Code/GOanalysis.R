#Load in the required packages and libraries

library(ggtext)
library(dplyr)
library(ggplot2)

#Load in the required datasets

B16GO <- read.csv("~/Desktop/Scott Senescence/GO_Analysis.csv")

#Visualization of GO terms 

label <- c(B16 = "T<sub>TP-B16</sub>", MC38 = "T<sub>TP-MC38</sub>", LCMV = "T<sub>CM-LCMV</sub>")

bubble <- ggplot(GO_analysis, aes(x= `Cell Line`, y = `Pathway`, size = `Number of genes`, color = `LogP`))+
  geom_point(data = GO_analysis)+
  scale_size(range = c(2,12), name = "Number of Genes")+
  scale_y_discrete(limits = rev(levels(GO_analysis$Pathway)), labels = function(x) str_wrap(x, width = 40))+
  scale_x_discrete(limits = c("B16", "MC38", "LCMV"), labels = label)+
  scale_color_gradient(low = "#a50f15", high = "#fcae91")+
  labs(
    y= "",
    x=""
  )+
  theme(panel.border = element_rect(color = "gray20", fill = NA, size = 1),panel.background = element_rect(fill = 'transparent'), 
              plot.background = element_rect(fill = 'transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.background = element_rect(fill = 'transparent'), legend.key = element_rect(fill = 'transparent'), legend.title = element_text(size = "12"),
              legend.text = element_text(size = "10", family = "Arial"), axis.text.x = element_markdown(size = "12", family = "Arial", color = "gray20", angle = 45, vjust = 1, hjust = 1), axis.title = element_text(size = "12"), axis.title.x = element_markdown(vjust = -1.5),
        axis.text.y = element_text(color = "gray20", size = "12"),
              plot.margin = margin(25,25,15,25))

bubble

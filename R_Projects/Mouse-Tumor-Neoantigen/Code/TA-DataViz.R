#Load in the required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(colorspace)
library(ggbeeswarm)

#Load in the data frame to be plotted and remove the row numbers column

tmmfinal <- read.csv("~/Desktop/Mouse-Tumor-Neoantigen/Data/tmmfinal.csv")
tmmfinal <- select(tmmfinal, -1)

#Combine all read counts into one column for plotting
combined <- tmmfinal %>%
  pivot_longer(cols = names (tmmfinal)[1:8],
               values_to = "tmmcounts") %>%
  select(tmmcounts, countsfortmm.gene)

#Add Allele status to each tmm normalized count for plotting

allele <- c("WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT",
            "WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT",
            "WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT",
            "WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT","WT", "MT",
            "WT", "MT","WT", "MT","WT", "MT","WT", "MT")

combined$allele <- allele



# plot a stacked bar graph to show the difference in expression of wildtype and mutant TMM normalized counts

theme_set(theme_minimal())
theme_update(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(color = "grey80", linewidth = .4),
  axis.ticks.x = element_line(color = "grey80", linewidth = .4),
  axis.title.y = element_blank(),
  axis.title.x = element_text(vjust = -2),
  plot.margin = margin(10, 15, 20, 15)
)


ggplot(combined, aes(y=tmmcounts, x=countsfortmm.gene, color = allele))+
geom_boxplot(
  aes(color = allele,
      color = after_scale(darken(colour, .1, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4))),
  width = .42,
  outlier.shape = NA
)+
  labs(y  ="Normalized TMM Counts")+
  ggbeeswarm::geom_quasirandom(
    ## draw bigger points
    size = 1.5,
    ## add some transparency
    alpha = .4,
    ## control range of the beeswarm
    width = .2
  )+ 
  scale_color_brewer(palette = "Set1")+
  coord_flip()+
  theme(legend.box.background = element_rect(size =0.5, color = "black"))+
  guides(color = guide_legend(title = "Allele"))








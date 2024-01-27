#Required packages
library(tidyverse)
library(EnvStats)
library(ggpubr)
library(ggdist)
library(rstatix)
library(ggtext)
library(ggridges)
library(colorspace)
library(ggsci)

#Load data frame
TCGA.gbm <- read_tsv(file = "./Desktop/Turbo TCGA/GBM/TCGA_GBM_normal.tsv")
TCGA.gbm.copynum <- read_tsv(file = "./Desktop/Turbo TCGA/GBM/TCGA_GBM_copynum.tsv")
TCGA.gbm.clinical <- read_tsv(file = "./Desktop/Turbo TCGA/GBM/TCGA_GBM_clinical.tsv")

TCGA.ov.clinical <- TCGA.ov.clinical %>% filter(!is.na(clinical.stage))

TCGA.brca.clinical <- TCGA.ov.clinical %>% filter(clinical.stage != "Stage IC")

TCGA.gbm <- TCGA.gbm %>% filter(sample.type != "Recurrent Tumor")
#Annotate the comparisson groups
gbm.comparissons <- list(c("Normal Tissue","Primary Tumor"))

#Rename column titles to remove the underscores
names(TCGA.gbm)[names(TCGA.gbm) == "_sample_type"] <- "sample.type"
names(TCGA.gbm.copynum)[names(TCGA.gbm.copynum) == "DFNA5 (average)"] <- "DFNA5"
names(TCGA.brca.clinical)[names(TCGA.brca.clinical) == "DFNA5 (average)"] <- "DFNA5"
names(TCGA.brca.clinical)[names(TCGA.brca.clinical) == "tumor_stage.diagnoses"] <- "clinical.stage"
names(TCGA.gbm.clinical)[names(TCGA.gbm.clinical) == "neoplasm_histologic_grade"] <- "grade"
names(TCGA.brca.clinical)[names(TCGA.brca.clinical) == "ENSG00000105928.12"] <- "GSDME"
names(TCGA.gbm.clinical)[names(TCGA.gbm.clinical) == "DFNA5"] <- "GSDME"


#Creats a variable that stores the number of samples per sameple type
add_sample.gbm <- function(x){
  return(c(y=max(x) + .025,
           label = length(x)))
}

#Run some stats 
stat.test.gbm <- TCGA.gbm %>% pairwise_wilcox_test(DFNA5 ~ sample.type)

stats <- stat.test.gbm %>% add_xy_position(fun = "mean_se", x = "sample.type")

#define a color pallette
pal <- c("#A034F0", "#159090", "#B2DF8A", "#33a02c")

#Generate plots using ggplot2 for gene expression
ggplot(TCGA.gbm, aes(x = sample.type, y = DFNA5))+
  ggdist::stat_halfeye(
    aes(colour = sample.type,
        fill = after_scale(lighten(color, .5))),
    adjust =.5,
    width = .45,
    justification = -.4,
    .width = 0,
    point_color = NA
  )+
  geom_boxplot(
    aes(color = sample.type,
        color = after_scale(darken(colour, .1, space = "HLS")),
        fill = after_scale(desaturate(lighten(color, .8), .4))),
    width = .42,
    outlier.shape = NA
  )+
  geom_point(
    aes(color = sample.type,
        color = after_scale(darken(colour, .1, space = "HLS"))),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 2,
    position = position_jitter(seed = 1, width = .12)
  )+
  geom_point(
    aes(fill = sample.type),
    color = "transparent",
    shape = 21,
    stroke = .4,
    size = 2,
    alpha = .3, 
    position = position_jitter(seed =1, width =.12)
  )+ 
  stat_summary(
    geom = "text",
    fun.data = add_sample.gbm,
    aes(label = paste("n=", ..label..)),
    size = 5,
    hjust = -1.5,
    vjust = 0.25
  )+
  coord_flip(xlim = c(1.2, 2), clip = "off")+
  scale_y_continuous(
    limits = c(5,15),
    expand = c(0.1,0.1)
  )+
  scale_x_discrete(expand = c(0.3,0.3))+
  scale_color_manual(values = pal, guide = "none")+
  scale_fill_manual(values = pal, guide = "none")+
  labs(
    x = NULL,
    y = "GSDME Expression Level [Log<sub>2</sub>(FKPM-UQ+1)]",
    title = ""
  )+
  stat_pvalue_manual(stats, label = "p = {p}", tip.length = 0.01, 
                     coord.flip = TRUE, 
                     angle = 0, 
                     hjust = -.1, 
                     bracket.nudge.y = 1.5,
                     size = 5,
                     bracket.size = 0.5)+
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_markdown(size = 12),
    axis.ticks.x = element_line(),
    axis.text = element_text(size = 14, color = "black"),
    title = element_markdown(size = 21, vjust = 2),
    axis.line = element_line(color = "black")
  )

##graph for copy number variation
TCGA.gbm.copynum <- TCGA.gbm.copynum %>%
  mutate(mycolor = if_else(DFNA5>0, "#E64B35FF", "#4DBBD5FF"))


ggplot(TCGA.gbm.copynum, aes(x = reorder(sample, -DFNA5), y = DFNA5))+
  geom_segment(aes(x = reorder(sample, -DFNA5),
                   xend = reorder(sample, -DFNA5),
                   y = 0, yend = DFNA5,
                   color = mycolor))+
  geom_point(aes(color = mycolor),
             fill = "white",
             alpha = .7
  )+
  scale_color_discrete(labels = c("Copy Loss", "Copy Gain"))+
  labs(
    x = "Patients",
    y = "Gasdermin E<br>[log<sub>2</sub>(copy number/2)]</br>"
  )+
  theme_minimal()+
  theme(
    panel.grid.major.y = element_line(colour = "#E6E6E6"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor.x = element_line(),
    axis.title.y = element_markdown(size = 12),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = c(0.9,0.73),
    legend.title = element_blank(),
    legend.text = element_text(size =14),
    legend.key.size = unit(1, "cm")
  )

#plot of gene expression across stage and grade using ridge plots.
ggplot(TCGA.brca.clinical, aes(x = GSDME, y = clinical.stage, fill = clinical.stage))+
  geom_density_ridges(scale = 3)+
  scale_y_discrete(expand = c(0.0,0))+
  scale_x_continuous(expand = c(0.01,0))+
  scale_fill_npg()+
  coord_cartesian(clip = "off")+
  labs(
    x = "GSDME Expression Level <br>[Log<sub>2</sub>(FKPM-UQ+1)]</br>"
  )+
  theme_ridges(center_axis_labels = TRUE)+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_markdown(size = 11),
    axis.title.y = element_blank()
  )

#ridge plot of gene expression across grade
ggplot(TCGA.gbm.clinical, aes(x = GSDME, y = grade, fill = grade))+
  geom_density_ridges(scale = 3)+
  scale_y_discrete(expand = c(0.0,0))+
  scale_x_continuous(expand = c(0.01,0))+
  scale_fill_npg()+
  coord_cartesian(clip = "off")+
  labs(
    x = "GSDME Expression Level <br>[Log<sub>2</sub>(FKPM-UQ+1)]</br>"
  )+
  theme_ridges(center_axis_labels = TRUE)+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_markdown(size = 11),
    axis.title.y = element_blank()
  )

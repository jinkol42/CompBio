y <- DGEList(Total_counts)
calcNormFactors(y, method = "TMM")
tmm <- cpm(y)
tmm <- as.data.frame(tmm)

tmm$gene <- Total_counts$Gene

write.csv(tmm, file = "./tmmnormalizedcounts.csv")


list_df = list(MC38Ad9MerVAF, MC38ADFLVAF, MC38gp33byvariant, MC38WTVAF)
vaf <- list_df %>% Reduce(full_join, by = 'position')

alldata <- MC38WTVAF %>%
  left_join(MC38gp33byvariant, by='position') %>%
  left_join(MC38Ad9MerVAF, by='position') %>%
  left_join(MC38ADFLVAF, by='position')

count.vaf <- select(alldata, contig.x,position,refCount.x,altCount.x,refCount.y,altCount.y,refCount.x.x,altCount.x.x,refCount.y.y,altCount.y.y)

rownames(count.vaf) <- count.vaf[,1]
count.vaf[,-1]

count.vaf <- count.vaf %>%
  rename (
    "wt.ref" = "refCount.x",
    "wt.alt" = "altCount.x",
    "gp33.ref" = "refCount.y",
    "gp33.alt" = "altCount.y",
    "Ad9mer.ref" = "refCount.x.x",
    "Ad9mer.alt" = "altCount.x.x",
    "Adfl.ref" = "refCount.y.y",
    "Adfl.alt" = "altCount.y.y"
  )

subset.mc38 <- select(count.vaf, contig.x, wt.ref, wt.alt, gp33.ref, gp33.alt )
subset.mc38.2 <- select(count.vaf, contig.x, position, wt.ref, wt.alt, gp33.ref, gp33.alt )
subset.mc38.3 <- select(count.vaf, contig.x, wt.ref, gp33.ref)
subset.mc38.4 <- select(count.vaf, contig.x, wt.alt, gp33.alt)
subset.mc38.5 <- select(count.vaf, contig.x, position,  wt.ref, wt.alt, gp33.ref, gp33.alt, Ad9mer.ref, Ad9mer.alt, Adfl.ref, Adfl.alt)


subset.mc38 <- subset.mc38 %>% na.omit()
subset.mc38.2 <- subset.mc38.2 %>% na.omit()
subset.mc38.3 <- subset.mc38.3 %>% na.omit()
subset.mc38.4 <- subset.mc38.4 %>% na.omit()
subset.mc38.5 <- subset.mc38.5 %>% na.omit()


x <- DGEList(subset.mc38)
calcNormFactors(x, method = "TMM")
tmm.mc38 <- cpm(x)
tmm.mc38 <- as.data.frame(tmm.mc38)

y <- DGEList(subset.mc38.3)
calcNormFactors(y, method = "TMM")
tmm.mc38.ref <- cpm(y)
tmm.mc38.ref <- as.data.frame(tmm.mc38.ref)

z <- DGEList(subset.mc38.5)
calcNormFactors(z, method = "TMM")
tmm.mc38.final <- cpm(z)
tmm.mc38.final <- as.data.frame(tmm.mc38.final)

t <- DGEList(countsfortmm)
calcNormFactors(t, method = "TMM")
tmm.final <- cpm(t)
tmm.final <- as.data.frame(tmm.final)

tmm.final <- cbind(tmm.final, countsfortmm$gene)
write.csv(tmm.final, file = "./tmmfinal.csv")

tmm.mc38.final <- cbind(tmm.mc38.final, subset.mc38.5$contig.x)
tmm.mc38.final <- cbind(tmm.mc38.final, subset.mc38.5$position)

write.csv(tmm.mc38.final, file = "./tmmfinalcounts.csv")
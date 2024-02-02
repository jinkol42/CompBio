B16.wherry <- select(B16, c(NAME, B16.1:B16.4))
rownames(B16.wherry) = make.names(B16.wherry$NAME, unique = TRUE)
B16.wherry <- B16.wherry[,-1]

log.B16.wherry <- log2(B16.wherry)

rownames(TExh_Normcounts) = make.names(TExh_Normcounts$X, unique = TRUE)
Tex <- TExh_Normcounts[,-1]

Comp <- merge(Tex, log.B16.wherry,
              by = "row.names", all = TRUE)

Comp <- na.omit(Comp)

rownames(Comp) = make.names(Comp$Row.names, unique = TRUE)
Comp <- Comp[,-1]

cond.ex <- factor(c("Prog2", "Prog2", "Prog2", "Int", "Int", "Int", "Prog1", "Prog1", "Prog1", "Term", "Term", "Term", "B16", "B16", "B16", "B16"))

coldata.ex <- data.frame(row.names = colnames(Comp), cond.ex)
dds.ex <- DESeqDataSetFromMatrix(countData = round(Comp), colData = coldata.ex, design = ~cond.ex)
vts.ex <- vst(dds.ex, blind = FALSE, fitType = c("mean"))
pcafig.ex <- plotPCA(vts.ex, intgroup = "cond.ex")
datapca.ex <- as.data.frame(pcafig.ex)




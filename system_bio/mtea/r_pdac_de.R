suppressMessages(library(limma))
lines <- c("CAPAN2","CFPAC1","MiaPaca2","PANC1")
dir.create("out", showWarnings=FALSE)
for (ln in lines) {
  f <- sprintf("pdac/GSE238198_Pdac.%s.genes.normalizedCounts.txt.gz", ln)
  d <- read.delim(gzfile(f), check.names=FALSE, row.names=1)
  y <- log2(as.matrix(d) + 1)
  tp <- sapply(strsplit(colnames(y), "_"), function(z) z[2])
  tp <- factor(c(NT="NT",`3h`="h3",`12h`="h12",`24h`="h24")[tp], levels=c("NT","h3","h12","h24"))
  design <- model.matrix(~ 0 + tp); colnames(design) <- levels(tp)
  keep <- rowMeans(y) > log2(1 + 5)   # mean normalised count > 5
  fit <- lmFit(y[keep, ], design)
  cm <- makeContrasts(h3=h3-NT, h12=h12-NT, h24=h24-NT, levels=design)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend=TRUE, robust=TRUE)
  for (k in colnames(cm)) {
    tt <- topTable(fit2, coef=k, number=Inf, sort.by="none")
    tt$gene <- rownames(tt)
    write.table(tt, sprintf("out/pdac_de_%s_%s.tsv", ln, k), sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("%-9s %-4s genes %5d | padj<.05: %5d | |logFC|>=1 & padj<.05: %4d\n", ln, k, nrow(tt), sum(tt$adj.P.Val<.05), sum(tt$adj.P.Val<.05 & abs(tt$logFC)>=1)))
  }
}

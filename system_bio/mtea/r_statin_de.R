suppressMessages(library(limma))
d <- read.delim("pdac/statin_log2fpkm.tsv", row.names=1, check.names=FALSE)
y <- as.matrix(d)
parts <- strsplit(colnames(y), ".", fixed=TRUE)
cell <- sapply(parts, function(p) p[1])
grp <- sapply(parts, function(p) if (length(p) == 2) "Control" else paste(p[2], p[3], sep="_"))
for (ln in unique(cell)) {
  s <- cell == ln; yy <- y[, s]; g <- factor(grp[s]); g <- relevel(g, "Control")
  design <- model.matrix(~ 0 + g); colnames(design) <- levels(g)
  fit <- lmFit(yy, design)
  cons <- setdiff(levels(g), "Control")
  cm <- makeContrasts(contrasts=paste0(cons, "-Control"), levels=design)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend=TRUE, robust=TRUE)
  for (k in seq_along(cons)) {
    tt <- topTable(fit2, coef=k, number=Inf, sort.by="none"); tt$gene <- rownames(tt)
    write.table(tt, sprintf("out/statin_de_%s_%s.tsv", ln, cons[k]), sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("%-9s %-22s padj<.05: %5d | |logFC|>=0.5: %4d\n", ln, cons[k], sum(tt$adj.P.Val<.05), sum(tt$adj.P.Val<.05 & abs(tt$logFC)>=0.5)))
  }
}

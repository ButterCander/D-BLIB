suppressMessages({library(limma); library(edgeR)})
cnt <- read.csv(gzfile("GSE285616_counts.csv.gz"), check.names=FALSE)
cnt <- cnt[!duplicated(cnt$Symbol), ]; rownames(cnt) <- cnt$Symbol; cnt$Symbol <- NULL
cnt <- as.matrix(cnt[rowSums(cnt) >= 10, ])
stopifnot(ncol(cnt) == 24)
lev <- c("c0","PI","PD","TAK","PITAK","PIPD")
cond <- factor(lev[(0:23 %% 6) + 1], levels=lev); batch <- factor(paste0("b", (0:23 %/% 6) + 1))
design <- model.matrix(~ batch + cond)
y <- DGEList(counts=cnt); y <- calcNormFactors(y, method="TMM")
v <- voom(y, design, plot=FALSE)
tg <- read.delim("out/task_genes_long.tsv")
idx <- lapply(split(tg$gene, tg$task), function(g) which(rownames(cnt) %in% g))
idx <- idx[sapply(idx, length) >= 5 & sapply(idx, length) <= 200]
cat("tasks tested:", length(idx), " distinct sets:", length(unique(lapply(idx, sort))), "\n")
cn <- colnames(design); e <- function(...) { w <- c(...); out <- setNames(numeric(length(cn)), cn); for (n in names(w)) out[paste0("cond", n)] <- w[[n]]; out }
CON <- list(
  PI=e(PI=1), PD=e(PD=1), TAK=e(TAK=1), PITAK=e(PITAK=1), PIPD=e(PIPD=1),
  PIPD_vsPI=e(PIPD=1,PI=-1), PIPD_vsPD=e(PIPD=1,PD=-1), PIPD_dose=e(PIPD=1,PI=-.5,PD=-.5), PIPD_add=e(PIPD=1,PI=-1,PD=-1),
  PITAK_vsPI=e(PITAK=1,PI=-1), PITAK_vsTAK=e(PITAK=1,TAK=-1), PITAK_dose=e(PITAK=1,PI=-.5,TAK=-.5), PITAK_add=e(PITAK=1,PI=-1,TAK=-1))
res <- list()
for (nm in names(CON)) {
  cvec <- CON[[nm]]
  for (meth in c("camera_fixed","camera_est","fry")) {
    r <- switch(meth,
      camera_fixed = camera(v, idx, design, contrast=cvec, inter.gene.cor=0.01),
      camera_est   = camera(v, idx, design, contrast=cvec, inter.gene.cor=NA),
      fry          = fry(v, idx, design, contrast=cvec))
    r$task <- rownames(r); r$q <- p.adjust(r$PValue, "BH"); r$contrast <- nm; r$method <- meth
    res[[paste(nm, meth)]] <- r
    write.table(r, sprintf("out/r_%s_%s.tsv", meth, nm), sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("%-11s %-13s p<.05: %3d (down %3d) | BH<.05: %3d\n", nm, meth, sum(r$PValue<.05), sum(r$PValue<.05 & r$Direction=="Down"), sum(r$q<.05)))
  }
}

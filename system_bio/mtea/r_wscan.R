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

cat("
w-scan: AB - w*A - w*B
")
for (pair in list(c("PIPD","PI","PD"), c("PITAK","PI","TAK"))) {
  for (w in c(0.5,0.6,0.7,0.8,0.9,1.0)) {
    cvec <- e(setNames(c(1,-w,-w), pair))
    out <- sapply(c("camera_fixed","camera_est","fry"), function(m) {
      r <- switch(m, camera_fixed=camera(v, idx, design, contrast=cvec, inter.gene.cor=0.01),
                     camera_est=camera(v, idx, design, contrast=cvec, inter.gene.cor=NA),
                     fry=fry(v, idx, design, contrast=cvec))
      q <- p.adjust(r$PValue, "BH"); sprintf("p<.05 %3d (down %3d) BH %3d", sum(r$PValue<.05), sum(r$PValue<.05 & r$Direction=="Down"), sum(q<.05)) })
    cat(sprintf("%-6s w=%.1f | %s
", pair[1], w, paste(names(out), out, sep=": ", collapse=" | ")))
  }
}

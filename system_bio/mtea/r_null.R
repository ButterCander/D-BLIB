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
cn <- colnames(design); e <- function(...) { w <- c(...); out <- setNames(numeric(length(cn)), cn); for (n in names(w)) out[paste0('cond', n)] <- w[[n]]; out }

set.seed(7); B <- 40
fp <- list(); bhc <- list()
for (b in seq_len(B)) {
  cs <- unlist(lapply(0:3, function(k) sample(lev)))   # shuffle condition labels within each batch
  condp <- factor(cs, levels=lev)
  dp <- model.matrix(~ batch + condp); colnames(dp) <- sub("condp", "cond", colnames(dp))
  vp <- voom(y, dp, plot=FALSE)
  cvs <- list(single=e(PI=1), dose=e(PIPD=1,PI=-.5,PD=-.5), vsPI=e(PIPD=1,PI=-1))
  for (nm in names(cvs)) for (m in c("camera_fixed","camera_est","fry")) {
    r <- switch(m, camera_fixed=camera(vp, idx, dp, contrast=cvs[[nm]], inter.gene.cor=0.01),
                   camera_est=camera(vp, idx, dp, contrast=cvs[[nm]], inter.gene.cor=NA),
                   fry=fry(vp, idx, dp, contrast=cvs[[nm]]))
    k <- paste(m, nm); fp[[k]] <- c(fp[[k]], mean(r$PValue < .05)); bhc[[k]] <- c(bhc[[k]], sum(p.adjust(r$PValue,"BH") < .05))
  }
}
cat("
Label-shuffle null (", B, " shuffles within batch; 120 task sets): false-positive rate at p<.05 | mean BH<.05 rejections per null contrast
")
for (k in names(fp)) cat(sprintf("%-22s FPR %.3f | BH rejections %.1f
", k, mean(fp[[k]]), mean(bhc[[k]])))

suppressMessages({library(edgeR)})
cnt <- read.csv(gzfile("GSE285616_counts.csv.gz"), row.names=1, check.names=FALSE)
cond <- rep(c("ctrl","PI","PD","TAK","PITAK","PIPD"), 4)
batch <- factor(rep(1:4, each=6))
cond <- factor(cond, levels=c("ctrl","PI","PD","TAK","PITAK","PIPD"))
y <- DGEList(counts=as.matrix(cnt))
keep <- filterByExpr(y, group=cond)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~batch+cond)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
cat("genes tested:", nrow(y), "\n")
res <- list()
for (cn in c("condPI","condPD","condTAK","condPITAK","condPIPD")) {
  q <- glmQLFTest(fit, coef=cn)
  tt <- topTags(q, n=Inf)$table
  res[[cn]] <- tt
  deg <- sum(tt$FDR<0.05 & abs(tt$logFC)>1)
  up <- sum(tt$FDR<0.05 & tt$logFC>1); dn <- sum(tt$FDR<0.05 & tt$logFC< -1)
  cat(cn, "DEG(|lfc|>=1, FDR<.05):", deg, "up", up, "down", dn, "\n")
}
genes <- c("HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS","SCD","SREBF1","SREBF2","INSIG1","LDLR","SQLE","LSS","ASNS","PHGDH","PSAT1","ODC1","RRM1","RRM2","CTPS1","TDO2","IDO1")
out <- NULL
for (cn in names(res)) {
  tt <- res[[cn]]
  g <- genes[genes %in% rownames(tt)]
  out <- rbind(out, data.frame(cond=cn, gene=g, logFC=round(tt[g,"logFC"],2), FDR=signif(tt[g,"FDR"],2)))
}
w <- reshape(out[,c("cond","gene","logFC")], idvar="gene", timevar="cond", direction="wide")
print(w, row.names=FALSE)
w2 <- reshape(out[,c("cond","gene","FDR")], idvar="gene", timevar="cond", direction="wide")
print(w2, row.names=FALSE)

# Prepare data for genome-scale metabolic modeling (GEM)

library(my.utils)
library(gembox)
data("recon1")

load("../data/Fraietta.RData")
dat <- list(mock=mat.mock, stim=mat.stim)
dat <- lapply(dat, function(x) {
  x <- x[rowMeans(x==0)<0.9, ]
  prep.data(x, norm.method="quantile")
}) # log-transformed

exprs.int <- lapply(dat, function(x) {
  list(r=exprs2int(recon1, x[, responders]),
       nr=exprs2int(recon1, x[, !responders]))
}) # 236 genes not in data

phe <- data.table(grp=factor(responders, levels=c(FALSE,TRUE)))
de.res <- lapply(dat, function(x) de(x, phe, coef="grpTRUE", robust=TRUE, trend=TRUE))
dflux.int <- lapply(de.res, function(x) de.dt2dflux(recon1, x, topn=80, padj.cutoff=1.1))
de.res1 <- lapply(de.res, function(x) x[, log.fc:=-log.fc])
dflux.int.rev <- lapply(de.res1, function(x) de.dt2dflux(recon1, x, topn=80, padj.cutoff=1.1))

save(exprs.int, dflux.int, dflux.int.rev, file="data.for.gem.RData")



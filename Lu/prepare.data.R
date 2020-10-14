# Prepare the single-cell gene expression data of persistent and non-persistent T cell clones from Lu et al. for genome-scale metabolic modeling (GEM)

library(SingleCellExperiment)
library(scater)
library(scran)
library(my.utils)
library(gembox)
data("recon1")

load("../data/Lu.RData")


### process single cell data

# remove low quality cells
df <- perCellQCMetrics(sce)
reasons <- quickPerCellQC(df)
sum(reasons$discard) # 2
discard <- df$sum<2000 | df$detected<200
sum(discard) # [1] 21511
sce.sub <- sce[, !discard]

# normalization by deconvolution
set.seed(1)
clust <- quickCluster(sce.sub)
sf <- calculateSumFactors(sce.sub, cluster=clust)
sce.sub <- logNormCounts(sce.sub)


### get clone barcodes

# clone 9.1-NP, non-persistent, d12
id.np9 <- unique(unlist(lapply(tcr[c("PBL4095-d12-9-rxn-A TCR","PBL4095-d12-9-rxn-B TCR","PBL4095-d12-9-rxn-C TCR")], function(x) x[cdr3=="CASSLGEGRVDGYTF" & raw_clonotype_id %in% paste0("clonotype",1:2), barcode])))

# clone 9.2-P, persistent, d12
id.p9 <- unique(c(
  tcr[["PBL4095-d12-9-rxn-A TCR"]][cdr3=="CASSFGQSSTYGYTF" & raw_clonotype_id %in% paste0("clonotype",c(4,5,8,10)), barcode],
  tcr[["PBL4095-d12-9-rxn-B TCR"]][cdr3=="CASSFGQSSTYGYTF" & raw_clonotype_id %in% paste0("clonotype",c(3,5,7,9)), barcode],
  tcr[["PBL4095-d12-9-rxn-C TCR"]][cdr3=="CASSFGQSSTYGYTF" & raw_clonotype_id %in% paste0("clonotype",c(7,9,10,28)), barcode]
))

# clone 10, persistent, d12
id.p10 <- tcr[["PBL4095-d12-10 TCR"]][cdr3=="CASSDPGTEAFF" & raw_clonotype_id %in% paste0("clonotype",1:4), unique(barcode)]

# get expression data of the clones
get.med.expr <- function(ids, f0=0.9) {
  mat <- logcounts(sce.sub)[,colData(sce.sub)$Barcode %in% ids]
  print(ncol(mat))
  rownames(mat) <- rowData(sce.sub)$Symbol
  idx <- rowMeans(mat==0)<f0
  print(sum(idx))
  apply(mat[idx,], 1, function(x) median(x[x!=0]))
}
med.expr.np9 <- get.med.expr(id.np9) # 3999, 3829
med.expr.p9 <- get.med.expr(id.p9) # 333, 3706
med.expr.p10 <- get.med.expr(id.p10) # 1455, 4069

# data for GEM
expr.int.np9 <- exprs2int(recon1, med.expr.np9)
expr.int.p9 <- exprs2int(recon1, med.expr.p9)
expr.int.p10 <- exprs2int(recon1, med.expr.p10)

save(list=ls(pattern="^med\\.expr|expr\\.int"), file="data.for.gem.RData")


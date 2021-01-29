library(my.utils)
library(survival)

dat <- readRDS("../data/tcga.xena.tmm.log.cpm.by.cancer.type.RDS")
tcga.purity <- readRDS("../data/tcga.purity.RDS")
for (i in names(dat)) {
  if (toupper(i) %in% names(tcga.purity)) {
    dat[[i]]$pheno$purity <- tcga.purity[[toupper(i)]]$purity[match(dat[[i]]$pheno$sample.id, tcga.purity[[toupper(i)]]$sample.id)]
  }
}


### check the correlation between UCP2 gene expression and several T cell memory and stemness genes, after controlling for tumor purity, in each TCGA cancer type

check.cor <- function(gn1, gn2="UCP2") {
  dat1 <- dat[sapply(dat, function(x) "purity" %in% names(x$pheno))]
  res <- rbindlist(lapply(dat1, function(x) {
    #tmp <- cor.test(x$expr[x$geneid$symbol==gn1,], x$expr[x$geneid$symbol==gn2,], method="s")
    #data.table(rho=tmp$estimate, pval=tmp$p.value)
    dt1 <- data.table(g1=x$expr[x$geneid$symbol==gn1,], g2=x$expr[x$geneid$symbol==gn2,], purity=x$pheno$purity)
    res <- coef(summary(lm(g2 ~ g1 + purity, data=dt1)))
    res <- data.table(coef=res["g1", "Estimate"], pval=res["g1", "Pr(>|t|)"])
  }), idcol="cancer.type")
  res[, cancer.type:=toupper(cancer.type)]
  res[, padj:=p.adjust(pval, "BH")]
  res[order(padj,pval)]
}

gns <- cn("CD8A","TCF7","CD27","SELL")
res <- lapply(gns, check.cor)
saveRDS(res, file="Tm.gene.cor.RDS")


### check the association between UCP2 gene expression and cancer patient survival, in each TCGA cancer type

surv.res <- rbindlist(lapply(dat, function(x) {
  dt1 <- copy(x$pheno)
  dt1[, y:=x$expr[x$geneid$symbol=="UCP2",]]
  dt1[, stage1:=grepl("iii|iv", patho.stage)]
  dt1[, race1:=ifelse(race=="black or african american", "aa","non-aa")]
  dt1[is.na(race1), race1:="unknown"]
  dt1 <- dt1[!grepl("c4", subtype.immune)]
  tryCatch({
    cox.res <- coef(summary(coxph(Surv(os.days, os) ~ y + age + stage1 + strata(gender, race1), data=dt1)))
    data.table(coef=cox.res["y","coef"], pval=cox.res["y","Pr(>|z|)"])
  }, error=function(e) data.table(coef=NA, pval=NA))
}), idcol="cancer.type")
surv.res[, padj:=p.adjust(pval, "BH")]
surv.res <- surv.res[order(padj,pval)]

saveRDS(surv.res, file="surv.res.RDS")


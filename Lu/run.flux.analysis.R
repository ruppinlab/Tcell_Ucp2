# Run metabolic flux analysis for the persistent and non-persistent T cell clones from Lu et al.

library(gembox)
data("recon1")
nc <- 16L # number of cores to use, adjust as appropriate

load("data.for.gem.RData")

# run the iMAT algorithm to compute metabolic fluxes for the persisent clone (9.2-P) and non-persistent clone (9.1-NP), respectively
imat.res <- list(
  p9=imat(recon1, expr.int.p9, samp.pars=list(nc=nc, n.sample=5000)),
  np9=imat(recon1, expr.int.np9, samp.pars=list(nc=nc, n.sample=5000))
)
saveRDS(imat.res, file="imat.res.RDS")

# differential flux analysis of the persistent clone compared to the non-persistent clone
df.res <- get.diff.flux(imat.res$np9$result.model, imat.res$p9$result.model, nc=nc)

# perform metabolic pathway enrichment based on the differential flux result
rset <- subsystems2gsets(recon1)
df.gsea.res <- pathway.gsea(df.res, rset)
save(df.res, df.gsea.res, file="dflux.res.RData")


# Run flux analysis for responders and non-responders to anti-CD19 CAR-T therapy from Fraietta et al.

library(gembox)
data("recon1")
nc <- 16L

load("./data.for.gem.RData")

bm <- get.biomass.idx(recon1)
model <- set.rxn.bounds(recon1, bm, lbs=0.1, relative=TRUE)

# run the iMAT algorithm to compute metabolic fluxes in responders and non-responders
imat.res <- lapply(exprs.int$mock, function(x) imat(model, x, samp.pars=list(nc=nc, n.sample=5000)))
saveRDS(imat.res, file="imat.res.RDS")

# perform differential flux analysis for responders compared to non-responders
df.res <- get.diff.flux(imat.res$nr$result.model, imat.res$r$result.model, nc=nc)

# perform metabolic pathway enrichment based on the differential flux result
rset <- subsystems2gsets(recon1)
df.gsea.res <- pathway.gsea(df.res, rset)
save(df.res, df.gsea.res, file="dflux.res.RData")


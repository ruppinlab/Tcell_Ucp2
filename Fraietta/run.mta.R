# Run the MTA algorithm for predicting metabolic reaction knockouts that may result in non-responsiveness to anti-CD19 CAR-T therapy, using the Fraietta et al. dataset.

library(gembox)
data("recon1")
nc <- 16L # number of cores to use, adjust as appropriate

load("./data.for.gem.RData")

bm <- get.biomass.idx(recon1)
model <- set.rxn.bounds(recon1, bm, lbs=0.05, relative=TRUE)

imat.res <- lapply(exprs.int$mock, function(x) imat(model, x, samp.pars=list(nc=nc, n.sample=5000)))
vref <- rowMeans(imat.res$r$result.model$sample$pnts[, 2001:5000])
saveRDS(vref, file="vref.for.mta.RDS")
mta.raw.res <- rmetal(model, vref, dflux.int.rev$mock, nc=nc, detail=FALSE)

# clean up results from raw output from MTA
mta.res <- mta.raw.res$result.rmetal[!is.na(rTS), .(id, reaction=rxn, score=rTS, percent.rank=rank(-rTS, na.last="keep")/(.N-1)*100)]
mta.res <- mta.res[order(percent.rank)]

# metabolic pathway enrichment of top 10% hits
rset <- subsystems2gsets(recon1)
enr.res <- enrich.gsets(mta.res[id!=0 & percent.rank<10 & score>score[id==0], reaction], rset, recon1$rxns, nc=nc)

save(mta.res, enr.res, file="mta.res.RData")


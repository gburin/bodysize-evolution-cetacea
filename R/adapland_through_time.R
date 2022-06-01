library("bayou")
library("here")

here::i_am("R/adapland_time.R")

prior.per.branch <- function(chain, k){
    prior <- (1 - (chain$tree$edge.length / sum(chain$tree$edge.length)))^k
    return(1 - prior)
}

full.tree <- read.tree(here::here("data/final_tree_cetacea.tre"))

oligo.tree <- read.tree(here::here("data/adapland_timeslices/Oligocene_tree.tre"))
earlymio.tree <- read.tree(here::here("data/adapland_timeslices/EarlyMiocene_tree.tre"))
midmio.tree <- read.tree(here::here("data/adapland_timeslices/MidMiocene_tree.tre"))
latemio.tree <- read.tree(here::here("data/adapland_timeslices/LateMiocene_tree.tre"))
plio.tree <- read.tree(here::here("data/adapland_timeslices/Pliocene_tree.tre"))

fulldata <- read.csv(here::here("data/data_imputation.csv"), na.strings = "NA")
fulldata$species[grep("Kekenodon", fulldata$species)] <- "Kekenodon_onamata"
fulldata <- fulldata[-grep("_cf", fulldata$species), ]
fulldata <- fulldata[-grep("_indet", fulldata$species), ]
fulldata <- fulldata[-c(grep("_sp$", fulldata$species), grep("_sp[^i]", fulldata$species)), ]

## Fixing name mismatch between tree and data
fulldata$species[grep("Eodelphis_kabatensis", fulldata$species)] <- "Eodelphinus_kabatensis"
fulldata$species[grep("Rodhocetus_kasrani", fulldata$species)] <- "Rodhocetus_kasranii"
fulldata$species[grep("Lagenorhynchus_obliquidens", fulldata$species)] <- "Sagmatias_obliquidens"
fulldata[nrow(fulldata) + 1, ] <- NA
## Tursiops australis is missing from data, so we will imput it solely from phylogeny
fulldata[nrow(fulldata), 1] <- "Tursiops_australis"
fulldata[nrow(fulldata), ncol(fulldata)] <- "extant"
rownames(fulldata) <- NULL

logdata <- fulldata[, -ncol(fulldata)]
logdata[, 2:ncol(logdata)] <- log10(logdata[, 2:ncol(logdata)])

## Using fuzzy matching to check for typos in species names
library("stringdist")
dist.name <- adist(unique(logdata$species), full.tree$tip.label, partial = TRUE, ignore.case = TRUE)
min.name <- apply(dist.name, 1, min)

match.names <- c()
for(i in 1:nrow(dist.name)){
    s2.i <- match(min.name[i], dist.name[i,])
    s1.i <- i
    match.names <- rbind(data.frame(s2.i = s2.i, s1.i = s1.i, s2name = full.tree$tip.label[s2.i], s1name = unique(logdata$species)[s1.i], adist = min.name[i]), match.names)
}

sp.to.change.data <- match.names$s1name[match.names$adist == 1]
for(i in 1:length(sp.to.change.data)){
    logdata$species[which(!is.na(match(logdata$species, sp.to.change.data[i])))] <- match.names$s2name[match.names$s1name == sp.to.change.data[i]]
}

full.tree.pruned <- drop.tip(full.tree, full.tree$tip.label[!(full.tree$tip.label %in% unique(logdata$species))])
oligo.tree.pruned <- drop.tip(oligo.tree, oligo.tree$tip.label[!(oligo.tree$tip.label %in% unique(logdata$species))])
earlymio.tree.pruned <- drop.tip(earlymio.tree, earlymio.tree$tip.label[!(earlymio.tree$tip.label %in% unique(logdata$species))])
midmio.tree.pruned <- drop.tip(midmio.tree, midmio.tree$tip.label[!(midmio.tree$tip.label %in% unique(logdata$species))])
latemio.tree.pruned <- drop.tip(latemio.tree, latemio.tree$tip.label[!(latemio.tree$tip.label %in% unique(logdata$species))])
plio.tree.pruned <- drop.tip(plio.tree, plio.tree$tip.label[!(plio.tree$tip.label %in% unique(logdata$species))])


logdata.oligo <- logdata[logdata$species %in% oligo.tree.pruned$tip.label, ]
logdata.earlymio <- logdata[logdata$species %in% earlymio.tree.pruned$tip.label, ]
logdata.midmio <- logdata[logdata$species %in% midmio.tree.pruned$tip.label, ]
logdata.latemio <- logdata[logdata$species %in% latemio.tree.pruned$tip.label, ]
logdata.plio <- logdata[logdata$species %in% plio.tree.pruned$tip.label, ]

#langhian.tree.pruned <- paleotree::dropZLB(langhian.tree.pruned)
#rupelian.tree.pruned <- paleotree::dropZLB(rupelian.tree.pruned)

imput.data <- readRDS(here::here("output/imput_data_withZBL.RDS"))
#rownames(imput.data$anc_recon)[Ntip(slater.tree.imput)] <- "Tursiops_australis"

## Calculating standard error for species with multiple measurements
sd.bayou.per.sp <- aggregate(logdata$tl, by = list(logdata$species), FUN = sd, na.rm = TRUE)
n.per.sp <- aggregate(logdata$tl, by = list(logdata$species), FUN = length)

sd.bayou.per.sp$se <- NA
sd.bayou.per.sp$se[which(!is.na(sd.bayou.per.sp[,2]))] <- sd.bayou.per.sp[which(!is.na(sd.bayou.per.sp[,2])), 2]/sqrt(n.per.sp[which(!is.na(sd.bayou.per.sp[,2])),2])

## Calculating pooling variance to use as standard error for species without multiple measurements
pooled.var <- sum((sd.bayou.per.sp[, 2])^2 * (n.per.sp[, 2] - 1), na.rm = TRUE)/(sum(n.per.sp[, 2], na.rm = TRUE) - sum(!is.na(sd.bayou.per.sp[, 2])))

sd.bayou.per.sp$se[is.na(sd.bayou.per.sp$se)] <- sqrt(pooled.var)/sqrt(n.per.sp[which(is.na(sd.bayou.per.sp[,2])),2])

tl.full <- setNames(imput.data$anc_recon[match(full.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(full.tree.pruned$tip.label, rownames(imput.data$anc_recon))])
tl.oligo <- setNames(imput.data$anc_recon[match(oligo.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(oligo.tree.pruned$tip.label, rownames(imput.data$anc_recon))])
tl.earlymio <- setNames(imput.data$anc_recon[match(earlymio.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(earlymio.tree.pruned$tip.label, rownames(imput.data$anc_recon))])
tl.midmio <- setNames(imput.data$anc_recon[match(midmio.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(midmio.tree.pruned$tip.label, rownames(imput.data$anc_recon))])
tl.latemio <- setNames(imput.data$anc_recon[match(latemio.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(latemio.tree.pruned$tip.label, rownames(imput.data$anc_recon))])
tl.plio <- setNames(imput.data$anc_recon[match(plio.tree.pruned$tip.label, rownames(imput.data$anc_recon)), "tl"], rownames(imput.data$anc_recon)[match(plio.tree.pruned$tip.label, rownames(imput.data$anc_recon))])


se.full <- sd.bayou.per.sp$se[match(names(tl.full), sd.bayou.per.sp[, 1])]
names(se.full) <- names(tl.full)
se.oligo <- sd.bayou.per.sp$se[match(names(tl.oligo), sd.bayou.per.sp[, 1])]
names(se.oligo) <- names(tl.oligo)
se.earlymio <- sd.bayou.per.sp$se[match(names(tl.earlymio), sd.bayou.per.sp[, 1])]
names(se.earlymio) <- names(tl.earlymio)
se.midmio <- sd.bayou.per.sp$se[match(names(tl.midmio), sd.bayou.per.sp[, 1])]
names(se.midmio) <- names(tl.midmio)
se.latemio <- sd.bayou.per.sp$se[match(names(tl.latemio), sd.bayou.per.sp[, 1])]
names(se.latemio) <- names(tl.latemio)
se.plio <- sd.bayou.per.sp$se[match(names(tl.plio), sd.bayou.per.sp[, 1])]
names(se.plio) <- names(tl.plio)


priorOU.oligo <- make.prior(oligo.tree.pruned, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 2, kmax = 1000), dsb = list(bmax = Inf, prob = oligo.tree.pruned$edge.length), 
                                 dtheta = list(mean = mean(tl.full), sd = 1.5 * sd(tl.full))),
                      plot.prior = FALSE
                      )

mcmcOU.oligo <- bayou.makeMCMC(oligo.tree.pruned, tl.oligo, SE = se.oligo, prior = priorOU.oligo, file.dir = here::here("output/bayou/adapland/oligocene/bayou_oligo_fixedSE_wZBL"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.oligo$run(1000000) # Run the MCMC

priorOU.earlymio <- make.prior(earlymio.tree.pruned, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 2, kmax = 1000), dsb = list(bmax = Inf, prob = earlymio.tree.pruned$edge.length), 
                                 dtheta = list(mean = mean(tl.full), sd = 1.5 * sd(tl.full))),
                      plot.prior = FALSE
                      )

mcmcOU.earlymio <- bayou.makeMCMC(earlymio.tree.pruned, tl.earlymio, SE = se.earlymio, prior = priorOU.earlymio, file.dir = here::here("output/bayou/adapland/earlymiocene/bayou_earlymio_fixedSE_wZBL"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.earlymio$run(1000000) # Run the MCMC

priorOU.midmio <- make.prior(midmio.tree.pruned, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 2, kmax = 1000), dsb = list(bmax = Inf, prob = midmio.tree.pruned$edge.length), 
                                 dtheta = list(mean = mean(tl.full), sd = 1.5 * sd(tl.full))),
                      plot.prior = FALSE
                      )

mcmcOU.midmio <- bayou.makeMCMC(midmio.tree.pruned, tl.midmio, SE = se.midmio, prior = priorOU.midmio, file.dir = here::here("output/bayou/adapland/midmiocene/bayou_midmio_fixedSE_wZBL"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.midmio$run(1000000) # Run the MCMC

priorOU.latemio <- make.prior(latemio.tree.pruned, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 2, kmax = 1000), dsb = list(bmax = Inf, prob = latemio.tree.pruned$edge.length), 
                                 dtheta = list(mean = mean(tl.full), sd = 1.5 * sd(tl.full))),
                      plot.prior = FALSE
                      )

mcmcOU.latemio <- bayou.makeMCMC(latemio.tree.pruned, tl.latemio, SE = se.latemio, prior = priorOU.latemio, file.dir = here::here("output/bayou/adapland/latemiocene/bayou_latemio_fixedSE_wZBL"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.latemio$run(1000000) # Run the MCMC

priorOU.plio <- make.prior(plio.tree.pruned, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 2, kmax = 1000), dsb = list(bmax = Inf, prob = plio.tree.pruned$edge.length), 
                                 dtheta = list(mean = mean(tl.full), sd = 1.5 * sd(tl.full))),
                      plot.prior = FALSE
                      )

mcmcOU.plio <- bayou.makeMCMC(plio.tree.pruned, tl.plio, SE = se.plio, prior = priorOU.plio, file.dir = here::here("output/bayou/adapland/pliocene/bayou_plio_fixedSE_wZBL"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.plio$run(1000000) # Run the MCMC

#######

plot.shifts.coded(mcmcOU.oligo, 0.3, 2, 0.1, tiplab = TRUE, ptp.thresh = 5)
plot.shifts.coded(mcmcOU.earlymio, 0.3, 2, 0.1, tiplab = TRUE, ptp.thresh = 5)
plot.shifts.coded(mcmcOU.midmio, 0.3, 2, 0.1, tiplab = TRUE, ptp.thresh = 5)
plot.shifts.coded(mcmcOU.latemio, 0.3, 2, 0.07, tiplab = TRUE, ptp.thresh = 2)
plot.shifts.coded(mcmcOU.plio, 0.3, 2, 0.08, tiplab = TRUE, ptp.thresh = 2)

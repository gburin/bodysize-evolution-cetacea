if(!require("ape")){install.packages("ape")};library("ape")
if(!require("phytools")){install.packages("phytools")};library("phytools")
if(!require("remotes")){install.packages("remotes")};library("remotes")
if(!require("bayou")){install_github("uyedaj/bayou", upgrade = "never")};library("bayou")
if(!require("Rphylopars")){install.packages("Rphylopars")};library("Rphylopars")
## if(!require("here")){install.packages("here")};library("here")
knitr::opts_chunk$set(message=FALSE, warning=FALSE, fig.width = 20, fig.height = 15)

here::i_am("R/bayou_analyses/noarchaeo/bayou_noarchaeo_k50.R")

## slater.tree <- read.tree(here::here("data/safe_crowngroup_noanc.tre"))
slater.tree <- read.tree(here::here("data/final_tree_cetacea.tre"))
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
dist.name <- adist(unique(logdata$species), slater.tree$tip.label, partial = TRUE, ignore.case = TRUE)
min.name <- apply(dist.name, 1, min)

match.names <- c()
for(i in 1:nrow(dist.name)){
    s2.i <- match(min.name[i], dist.name[i,])
    s1.i <- i
    match.names <- rbind(data.frame(s2.i = s2.i, s1.i = s1.i, s2name = slater.tree$tip.label[s2.i], s1name = unique(logdata$species)[s1.i], adist = min.name[i]), match.names)
}

sp.to.change.data <- match.names$s1name[match.names$adist == 1]
for(i in 1:length(sp.to.change.data)){
    logdata$species[which(!is.na(match(logdata$species, sp.to.change.data[i])))] <- match.names$s2name[match.names$s1name == sp.to.change.data[i]]
}

slater.tree.pruned <- drop.tip(slater.tree, slater.tree$tip.label[!(slater.tree$tip.label %in% unique(logdata$species))])

slater.tree.imput <- paleotree::dropZLB(slater.tree.pruned)
slater.tree.noarchaeo <- extract.clade(slater.tree.imput, getMRCA(slater.tree.imput, c("Mystacodon_selenensis", "Delphinus_capensis")))
baleen.tree <- extract.clade(slater.tree.imput, getMRCA(slater.tree.imput, c("Mystacodon_selenensis", "Aetiocetus_polydentatus")))
toothed.tree <- extract.clade(slater.tree.imput, getMRCA(slater.tree.imput, c("Ashleycetus_planicapitis", "Delphinus_capensis")))
extant.tree <- drop.fossil(slater.tree)

logdata.clean <- logdata[logdata$species %in% slater.tree.noarchaeo$tip.label, ]
logdata.clean.baleen <- logdata[logdata$species %in% baleen.tree$tip.label, ]
logdata.clean.toothed <- logdata[logdata$species %in% toothed.tree$tip.label, ]
logdata.extant <- logdata[logdata$species %in% extant.tree$tip.label, ]

## imput.data <- phylopars(logdata.clean, slater.tree.imput, model = "BM", pheno_error = TRUE)
## saveRDS(imput.data, file = here::here("output/imput_data.RDS"))

imput.data <- readRDS(here::here("output/imput_data.RDS"))
#rownames(imput.data$anc_recon)[Ntip(slater.tree.imput)] <- "Tursiops_australis"

tl.bayou <- setNames(imput.data$anc_recon[1:Ntip(slater.tree.noarchaeo)], rownames(imput.data$anc_recon)[1:Ntip(slater.tree.noarchaeo)])

## Calculating standard error for species with multiple measurements
sd.bayou.per.sp <- aggregate(logdata.clean$tl, by = list(logdata.clean$species), FUN = sd, na.rm = TRUE)
n.per.sp <- aggregate(logdata.clean$tl, by = list(logdata.clean$species), FUN = length)

sd.bayou.per.sp$se <- NA
sd.bayou.per.sp$se[which(!is.na(sd.bayou.per.sp[,2]))] <- sd.bayou.per.sp[which(!is.na(sd.bayou.per.sp[,2])), 2]/sqrt(n.per.sp[which(!is.na(sd.bayou.per.sp[,2])),2])

## Calculating pooling variance to use as standard error for species without multiple measurements
pooled.var <- sum((sd.bayou.per.sp[, 2])^2 * (n.per.sp[, 2] - 1), na.rm = TRUE)/(sum(n.per.sp[, 2], na.rm = TRUE) - sum(!is.na(sd.bayou.per.sp[, 2])))

sd.bayou.per.sp$se[is.na(sd.bayou.per.sp$se)] <- sqrt(pooled.var)/sqrt(n.per.sp[which(is.na(sd.bayou.per.sp[,2])),2])

se.bayou <- sd.bayou.per.sp$se[match(names(tl.bayou), sd.bayou.per.sp[, 1])]
names(se.bayou) <- names(tl.bayou)

## Setting priors - Full tree

priorOU.noarchaeo.k50 <- make.prior(slater.tree.noarchaeo, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = 50, kmax = 1000), dsb = list(bmax = Inf, prob = slater.tree.noarchaeo$edge.length), 
                                 dtheta = list(mean = mean(tl.bayou), sd = 1.5 * sd(tl.bayou))),
                      plot.prior = FALSE
                      )
## startpars.noarchaeo.k50 <- priorSim(priorOU.noarchaeo.k50, slater.tree.noarchaeo, plot = FALSE)$pars[[1]]


## Setting up and running MCMC - Full tree

mcmcOU.noarchaeo.k50 <- bayou.makeMCMC(slater.tree.noarchaeo, tl.bayou, SE = se.bayou, prior = priorOU.noarchaeo.k50, file.dir = here::here("output/bayou/noarchaeo/bayou_noarchaeo_k50_fixedSE"), outname = "modelOU_imput_r001", plot.freq = NULL, samp = 1000) # Set up the MCMC

## mcmcOU.noarchaeo.k50$run(1000000) # Run the MCMC

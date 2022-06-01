library("ape")
library("ggplot2")
library("cowplot")
library("ggtree")
library("ggforce")
library("deeptime")
library("paleobioDB")
library("rgbif")
library("bayou")
library("strap")
source(here::here("R/branching.times.with.extinction.R"))
source(here::here("R/plotting_functions.R"))

load(here::here("output/bayou/chains/all_chains.RData"))
load(here::here("output/bayou/chains/all_mcmc.RData"))
load(here::here("output/bayou/chains/all_tl.RData"))
load(here::here("output/bayou/chains/theta_per_branch.RData"))
load(here::here("output/bayou/chains/theta_hpd_per_branch.RData"))

here::i_am("R/fig_phylo.R")

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

## Extracting shift information

chainOU <- set.burnin(chainOU.full.wZBL.k15, 0.3)
capture.output(chainOU.summ <- summary(chainOU))

shifts.full <- chainOU.summ$branch.posteriors
root.theta <- chainOU.summ$statistics[7, 1]
shifts.full$prior <- prior.per.branch(mcmcOU.full.wZBL.k15, 15)
shifts.full <- shifts.full[order(as.numeric(rownames(shifts.full))),]
shifts.full$posterior.to.prior <- shifts.full$pp / shifts.full$prior
shifts.full <- shifts.full[order(shifts.full$pp, decreasing = TRUE), ]
shifts.full.sig <- shifts.full[shifts.full$pp >= 0.1, ]
if(any(shifts.full.sig$posterior.to.prior == Inf)){
    shifts.full.sig$posterior.to.prior[shifts.full.sig$posterior.to.prior == Inf] <- max(shifts.full.sig$posterior.to.prior[shifts.full.sig$posterior.to.prior != Inf]) * 10
}
shifts.full.sig$ptp.scale <- 0.2
shifts.full.sig$ptp.scale[which(shifts.full.sig$posterior.to.prior >= 10)]  <- 0.5
shifts.full.sig$ptp.scale[which(shifts.full.sig$posterior.to.prior >= 100)]  <- 1
if(!is.null(ptp.thresh)){
    shifts.full.sig <- shifts.full.sig[shifts.full.sig$posterior.to.prior >= ptp.thresh,]
}

shifts.nodes <- data.frame(edge = sort(as.integer(rownames(shifts.full.sig))),
                           anc.node = setNames(mcmcOU.full.wZBL.k15$tree$edge[sort(as.integer(rownames(shifts.full.sig))), 2], NULL))
shifts.nodes$parent.shift <- NA
for(i in 1:(length(shifts.nodes$anc.node) - 1)){
    parent.shifts <- intersect(shifts.nodes$edge, match(phylobase::ancestors(as(mcmcOU.full.wZBL.k15$tree, "phylo4"), shifts.nodes$anc.node[i]), mcmcOU.full.wZBL.k15$tree$edge[,2]))
    if(length(parent.shifts) != 0){
        shifts.nodes$parent.shift[i] <- min(parent.shifts)
    } else {
        shifts.nodes$parent.shift[i] <- 0
    }
}
shifts.nodes$parent.shift[is.na(shifts.nodes$parent.shift)] <- 0

root.node <- getMRCA(mcmcOU.full.wZBL.k15$tree, mcmcOU.full.wZBL.k15$tree$tip.label)
shift.colours <- c("#ff000080", "#0000ff80")
shift.pch <- c(24, 25)
shifts.full.sig$direction <- NA
shifts.full.sig$dir.pch <- NA
shifts.full.sig$change.mag <- NA
for(i in 1:(nrow(shifts.full.sig))){
    focal.theta <- shifts.full.sig$magnitude.of.theta2[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]]
    if(shifts.nodes$parent.shift[i] != 0){
        ancest.theta <- shifts.full.sig$magnitude.of.theta2[as.integer(rownames(shifts.full.sig)) == shifts.nodes$parent.shift[i]]
    } else {
        ancest.theta <- root.theta
    }
    change <- ifelse(focal.theta > ancest.theta, 1, 2)
    shifts.full.sig$change.mag[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- max(focal.theta, ancest.theta)/min(focal.theta, ancest.theta)
    shifts.full.sig$dir.pch[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- shift.pch[change]
    shifts.full.sig$direction[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- ifelse(change == 1, rgb(1, 0, 0, shifts.full.sig$ptp.scale[i]), rgb(0, 0, 1, shifts.full.sig$ptp.scale[i]))
}

edge.labs <- rep(FALSE, nrow(mcmcOU.full.wZBL.k15$tree$edge))
edge.labs[match(rownames(mcmcOU.full.wZBL.k15$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(mcmcOU.full.wZBL.k15$tree)$edge))] <- TRUE

p <- ggtree(slater.tree.pruned) %<+% data.frame(edges = edge.labs)

p.vertical <-
    ggtree(slater.tree.pruned, position = position_nudge(x = -max(branching.times.with.extinct(slater.tree.pruned)))) +
    scale_x_continuous(breaks = seq(-50, 0, 10), labels = abs(seq(-50, 0, 10))) +
    scale_y_reverse() +
    coord_geo(xlim = c(-max(branching.times.with.extinct(slater.tree.pruned)), 0), ylim = c(-2, Ntip(slater.tree.pruned)), abbrv = list(TRUE, FALSE, FALSE), dat = list("stages", "epochs", "periods"), pos = list("bottom", "bottom", "bottom"), center_end_labels = TRUE, size = "auto", fittext_args = list(size = 20), neg = TRUE, expand = TRUE, clip = "off", skip = list("Paleocene", "Thanetian")) +
    annotate("rect", xmin = -2.58, xmax = 0, ymin = -2, ymax = Ntip(slater.tree.pruned), fill = "lightgrey", alpha = 0.25) +
    annotate("rect", xmin = -23.03, xmax = -5.333, ymin = -2, ymax = Ntip(slater.tree.pruned), fill = "lightgrey", alpha = 0.25) +
    annotate("rect", xmin = -56, xmax = -33.9, ymin = -2, ymax = Ntip(slater.tree.pruned), fill = "lightgrey", alpha = 0.25) +
    geom_strip("Balaenoptera_brydei", "Balaenoptera_bonaerensis", label = "Balaenidae", color = "red", offset.text = 0.5, barsize = 2) +
    theme_tree2(axis.text.x = element_text(size = 14)) +
    theme(plot.margin = margin(t = 10, b = 10, r = 5, l = 5, unit = "mm"))


edge <- data.frame(slater.tree.pruned$edge, edge_num = 1:nrow(slater.tree.pruned$edge))
colnames(edge) = c("parent", "node", "edge_num")


ggtree(slater.tree.pruned) + geom_label(aes(x = branch, label = c(edge.labs, FALSE)))
    geom_nodelab(aes(subset = (node == 600), label = "BRANCH"), shape = 19, nudge_x = 0, nudge_y = 0, size = 5, col = "#fe4a49")





source(here::here("R/phenogram_density.R"))

## Palette from https://stackoverflow.com/a/46810812

manualcolors<-c('black',
                'forestgreen',
                'red2',
                'orange',
                'cornflowerblue', 
                #'magenta',
                'darkolivegreen4',
                'indianred1',
                'tan4',
                'darkblue', 
                #'mediumorchid1',
                'firebrick4',
                'yellowgreen',
                'lightsalmon',
                'tan3',
                'tan1',
                'darkgray',
                #'wheat4',
                '#DDAD4B',
                'chartreuse', 
                #'seagreen1',
                #'moccasin',
                'mediumvioletred',
                'seagreen',
                'cadetblue1',
                'darkolivegreen1',
                'tan2',
                'tomato3' ,
                '#7CE3D8',
                'gainsboro')


plot.cols <- c("cornflowerblue", #1
               "firebrick4", #2
               "orange", #3
               "indianred1", #4
               "tomato3", #5
               "darkblue", #6
               "tan3", #7
               "yellowgreen", #8
               "darkolivegreen4", #9
               "cadetblue1", #10
               "tan4", #11
               "seagreen", #12
               "red2", #13
               "darkgray", #14
               "tan2", #15
               "chartreuse") #16


par(fig = c(0.07, 0.375, 0.25, 0.65), new = T)

#plot.cols <- manualcolors[sample(1:length(manualcolors), 16)]
phenogram.density.new(mcmcOU.full.wZBL.k15$tree, tl.full.wZBL, burnin = 0.3, chain = chainOU.full.wZBL.k15, colors = plot.cols, pp.cutoff = 0.1, lwd = 3, fsize = 0.01)

## Final figure

pdf(here::here("output/fig2_draft_w_legend.pdf"), width = 180/25, height = 225/25)
#par(xpd = NA, mar = c(4, 4, 4, 4))

th <- max(branching.times.with.extinct(slater.tree.pruned))
slater.tree.pruned$root.time <- th

geoscalePhylo(slater.tree.pruned, units = c("Period", "Epoch"), boxes = "Epoch", edge.col = NA, cex.tip = 0.01, erotate = 0, cex.ts = 1.5, quat.rm = TRUE, cex.age = 2)
alpha.bg <- 50

#segments(x0 = seq(0, th, 5), y0 = -5, y1 = 350, lty = 2, col = "lightgrey")

par(new = TRUE)
plot.shifts.coded.color(mcmcOU.full.wZBL.k15,
                  burnin = 0.3,
                  nshift = 15,
                  pp.cutoff = 0.1,
                  tiplab = FALSE,
                  change.cex = FALSE,
                  lwd = 1.5,
                  edge.col = "#00000080",
                  colour = plot.cols,
                  theta.value = FALSE,
                  shift.numb = TRUE)

legend(x = 0,
       y = 250,
       legend = c("Increse in Theta", "Decrease in Theta"),
       pch = c(24, 25),
       col = "black",
       bg = "white"
       )

legend(x = 2.5,
       y = 200,
       legend = c("Pakicetids",
                  "Protocedits",
                  "Shift 3",
                  "Shift 4",
                  "Shift 5",
                  "Shift 6",
                  "Shift 7",
                  "Shift 8",
                  "Shift 9",
                  "Shift 10",
                  "Shift 11",
                  "Shift 12",
                  "Shift 13",
                  "Shift 14"),
       pch = as.character(1:14),
       col = "black",
       bg = "white"
       )

dev.off()

pdf(here::here("output/phenogram_k15.pdf"), width = 88/12.5, height = 180/12.5)
th <- max(branching.times.with.extinct(slater.tree.pruned))
slater.tree.pruned$root.time <- th
geoscalePhylo(slater.tree.pruned, units = c("Period", "Epoch"), boxes = "Epoch", edge.col = NA, cex.tip = 0.01, erotate = 0, cex.ts = 1.5, quat.rm = TRUE, cex.age = 2)
alpha.bg <- 50
par(new = TRUE)
phenogram.density.new(mcmcOU.full.wZBL.k15$tree, tl.full.wZBL, burnin = 0.3, chain = chainOU.full.wZBL.k15, colors = plot.cols, pp.cutoff = 0.1, lwd = 3, fsize = 0.01, ylab = "", xlab = "")
dev.off()


### Phylo with bars

pdf(here::here("output/phylo_with_bars.pdf"), width = 180/25, height = 225/25);plotTree.wBars(ladderize(mcmcOU.full.wZBL.k15$tree), tl.full.wZBL);dev.off()

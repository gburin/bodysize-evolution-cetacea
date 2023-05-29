plot.bayou <- function(chain, n.shifts, pal, phy, plot.title = "Plot", edge.label = FALSE, circle = TRUE, add = FALSE, plot = TRUE, pal.pack = "rcb", cex.circle = 1){
    post <- chain$sb[round(0.3 * length(chain$sb), 0):length(chain$sb)]
    high.pp <- as.numeric(names(sort(table(unlist(post))))[(length(table(unlist(post))) - n.shifts + 1):length(table(unlist(post)))])
    pp <- sort(table(unlist(post)))/length(post)
    #L <- Lposterior(chain, tree = phy, burnin = 0.3)
    pars <- list()
    pars$sb <- as.numeric(high.pp)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb) + 1
    pars$loc <- 0.5 * phy$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb) + 1)
    tr <- ladderize(pars2simmap(pars, phy)$tree)
    regimes <- as.numeric(sapply(tr$maps, function(x){names(x)[length(x)]}))
    branch.cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(pars$ntheta)[regimes[match(rownames(ladderize(phy)$edge), rownames(phy$edge))]]
    ## if(n.shifts + 1 <= 11){
    ##     cols <- RColorBrewer::brewer.pal(n.shifts + 1, pal)
    ## } else {
    ##     cols <- RColorBrewer::brewer.pal(11, pal)
    ## }
    if(pal.pack == "rcb"){
        cols <- colorRampPalette(RColorBrewer::brewer.pal(11, pal))(pars$ntheta)
    } else {
        cols <- colorRampPalette(MetBrewer::met.brewer(name = pal, n = 11))(pars$ntheta)
    }
    if(plot == TRUE){
        plot(ladderize(tr), edge.col = branch.cols, edge.width = 2, cex = 0.5, mar = c(5.1, 1.1, 2.1, 0))
        title(plot.title)
        if(circle == TRUE){
            edgelabels(edge = as.numeric(rownames(tr$edge)[high.pp]), pie = rep(1, length(high.pp)), frame = "none", bg = NULL, piecol = "#FF000075", cex = cex.circle)
        }
        if(edge.label == TRUE){
            edgelabels(edge = as.numeric(rownames(tr$edge)[high.pp]), text = round(pp[(length(pp) - length(high.pp)):length(pp)], 3), frame = "none", bg = NULL, cex = 1.2, adj = c(0.5, -1))
        }
    }
    return(ladderize(tr))
}



plotTree.barplot.match <- function (tree, x, args.plotTree = list(), args.barplot = list(), chain, n.shifts, pal = "Set3", phy = mcmcOU$tree, edge.label = TRUE, circle = TRUE, prop = c(3, 1), plot.title = NULL, ...) 
{
    if(n.shifts + 1 <= 11){
    colors <- setNames(RColorBrewer::brewer.pal(n.shifts + 1, pal), 1:(n.shifts + 1))
} else {
    colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(11, pal))(n.shifts + 1), 1:(n.shifts + 1))
}
    tree.map <- plot.bayou(chain = chain, n.shifts = n.shifts, pal = pal, phy = phy, edge.label = TRUE, circle = TRUE, plot = FALSE)
    ss <- getStates(tree.map, "tips")
    barcols <- setNames(sapply(ss, function(x, y){y[which(names(y) == x)]}, y = colors), names(ss))
    if (hasArg(add)) 
        add <- list(...)$add
    else add <- FALSE
    if (hasArg(args.axis)) 
        args.axis <- list(...)$args.axis
    else args.axis <- list()
    args.axis$side <- 1
    cw <- reorder(tree.map)
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        x <- x[cw$tip.label, ]
        x <- t(x)
    }
    args.barplot$height <- if (is.matrix(x)) 
        x[, cw$tip.label]
    else x[cw$tip.label]
    args.barplot$plot <- FALSE
    args.barplot$horiz <- TRUE
    args.barplot$axes <- FALSE
    args.barplot$names.arg <- rep("", Ntip(cw))
    args.barplot$col = barcols
    args.barplot$border = barcols
    if (is.null(args.barplot$beside)) 
        args.barplot$beside <- FALSE
    if (is.null(args.barplot$space)) 
        args.barplot$space <- if (args.barplot$beside) 
            c(0, 1)
        else 0.7
    if (is.null(args.barplot$mar)) 
        args.barplot$mar <- c(5.1, 0, 2.1, 1.1)
    else args.barplot$mar[2] <- 0.1
    obj <- as.matrix(do.call(barplot, args.barplot))
    if (dim(obj)[2] == 1) 
        obj <- t(obj)
    args.plotTree$tips <- setNames(colMeans(obj), cw$tip.label)
    args.barplot$plot <- TRUE
    args.barplot$ylim <- range(args.plotTree$tips)
    args.plotTree$tree <- cw
    args.barplot$xlim <- c(3, 4.5)
    if (is.null(args.plotTree$mar)) 
        args.plotTree$mar <- c(5.1, 1.1, 2.1, 0)
    else {
        args.plotTree$mar[4] <- 0.1
    }
    if (args.plotTree$mar[1] != args.barplot$mar[1]) 
        args.plotTree$mar[1] <- args.barplot$mar[1]
    if (args.plotTree$mar[3] != args.barplot$mar[3]) 
        args.plotTree$mar[3] <- args.barplot$mar[3]
    if (is.null(args.plotTree$ftype)) 
        args.plotTree$ftype <- "i"
    if (is.null(args.plotTree$lwd)) 
        args.plotTree$lwd <- 1
    if (!add) 
        #par(mfrow = c(1, 2))
    layout(matrix(c(1, 2), ncol = 2), widths = prop, heights = c(1, 1))
    #do.call(plotTree, args.plotTree)
    plot.bayou(chain = chain, n.shifts = n.shifts, pal = pal, phy = phy, edge.label = TRUE, circle = TRUE, plot.title = ifelse(is.null(plot.title), paste0("\n\n", n.shifts, " Shifts - No Zero Branch Lengths"), plot.title))
    if (!is.null(args.plotTree$plot) && args.plotTree$plot == 
        FALSE) 
        par(new = TRUE)
    par(mar = args.barplot$mar)
    obj <- do.call(barplot, args.barplot)
    if (!is.null(args.barplot$xlab)) 
        args.axis$xlab <- args.barplot$xlab
    else args.axis$xlab <- "x"
    do.call(axis, args.axis)
    invisible(obj)
}


prior.sim <- function(n, k, tree, trait){
    priorOU <- make.prior(tree, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                 dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.01), dsig2 = list(scale = 0.01),
                                 dk = list(lambda = k, kmax = 2 * Ntip(tree) - 2), dsb = list(bmax = 1, prob = 1), 
                                 dtheta = list(mean = mean(trait), sd = 1.5 * sd(trait))),
                      plot.prior = FALSE
                      )
    startpars <- priorSim(priorOU, slater.tree.imput, plot = FALSE)$pars[[1]]
    return(priorOU(startpars))
}



plotTree.barplot.match2 <- function (tree, x, args.plotTree = list(), args.barplot = list(), chain, pal.pack = "rcb", pal = "Set3", pp.cutoff = 0.3, phy = mcmcOU$tree, edge.label = TRUE, circle = TRUE, prop = c(3, 1), plot.title = NULL, tiplab = FALSE, ...) 
{
    L <- Lposterior(chain, phy)
    L <- L[order(L[,1], decreasing = TRUE),]
    n.shifts <- sum(L[,1] > pp.cutoff)
##     if(n.shifts + 1 <= 11){
##     colors <- setNames(RColorBrewer::brewer.pal(n.shifts + 1, pal), 1:(n.shifts + 1))
## } else {
##     colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(11, pal))(n.shifts + 1), 1:(n.shifts + 1))
    ## }
    if(pal.pack == "rcb"){
        colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(11, pal))(n.shifts + 1), 1:(n.shifts + 1))
    } else {
        colors <- setNames(colorRampPalette(MetBrewer::met.brewer(name = pal, n = 11))(n.shifts + 1), 1:(n.shifts + 1))
        }
    tree.map <- plot.bayou(chain = chain, n.shifts = n.shifts, pal = pal, pal.pack = pal.pack, phy = phy, edge.label = TRUE, circle = TRUE, plot = FALSE)
    ss <- getStates(tree.map, "tips")
    barcols <- setNames(sapply(ss, function(x, y){y[which(names(y) == x)]}, y = colors), names(ss))
    if (hasArg(add)) 
        add <- list(...)$add
    else add <- FALSE
    if (hasArg(args.axis)) 
        args.axis <- list(...)$args.axis
    else args.axis <- list()
    args.axis$side <- 1
    cw <- reorder(tree.map)
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        x <- x[cw$tip.label, ]
        x <- t(x)
    }
    args.barplot$height <- if (is.matrix(x)) 
        x[, cw$tip.label]
    else x[cw$tip.label]
    args.barplot$plot <- FALSE
    args.barplot$horiz <- TRUE
    args.barplot$axes <- FALSE
    args.barplot$names.arg <- rep("", Ntip(cw))
    args.barplot$col = barcols
    args.barplot$border = barcols
    if (is.null(args.barplot$beside)) 
        args.barplot$beside <- FALSE
    if (is.null(args.barplot$space)) 
        args.barplot$space <- if (args.barplot$beside) 
            c(0, 1)
        else 0.7
    if (is.null(args.barplot$mar)) 
        args.barplot$mar <- c(5.1, 0, 2.1, 1.1)
    else args.barplot$mar[2] <- 0.1
    obj <- as.matrix(do.call(barplot, args.barplot))
    if (dim(obj)[2] == 1) 
        obj <- t(obj)
    args.plotTree$tips <- setNames(colMeans(obj), cw$tip.label)
    args.barplot$plot <- TRUE
    args.barplot$ylim <- range(args.plotTree$tips)
    args.plotTree$tree <- cw
    args.barplot$xlim <- c(3, 4.5)
    if (is.null(args.plotTree$mar)) 
        args.plotTree$mar <- c(5.1, 1.1, 2.1, 0)
    else {
        args.plotTree$mar[4] <- 0.1
    }
    if (args.plotTree$mar[1] != args.barplot$mar[1]) 
        args.plotTree$mar[1] <- args.barplot$mar[1]
    if (args.plotTree$mar[3] != args.barplot$mar[3]) 
        args.plotTree$mar[3] <- args.barplot$mar[3]
    if (is.null(args.plotTree$ftype)) 
        args.plotTree$ftype <- "i"
    if (is.null(args.plotTree$lwd)) 
        args.plotTree$lwd <- 1
    if (!add) 
        #par(mfrow = c(1, 2))
    layout(matrix(c(1, 2), ncol = 2), widths = prop, heights = c(1, 1))
    do.call(plotTree, args.plotTree)
    #plotSimmap.mcmc(chain = chain, pp.cutoff = pp.cutoff, edge.type = "none", circles = circle, lwd = 2.5, show.tip.label = tiplab)
    if (!is.null(args.plotTree$plot) && args.plotTree$plot == 
        FALSE) 
        par(new = TRUE)
    par(mar = args.barplot$mar)
    obj <- do.call(barplot, args.barplot)
    if (!is.null(args.barplot$xlab)) 
        args.axis$xlab <- args.barplot$xlab
    else args.axis$xlab <- "x"
    do.call(axis, args.axis)
    invisible(obj)
}





### Colour-coding shifts based on the direction of change

plot.shifts.coded <- function(chain, burnin, nshift, pp.cutoff, type = "phylogram", tiplab = FALSE, cex.lab = 0.5, change.cex = FALSE, lwd = 2, edge.col = "grey", ptp.thresh = NULL, time.axis = FALSE, theta.value = FALSE){
    prior.per.branch <- function(chain, k){
    if(any(chain$tree$edge.length == 0)){
        chain$tree$edge.length <- chain$tree$edge.length + 0.01
    }
    prior <- (1 - (chain$tree$edge.length / sum(chain$tree$edge.length)))^k
    return(1 - prior)
    }
    
    chainOU <- chain$load()
    chainOU <- set.burnin(chainOU, burnin)
    capture.output(chainOU.summ <- summary(chainOU))

    shifts.full <- chainOU.summ$branch.posteriors
    root.theta <- chainOU.summ$statistics[7, 1]
    shifts.full$prior <- prior.per.branch(chain, nshift)
    shifts.full <- shifts.full[order(as.numeric(rownames(shifts.full))),]
    shifts.full$posterior.to.prior <- shifts.full$pp / shifts.full$prior
    shifts.full <- shifts.full[order(shifts.full$pp, decreasing = TRUE), ]
    shifts.full.sig <- shifts.full[shifts.full$pp >= pp.cutoff, ]
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
                               anc.node = setNames(chain$tree$edge[sort(as.integer(rownames(shifts.full.sig))), 2], NULL))
    shifts.nodes$parent.shift <- NA
    for(i in 1:(length(shifts.nodes$anc.node) - 1)){
        parent.shifts <- intersect(shifts.nodes$edge, match(phylobase::ancestors(as(chain$tree, "phylo4"), shifts.nodes$anc.node[i]), chain$tree$edge[,2]))
        if(length(parent.shifts) != 0){
            shifts.nodes$parent.shift[i] <- min(parent.shifts)
        } else {
            shifts.nodes$parent.shift[i] <- 0
        }
    }
    shifts.nodes$parent.shift[is.na(shifts.nodes$parent.shift)] <- 0

    root.node <- getMRCA(chain$tree, chain$tree$tip.label)
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
        if(change.cex == TRUE){
            shifts.full.sig$direction[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- ifelse(change == 1, rgb(1, 0, 0, shifts.full.sig$ptp.scale[i]), rgb(0, 0, 1, shifts.full.sig$ptp.scale[i]))
        } else {
            shifts.full.sig$direction[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- ifelse(change == 1, rgb(1, 0, 0, 1), rgb(0, 0, 1, 1))
        }
    }

    plot(ladderize(chain$tree), type = type, edge.width = lwd, edge.col = edge.col, cex = cex.lab, show.tip.label = tiplab)
    if(type == "phylogram"){
        if(time.axis == TRUE){
            axisPhylo()
        }
    }
    if(change.cex == TRUE){
        edgelabels(edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)), pch = shifts.full.sig$dir.pch, col = shifts.full.sig$direction, frame = "none", bg = shifts.full.sig$direction, cex = exp(shifts.full.sig$change.mag * 2)/3)
    } else {
        edgelabels(edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)), pch = shifts.full.sig$dir.pch, col = gsub("80", "ff", shifts.full.sig$direction), frame = "none", bg = gsub("80", "ff", shifts.full.sig$direction), cex = 2)
    }
    if(theta.value == TRUE){
        edgelabels(text = round(shifts.full.sig$magnitude.of.theta2, 4), edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)), frame = "none", bg = NA, adj = c(0.5, 1.5), cex = 0.75)
        nodelabels(text = round(root.theta, 4), node = getMRCA(chain$tree, chain$tree$tip.label), frame = "none", bg = NA, adj = c(0.75, 0.5), cex = 0.75)
    }
}


        

### Plotting branches coloured by theta values
plot.branch.theta <- function(tree, theta, hpd = NULL, cex.lab = 0.75, cex.theta = 0.75){
    gr <- .bincode(theta, seq(min(theta), max(theta), length.out = length(theta)), include.lowest = TRUE)
    col <- colorRampPalette(heat.colors(3, rev = TRUE)[1:3])(length(theta))[gr]
    if(!is.null(hpd)){
        edgelab <- paste0(round(theta, 2), " [", round(hpd[1,], 2), " - ", round(hpd[2,], 2), "]")
        plot(ladderize(tree$tree), edge.col = col, edge.width = 2, cex = cex.lab)
        edgelabels(text = edgelab, frame = "none", adj = c(0.5, -0.5), cex = cex.theta)
    } else {
        plot(ladderize(tree$tree), edge.col = col, edge.width = 2, cex = cex.lab)
        edgelabels(text = round(theta, 2), frame = "none", adj = c(0.5, -0.5), cex = cex.theta)
    }
}


### Extracting theta median and hpd from bayou
get_branch_par <- function(i, tree, dt, chain, par, stat = "median"){
    branchRegime <- function (branch, abranches, chain, parameter, seqx, summary = FALSE, stat = stat){
        ancs <- c(branch, abranches[[branch]])
        ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == 
                                                                                  ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
        ancshifts <- sapply(ancshifts, function(x) ifelse(length(x) == 0, 1, x))
        ests <- sapply(1:length(ancshifts), function(x) chain[[parameter]][[seqx[x]]][ancshifts[x]])
        res <- cbind(ests)
        if (summary) {
            if(stat == "median"){
                return(apply(res, 2, stats::median))
            } else if(stat == "hpd"){
                return(apply(res, 2, stats::quantile, probs = c(0.025, 0.975)))
            } else if(stat == "dist"){
                return(res)
            }
        }
        else {
            return(res)
        }
    }
    
    inter_mat <- matrix(nrow = length(chain$gen), ncol = length(tree$edge.length))
    variable <- par
    cache <- bayou:::.prepare.ou.univariate(tree, dt)
    abranches <- lapply(1:nrow(tree$edge), bayou:::.ancestorBranches, cache = cache)
    seq1 <- i #The set of samples you want to use to calculate the median value
    suppressWarnings(sapply(1:nrow(tree$edge), function(x){branchRegime(x, abranches, chain, variable, seq1, summary = TRUE, stat = stat)}))
}



prior.per.branch <- function(chain, k){
    if(any(chain$tree$edge.length == 0)){
        chain$tree$edge.length <- chain$tree$edge.length + 0.01
    }
    prior <- (1 - (chain$tree$edge.length / sum(chain$tree$edge.length)))^k
    return(1 - prior)
}

## Modified version of the bayou function for plotting the phenogram, with a small correction to deal with user-defined colors
phenogram.density.new <- 
    function (tree, dat, burnin = 0, chain, colors = NULL, pp.cutoff = NULL, 
              K = NULL, ...) 
{
    tree <- reorder(tree, "postorder")
    dat <- dat[tree$tip.label]
    postburn <- round(length(chain$gen) * burnin, 0):length(chain$gen)
    chain2 <- lapply(chain, function(x) x[postburn])
    theta <- chain2$theta
    no.theta <- lapply(theta, length)
    min.theta <- min(unlist(theta))
    max.theta <- max(unlist(theta))
    if (is.null(K)) {
        K <- list(unique(unlist(no.theta)))
    }
    if (!is.null(pp.cutoff)) {
        L <- Lposterior(chain2, tree)
        pp <- L$pp
        pars <- list()
        pars$sb <- which(pp > pp.cutoff)
        pars$k <- length(pars$sb)
        pars$ntheta <- length(pars$sb) + 1
        pars$loc <- L$rel.location[pars$sb] * tree$edge.length[pars$sb]
        pars$t2 <- 2:(length(pars$sb) + 1)
        if (length(pars$sb) > 0) {
            tr <- pars2simmap(pars, tree)
            tree <- tr$tree
            colors <- colors
            names(colors) <- 1:length(colors)
        }
        else {
            tr <- tree
            colors <- 1
            names(colors) <- 1
        }
    }
    if (is.null(colors)) {
        ntheta <- length(unique(names(unlist(tree$maps))))
        colors <- rainbow(ntheta)
        names(colors) <- 1:ntheta
    }
    nH <- max(nodeHeights(tree))
    #plot(c(0, nH + 0.3 * nH), c(min(dat) - 0.25, max(dat) + 0.25), 
    #     type = "n", xlab = "Time", ylab = "Phenotype")
    phenogram(tree, dat, colors = colors, spread.labels = FALSE, xaxt = "n",
              ...)
    dens.theta <- lapply(1:length(K), function(x) density(unlist(theta[no.theta %in% 
                                                                       K[[x]]])))
    ## tmp <- sapply(1:length(dens.theta), function(Q) {
    ##     lines(nH + dens.theta[[Q]]$y * (0.3 * nH)/max(dens.theta[[Q]]$y), 
    ##           dens.theta[[Q]]$x, col = Q + 1)
    ## })
}




plotSimmap.mcmc.lad <- 
    function (chain, burnin = NULL, lwd = 1, edge.type = c("regimes", 
                                                           "theta", "none", "pp"), pal = rainbow, pp.cutoff = 0.3, circles = TRUE, 
              circle.cex.max = 3, circle.col = "red", circle.pch = 21, 
              circle.lwd = 0.75, circle.alpha = 100, pp.labels = FALSE, 
              pp.col = 1, pp.alpha = 255, pp.cex = 0.75, edge.color = 1, 
              parameter.sample = 1000, ...) 
{
    tree <- attributes(chain)$tree
    edge.type <- match.arg(edge.type, c("regimes", "theta", "none", 
                                        "pp"))
    cache <- bayou:::.prepare.ou.univariate(tree, attributes(chain)$dat)
    tree <- cache$phy
    if (is.null(burnin)) 
        burnin = attributes(chain)$burnin
    if (is.null(burnin)) 
        burnin = 0
    if (burnin == 0) 
        postburn <- 1:length(chain$gen)
    else {
        postburn <- round(burnin * length(chain$gen), 0):length(chain$gen)
    }
    L <- Lposterior(chain, tree, burnin = burnin)
    if (!is.null(pp.cutoff)) {
        pp <- L$pp
        pars <- list()
        pars$sb <- which(pp > pp.cutoff)
        pars$k <- length(pars$sb)
        pars$ntheta <- length(pars$sb) + 1
        pars$loc <- L$rel.location[pars$sb] * tree$edge.length[pars$sb]
        pars$t2 <- 2:(length(pars$sb) + 1)
        if (length(pars$sb) > 0) {
            tr <- pars2simmap(pars, tree)$tree
            colors <- NULL
        }
        else {
            tr <- tree
            colors <- NULL
        }
    }
    else {
        tr <- tree
        tr$maps <- lapply(tr$edge.length, function(x) setNames(x, 
                                                               1))
        colors <- setNames(1, 1)
    }
    .colorRamp <- function(trait, .pal, nn) {
        strait <- (trait - min(trait))/(max(trait - min(trait)))
        itrait <- floor(strait * (nn - 1)) + 1
        if (!is.null(.pal)) {
            return(.pal(nn + 1)[itrait])
        }
        else {
            return(itrait)
        }
    }
    if (edge.type == "none") {
        ape::plot.phylo(ladderize(tr), edge.color = edge.color, lwd = lwd, 
                        ...)
    }
    if (edge.type == "regimes") {
        plotRegimes(ladderize(tr), col = colors, lwd = lwd, pal = pal, ...)
    }
    if (edge.type == "theta") {
        plotBranchHeatMap(ladderize(tree), chain, "theta", burnin = burnin, 
                          pal = heat.colors, ...)
    }
    if (edge.type == "pp") {
        ape::plot.phylo(ladderize(tree), edge.color = .colorRamp(L$pp, pal, 
                                                      100), ...)
    }
    if (circles) {
        circle.cexs <- seq(0, circle.cex.max, length.out = 100)[.colorRamp(L$pp, 
                                                                           NULL, 100)]
        edgelabels(pch = circle.pch, lwd = circle.lwd, bg = makeTransparent(circle.col, 
                                                                            circle.alpha), cex = circle.cexs)
    }
    if (pp.labels) {
        edgelabels(round(L$pp, 2), col = makeTransparent(pp.col, 
                                                         pp.alpha), cex = pp.cex, frame = "none")
    }
}



















### Colour-coding shifts based on the direction of change

plot.shifts.coded.color <- function(chain, burnin, nshift, pp.cutoff, type = "phylogram", tiplab = FALSE, cex.lab = 0.5, change.cex = FALSE, lwd = 2, edge.col = "grey", ptp.thresh = NULL, time.axis = FALSE, theta.value = FALSE, colours = NULL, shift.numb = FALSE){
    prior.per.branch <- function(chain, k){
    if(any(chain$tree$edge.length == 0)){
        chain$tree$edge.length <- chain$tree$edge.length + 0.01
    }
    prior <- (1 - (chain$tree$edge.length / sum(chain$tree$edge.length)))^k
    return(1 - prior)
    }
    chainOU <- chain$load()
    chainOU <- set.burnin(chainOU, burnin)
    capture.output(chainOU.summ <- summary(chainOU))
    shifts.full <- chainOU.summ$branch.posteriors
    root.theta <- chainOU.summ$statistics[7, 1]
    shifts.full$prior <- prior.per.branch(chain, nshift)
    shifts.full <- shifts.full[order(as.numeric(rownames(shifts.full))),]
    shifts.full$posterior.to.prior <- shifts.full$pp / shifts.full$prior
    shifts.full <- shifts.full[order(shifts.full$pp, decreasing = TRUE), ]
    shifts.full.sig <- shifts.full[shifts.full$pp >= pp.cutoff, ]
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
                               anc.node = setNames(chain$tree$edge[sort(as.integer(rownames(shifts.full.sig))), 2], NULL))
    shifts.nodes$parent.shift <- NA
    for(i in 1:(length(shifts.nodes$anc.node) - 1)){
        parent.shifts <- intersect(shifts.nodes$edge, match(phylobase::ancestors(as(chain$tree, "phylo4"), shifts.nodes$anc.node[i]), chain$tree$edge[,2]))
        if(length(parent.shifts) != 0){
            shifts.nodes$parent.shift[i] <- min(parent.shifts)
        } else {
            shifts.nodes$parent.shift[i] <- 0
        }
    }
    shifts.nodes$parent.shift[is.na(shifts.nodes$parent.shift)] <- 0
    root.node <- getMRCA(chain$tree, chain$tree$tip.label)
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
        if(change.cex == TRUE){
            shifts.full.sig$direction[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- ifelse(change == 1, rgb(1, 0, 0, shifts.full.sig$ptp.scale[i]), rgb(0, 0, 1, shifts.full.sig$ptp.scale[i]))
        } else {
            shifts.full.sig$direction[as.integer(rownames(shifts.full.sig)) == shifts.nodes$edge[i]] <- ifelse(change == 1, rgb(1, 0, 0, 1), rgb(0, 0, 1, 1))
        }
    }
    post <- chainOU$sb[round(0.3 * length(chainOU$sb), 0):length(chainOU$sb)]
    high.pp <- rownames(shifts.full.sig)
    pp <- sort(table(unlist(post)))/length(post)
    pars <- list()
    pars$sb <- as.numeric(high.pp)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb) + 1
    pars$loc <- 0.5 * chain$tree$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb) + 1)
    tr <- pars2simmap(pars, chain$tree)$tree
    regimes <- as.numeric(sapply(tr$maps, function(x){names(x)[length(x)]}))
    branch.cols <- colours[c(1,#1
                             3,#2
                             14,#3
                             6,#4
                             12,#5
                             10,#6
                             15,#7
                             8,#8
                             13,#9
                             9,#10
                             7,#11
                             11,#12
                             2,#13
                             4,#14
                             5)][regimes[match(rownames(ladderize(chain$tree)$edge), rownames(chain$tree$edge))]]
    plot(ladderize(chain$tree), edge.col = branch.cols, edge.width = lwd, show.tip.label = tiplab)
    if(type == "phylogram"){
        if(time.axis == TRUE){
            axisPhylo()
        }
    }
    if(change.cex == TRUE){
        edgelabels(edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)),
                   pch = shifts.full.sig$dir.pch,
                   col = shifts.full.sig$direction,
                   frame = "none",
                   bg = shifts.full.sig$direction,
                   cex = exp(shifts.full.sig$change.mag * 2)/3)
    } else {
        edgelabels(edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)),
                   pch = shifts.full.sig$dir.pch,
                   col = "black",
                   frame = "none",
                   bg = "white",
                   cex = 2.5)
    }
    if(shift.numb == TRUE){
        edgelabels(text = 1:nrow(shifts.full.sig), edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)), frame = "n", bg = NA, cex = 1.5, adj = c(0.5, -0.75))
    }
    if(theta.value == TRUE){
        edgelabels(text = round(shifts.full.sig$magnitude.of.theta2, 4), edge = match(rownames(chain$tree$edge)[as.numeric(rownames(shifts.full.sig))], rownames(ladderize(chain$tree)$edge)), frame = "none", bg = NA, adj = c(0.5, 1.5), cex = 0.75)
        nodelabels(text = round(root.theta, 4), node = getMRCA(chain$tree, chain$tree$tip.label), frame = "none", bg = NA, adj = c(0.75, 0.5), cex = 0.75)
    }
}

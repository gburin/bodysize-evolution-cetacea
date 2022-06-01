## Script to get the per branch estimates for the optima of body size.
## This will run in the servers for speed.

library( bayou )
library( ape )
library( geiger )

## dt <- read.csv(file = "../../../data/Proboscidean_tuks_evolution_31032021_trait_data.csv")
dt <- read.csv(file = "Proboscidean_tuks_evolution_31032021_trait_data.csv")
body_mass <- setNames(object = dt$Body.mass..kg., nm = dt$Species..name.match.phylogeny.)
males <- setNames(object = dt$Sex, nm = dt$Species..name.match.phylogeny.)
is_male <- males == "M" & !is.na( males )
male_body_mass <- body_mass[is_male]
spp <- unique( names( male_body_mass ) )
avg_male_body_mass <- vector(mode = "numeric", length = length(spp))
for( i in 1:length(spp) ){
    avg_male_body_mass[i] <- mean( male_body_mass[names(male_body_mass) == spp[i]] )
}
names( avg_male_body_mass ) <- spp
avg_male_body_mass <- log( avg_male_body_mass )
## tree <- read.tree( "../../../phylogeny/dating_phylogeny/full_elepha3_simple_dated.tre" )
tree <- read.tree( "full_elepha3_simple_dated.tre" )
tree_1 <- reorder(tree[[1]], "postorder")
tree_1_bin <- multi2di(phy = tree_1)
tree_1_bin$edge.length[tree_1_bin$edge.length == 0.0] <- min(tree_1_bin$edge.length[tree_1_bin$edge.length > 0.0]) / 100
dt_list <- treedata(phy = tree_1_bin, data = avg_male_body_mass, sort = TRUE)

chainOU <- readRDS(file = "bayou_mass_post.rds")

get_branch_par <- function(i, tree, dt, chain, par, stat = "median"){
    .branchRegime <- function (branch, abranches, chain, parameter, seqx, summary = FALSE, stat = stat){
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
    suppressWarnings(sapply(1:nrow(tree$edge), function(x) bayou:::.branchRegime(x, abranches, chain, variable, seq1, summary = TRUE, stat = stat)))
}

library( parallel )
body_mass_list <- mclapply(1:length(chainOU$gen), function(x) get_branch_par(i = x, tree = dt_list$phy, dt = dt_list$data[,1], chain = chainOU, par = "theta"), mc.cores = 50)
saveRDS(object = body_mass_list, file = "body_mass_per_branch.rds")

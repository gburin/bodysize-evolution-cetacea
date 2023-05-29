library("ape")
library("phytools")

here::i_am("R/fabric_model.R")

load(here::here("output/bayou/chains/all_tl.RData"))

slater.tree <- ladderize(read.tree(here::here("data/final_tree_cetacea.tre")))
slater.tree.pruned <- drop.tip(slater.tree, slater.tree$tip.label[!(slater.tree$tip.label %in% names(tl.full.wZBL))])

d <- read.csv(here::here("output/fabric_model/out.csv"))
head(d)
signodes <- d$Md5.Sum[which(d$No.Sig.Beta > 0)]
x <- read.csv(here::here("output/fabric_model/cetaceans.txt.VarRates.txt.csv"))

write.csv(x[match(signodes, x$Md5.Sum),], here::here("output/fabric_model/sigclades.csv"))

fabric.summary <- read.csv(here::here("output/fabric_model/sigclades.csv"))


beta.change <- data.frame(node = vector(mode = "numeric", length = nrow(fabric.summary)),
                          change = vector(mode = "numeric", length = nrow(fabric.summary))
                          )

for(i in 1:nrow(fabric.summary)){
    if(fabric.summary[i, "No.Taxa"] == 1){
        beta.change[i, "node"] <- match(fabric.summary[i, 43:(42 + fabric.summary[i, "No.Taxa"])], slater.tree.pruned$tip.label)
        beta.change[i, "change"] <- ifelse(fabric.summary[i, "Mean..Beta...BL."] > 0, 24, 25)
        } else {
            beta.change[i, "node"] <- getMRCA(slater.tree.pruned, tip = na.omit(unlist(fabric.summary[i, 43:(42 + fabric.summary[i, "No.Taxa"])])))
            beta.change[i, "change"] <- ifelse(fabric.summary[i, "Mean..Beta...BL."] > 0, 24, 25)
        }
}

beta.change[, "col"] <- ifelse(beta.change[, "change"] == 24, "red", "blue")

pdf(width = 20, height = 40)
plot(slater.tree.pruned, show.tip.label = TRUE, cex = 0.5);axisPhylo()
#nodelabels(node = beta.change[, "node"], pch = beta.change[, "change"], bg = "black")
edgelabels(edge = match(beta.change[, "node"], slater.tree.pruned$edge[,2]), pch = beta.change[, "change"], bg = beta.change[, "col"], col = beta.change[, "col"], cex = 2)
dev.off()

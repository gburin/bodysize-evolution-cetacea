############################################
####                                    ####
####    Daniel Varajao de Latorre       ####
####                                    ####
####            Função                  ####
####                                    ####
############################################

####    BRANCHING.TIMES COMPATIVEL COM EXTINÇÃO     ####
branching.times.with.extinct<-function(phy)
{
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    xx <- numeric(phy$Nnode)
    interns <- which(phy$edge[, 2] > n)
    for (i in interns)
        xx[phy$edge[i, 2] - n] <- xx[phy$edge[i,1] - n] + phy$edge.length[i]
    ##parte nova, novo algoritimo para calcular qual o tempo maximo da arvore
    ntip<-Ntip(phy)
    sum<-numeric(ntip)
    index<-numeric()
    x<-numeric()
    for(i in 1:ntip)
    {
        node<-i
        while(node!=ntip+1)
        {
            index<-which(phy$edge[,2]==node)
            sum[i]<-sum[i] + phy$edge.length[index]
            node<-phy$edge[index,1]
        }
    }
    depth <- max(sum)
    ##volta a ser usado o codigo da função original
    xx <- depth - xx
    names(xx) <- if (is.null(phy$node.label)) 
        (n + 1):(n + phy$Nnode)
    else phy$node.label
    xx
    
}
library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)

library(ggradar)
library(dplyr)
library(scales)
library(tibble)

method = 'gamma'
dir = 'C:/Liang/Code/Pro3/data1_20181220/Mevent'

# Compute the colless values and gamma values for a given phylo tree.
foo <- function(x, metric = "colless") {
  if (metric == "colless") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    colless(xx, "yule")  # calculate colless' metric
  } else if (metric == "gamma") {
    gammaStat(x)
  } else stop("metric should be one of colless or gamma")
}


colless_alldf = data.frame()
count = 0
for(i in c(1:6)){
  for(j in c(1:6)){
    count = count+1
    # numtrees <- 1  # lets simulate 1000 trees
    # trees <- pbtree(n = 50, nsim = numtrees, ape = F)  # simulate 500
    comb=paste0(i,j)
    file =  paste0(dir,comb,'.Rdata')
    # file = paste("/home/p274981/cluster/",n,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    source(file = file)
    L = event2L(result = events)
    phy = DDD::L2phylo(L,dropextinct = F)
    
    # Generate a 'multiphylo' class to contain all trees generated under the same paras
    # tree1 = pbtree(n = 50,nsim = 1)
    # tree2 = pbtree(n = 50, nsim = 1)
    # treelist = list(tree1,tree2)
    # class(treelist) = 'multiphylo'
    
    colless_df <- foo(trees, metric = 'colless')  # calculate metric for each tree
    gamma_df = foo(trees, metric = 'gamma')
    richness = trees$Nnode+1
    
    colless_df = cbind(colless_df,gamma_df,richness,i,j)
    colless_alldf = rbind(colless_alldf,colless_df)
    
  }
}

colnames(colless_alldf) = c('Colless','Gamma','Richness','psi','phi')

tdf = colless_alldf
tdf %>%
  mutate_all(funs(rescale))%>%
  tail(6) -> tdf1

tdf <- cbind(c(1:nrow(tdf)),tdf)

ggradar2(tdf)

windowsFonts(Times=windowsFont("TT Times New Roman"))

tdf1 = cbind(c(1:6),tdf1)
colnames(tdf1) = c('group',colnames(tdf1)[-1])
ggradar2(tdf1,legend.text.size = 8,style = 'straight',
        grid.label.size = 4,
        axis.label.size = 4,polygonfill = FALSE)


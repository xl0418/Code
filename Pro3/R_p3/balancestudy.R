library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
foo <- function(x, metric = "colless") {
  if (metric == "colless") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    colless(xx, "yule")  # calculate colless' metric
  } else if (metric == "gamma") {
    gammaStat(x)
  } else stop("metric should be one of colless or gamma")
}

theme_myblank <- function() {
  stopifnot(require(ggplot2))
  theme_blank <- ggplot2::theme_blank
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), plot.background = element_blank(), 
                 axis.title.x = element_text(colour = NA), axis.title.y = element_blank(), 
                 axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_blank(), 
                 axis.ticks = element_blank())
}

colless_alldf = data.frame()
count = 0
for(i in c(1:5)){
  for(j in c(1:5)){
    count = count+1
    numtrees <- 100  # lets simulate 1000 trees
    trees <- pbtree(n = 50, nsim = numtrees, ape = F)  # simulate 500
    
    # Generate a 'multiphylo' class to contain all trees generated under the same paras
    # tree1 = pbtree(n = 50,nsim = 1)
    # tree2 = pbtree(n = 50, nsim = 1)
    # treelist = list(tree1,tree2)
    # class(treelist) = 'multiphylo'
    
    colless_df <- ldply(trees, foo, metric = "colless")  # calculate metric for each tree
    colless_df = cbind(colless_df,count,i,j)
    colless_alldf = rbind(colless_alldf,colless_df)
   
  }
}
colnames(colless_alldf) = c('colvalue','group','psi','phi')


plothis <- ggplot(colless_alldf, aes(colvalue)) +  # plot histogram of distribution of values
  geom_histogram(binwidth = 0.2) + 
  theme_bw(base_size=18) + 
  scale_x_continuous(limits=c(-3,3), breaks=c(-3,-2,-1,0,1,2,3)) + 
  geom_vline(xintercept = -0.7, colour="red", linetype = "longdash") +
  geom_vline(xintercept = 0.7, colour="red", linetype = "longdash") +
  labs(x=expression(" "),y=expression(" ")) +
  background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
  panel_border()

plothis+facet_grid(vars(psi),vars(phi))

wholeplot = plothis+facet_grid(vars(psi),vars(phi))
y.grob <- textGrob("Frequency", 
                   gp=gpar(fontface="bold", col="blue", fontsize=15), rot=90)

x.grob <- textGrob("Colles index", 
                   gp=gpar(fontface="bold", col="blue", fontsize=15))

#add to plot

grid.arrange(arrangeGrob(wholeplot, left = y.grob, bottom = x.grob))




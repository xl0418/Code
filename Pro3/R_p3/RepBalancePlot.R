library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
method = 'gamma'
dir <- 'C:/Liang/Googlebox/Research/Project3/batchsim_results/'
scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')
jclabel = c(0,0.5,1)
plabel = c(0,1e2,1e4,1e6,1e8,-1)
i_n=1
sce = scenario[i_n]
f.name = sce.short[i_n]
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
for(i in c(1:3)){
  for(j in c(1:6)){
    comb=paste0(i,j)
    multitreefile <- paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/multitree',f.name,comb,'.tre')
    trees <- read.tree(multitreefile)
    colless_df <- ldply(trees, foo, metric = method)  # calculate metric for each tree
    colless_df = cbind(colless_df,count,i,j)
    colless_alldf = rbind(colless_alldf,colless_df)
    
  }
}
colnames(colless_alldf) = c('colvalue','group','psi','phi')


plothis <- ggplot(colless_alldf, aes(colvalue)) +  # plot histogram of distribution of values
  geom_histogram(binwidth = 0.2)+ theme_wsj() +
  # theme_bw(base_size=18) + 
  scale_x_continuous(limits=c(-3,3), breaks=c(-3,-2,-1,0,1,2,3)) + 
  geom_vline(xintercept = -0.7, colour="red", linetype = "longdash") +
  geom_vline(xintercept = 0.7, colour="red", linetype = "longdash") +
  labs(x=expression(" "),y=expression(" ")) 
# background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
# panel_border()

plothis+facet_grid(vars(psi),vars(phi))

wholeplot = plothis+facet_grid(vars(psi),vars(phi))
y.grob <- textGrob("Frequency", 
                   gp=gpar(fontface="bold", col="black", fontsize=13), rot=90)

x.grob <- textGrob(paste(method,"index"), 
                   gp=gpar(fontface="bold", col="black", fontsize=13))

#add to plot

grid.arrange(arrangeGrob(wholeplot, left = y.grob, bottom = x.grob))
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
jclabel = c('0','0.5','1')
plabel = c('0','1e2','1e4','1e6','1e8','Inf')
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
for(i in c(1:3)){
  for(j in c(1:6)){
    comb=paste0(i,j)
    multitreefile <- paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/multitree',f.name,comb,'.tre')
    trees <- read.tree(multitreefile)
    colless_df <- ldply(trees, foo, metric = method)  # calculate metric for each tree
    colless_df = cbind(colless_df,jclabel[i],plabel[j])
    colless_alldf = rbind(colless_alldf,colless_df)
    
  }
}
colnames(colless_alldf) = c('colvalue','psi','phi')

value.min <- min(colless_alldf$colvalue)
value.max <- max(colless_alldf$colvalue)


plothis <- ggplot(colless_alldf, aes(colvalue)) +  # plot histogram of distribution of values
  geom_histogram(binwidth = .5)+ theme_wsj() +
  # theme_bw(base_size=18) + 
  scale_x_continuous(limits=c(3,15), breaks=c(3,6,9,12,15)) + 
  # geom_vline(xintercept = -0.7, colour="red", linetype = "longdash") +
  # geom_vline(xintercept = 0.7, colour="red", linetype = "longdash") +
  labs(x=expression(" "),y=expression(" ")) 
# background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
# panel_border()



wholeplot = plothis+facet_grid(psi~phi)+
            stat_bin(bins=20)  +
  theme(strip.text.x = element_text(size=13),
        strip.text.y = element_text(size=13))

y.grob <- textGrob("Frequency", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), rot=90)

x.grob <- textGrob(paste(method,"index"), 
                   gp=gpar(fontface="bold", col="black", fontsize=10))

#add to plot

grid.arrange(arrangeGrob(wholeplot, left = y.grob, bottom = x.grob))
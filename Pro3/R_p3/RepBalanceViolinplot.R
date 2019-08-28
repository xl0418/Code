library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
method = 'colless'
dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')
jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('0','1e2','1e4','1e6','1e8','Inf')
wholeplot = list()
xlabels <- c(~list(psi==0),~list(psi==0.25),~list(psi==0.5),
             ~list(psi==0.75),~list(psi==1))
subtitles = c('Low distance','Intermediate distance','High distance')
for(i_n in c(1:3)){
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
  
  i.true.order = 1
  
  colless_alldf = data.frame()
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      if(i==1){
        comb=paste0(i,i)
      }else{
        comb=paste0(i,j)
      }    
      multitreefile <- paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/multitree',f.name,comb,'.tre')
      trees <- read.tree(multitreefile)
      colless_df <- ldply(trees, foo, metric = method)  # calculate metric for each tree
      colless_df = cbind(colless_df,jclabel[i.true.order],plabel[j])
      colless_alldf = rbind(colless_alldf,colless_df)
      
    }
    i.true.order=i.true.order+1
    
  }
  colnames(colless_alldf) = c('colvalue','psi','phi')
  
  if(method=='gamma'){
    value.min <- 0
    value.max <- 15
    gap <- 3
  }else{
    value.min <- 5
    value.max <- 40
    gap <- 10
  }
  
  colless_alldf$psi1 <- factor(colless_alldf$psi, labels = c('psi==0','0.25','0.5','0.75','1'))
  colless_alldf$phi1 <- factor(colless_alldf$phi, labels = c('phi==Inf','phi==10^{-2}',
                                                             'phi==10^{-4}','phi==10^{-6}',
                                                             'phi==10^{-8}','phi==0'))
  
  
  plothis <- ggplot(colless_alldf, aes(factor(psi),colvalue)) +  # plot histogram of distribution of values
    geom_violin(fill = "#A8D8B9", colour = "#227D51",width=0.8)+ theme_gdocs() +
    labs(x='',y="") +ylim(c(value.min,value.max))+
    scale_x_discrete(labels=xlabels)+
    ggtitle(subtitles[i_n])+theme(plot.title = element_text(size = 10, face = "bold"))
  
  
  wholeplot[[i_n]] = plothis+facet_grid( rows=vars(phi1),
                                  labeller = label_parsed)+
    theme(strip.text.x = element_text(size=13),
          strip.text.y = element_text(size=13))
  
}

m = matrix(1:3,1,3)

grob1 = textGrob("Values", gp=gpar(fontsize=16),rot = 90)

grob2 = arrangeGrob(grobs = wholeplot, layout_matrix = m,respect=TRUE)
grob4 = textGrob(method, gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")



grid.arrange(grob1,grob2,grob3,grob3,grob4,grob3,ncol = 3, widths = c(1,24,1),heights = c(20,1))



savefilename <- paste0(dir,f.name,'_',method,'.pdf')
ggsave(savefilename,wholeplot,width = 15,height = 10)

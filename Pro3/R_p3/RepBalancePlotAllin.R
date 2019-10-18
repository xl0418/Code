library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
source('C:/Liang/Code/Pro3/R_p3/deltar.R', echo=TRUE)
method = 'deltar'
dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')
jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('0','1e2','1e4','1e6','1e8','Inf')
plotl_est <- list()

for(tens1 in c(1:3)){
  sce = scenario[tens1]
  f.name = sce.short[tens1]
  # Compute the colless values and gamma values for a given phylo tree.
  foo <- function(x, metric = "colless") {
    if (metric == "colless") {
      xx <- as.treeshape(x)  # convert to apTreeshape format
      colless(xx, "yule")  # calculate colless' metric
    } else if (metric == "gamma") {
      gammaStat(x)
    } else if (metric == "deltar") {
      deltar(x)
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
  }else if(method=='deltar'){
    value.min <- -1
    value.max <- 1
    gap <- 0.1
    
  }else{
    value.min <- 5
    value.max <- 40
    gap <- 10
  }
  
  colless_alldf$psi1 <- factor(colless_alldf$psi, labels = c('psi==0','0.25','0.5','0.75','1'))
  colless_alldf$phi1 <- factor(colless_alldf$phi, labels = c('sigma[phi]==0','10^2',
                                                             '10^4','10^6','10^8','Inf'))
  
  

  plotl_est[[tens1]] <-ggplot(colless_alldf, aes(x=phi,y=colvalue)) +  # plot histogram of distribution of values
    geom_boxplot(aes(fill=psi), position=position_dodge(.9))+ theme_gdocs() +
    scale_fill_manual(values=c("#080808", "#4A225D", "#E83015","#F05E1C","#FFC408"))+
    labs(x=method,y="Frequency") +
    theme(legend.position = "none")


}

legend_plot = ggplot(colless_alldf, aes(x=phi,y=colvalue)) +  # plot histogram of distribution of values
  geom_boxplot(aes(fill=psi), position=position_dodge(.9))+ theme_gdocs() +
  scale_fill_manual(values=c("#080808", "#4A225D", "#E83015","#F05E1C","#FFC408"))+
  labs(x=method,y="Frequency") 

mylegend = g_legend(legend_plot)

sigdisp1 <- textGrob(substitute(paste(sigma[disp],'=0.1')), gp=gpar(fontsize=16),rot = 0)
sigdisp2 <- textGrob(substitute(paste(sigma[disp],'=1')), gp=gpar(fontsize=16),rot = 0)
sigdisp3 <- textGrob(substitute(paste(sigma[disp],'=10')), gp=gpar(fontsize=16),rot = 0)

sigjc1 <- textGrob(substitute(paste(sigma[JC],'=0.1')), gp=gpar(fontsize=16),rot = 90)
sigjc2 <- textGrob(substitute(paste(sigma[JC],'=1')), gp=gpar(fontsize=16),rot = 90)
sigjc3 <- textGrob(substitute(paste(sigma[JC],'=10')), gp=gpar(fontsize=16),rot = 90)

m = matrix(1:3,1,3)

grob.sigdisp <- arrangeGrob(sigdisp1, sigdisp2,sigdisp3,ncol = 3)
grob.sigjc <- arrangeGrob(sigjc1,ncol = 1)
grob2 = arrangeGrob(grobs = plotl_est, layout_matrix = m)
grob3 = textGrob("")


grid.arrange(grob3,grob.sigdisp,grob3,grob.sigjc,grob2,mylegend,ncol = 3,
             widths = c(1,16,2.5),heights = c(1,16))









savefilename <- paste0(dir,f.name,'_',method,'.pdf')
ggsave(savefilename,wholeplot,width = 15,height = 10)
# 
# y.grob <- textGrob("Frequency", 
#                    gp=gpar(fontface="bold", col="black", fontsize=10), rot=90)
# 
# x.grob <- textGrob(paste(method,"index"), 
#                    gp=gpar(fontface="bold", col="black", fontsize=10))
# 
# #add to plot
# 
# grid.arrange(arrangeGrob(wholeplot, left = y.grob, bottom = x.grob))

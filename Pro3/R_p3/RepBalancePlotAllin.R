library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
library(phyloTop)
source('C:/Liang/Code/Pro3/R_p3/deltar.R', echo=TRUE)
method = 'colless'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = rev(c('H','M','L'))
scenario = NULL
sce.short.comb.vec = NULL
for(i.letter in sce.short){
  for(j.letter in sce.short){
    sce.folder = paste0('sce',i.letter,j.letter)
    scenario = c(scenario,sce.folder)
    sce.short.comb = paste0(i.letter,j.letter)
    sce.short.comb.vec = c(sce.short.comb.vec,sce.short.comb)
  }
}

# Compute the colless values and gamma values for a given phylo tree.
foo <- function(x, metric = "colless") {
  if (metric == "colless") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    num.tips <- x$Nnode+1
    apTreeshape::colless(xx, "yule")*2/(num.tips-1)/(num.tips-2)  # calculate colless' metric
  } else if (metric == "gamma") {
    gammaStat(x)
  } else if (metric == "deltar") {
    deltar(x)
  } else if (metric == "sackin") {
    sackin.phylo(x)
  } else stop("metric should be one of colless or gamma")
}

jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('1','1e-2','1e-4','1e-6','1e-8','0')
x.ticks.labels = c(
  ~ "1",
  ~ list(10^{-2}),
  ~ list(10^{-4}),
  ~ list(10^{-6}),
  ~ list(10^{-8}),
  ~ "0"
)
plotl_est <- list()

for(tens1 in c(1:9)){
  scefolder = scenario[tens1]
  letter.comb = sce.short.comb.vec[tens1]
  i.true.order = 1
  print(scefolder)
  
  colless_alldf = data.frame()
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb=paste0(i,j)
      
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')
      
      trees <- read.tree(multitreefile)
      colless_df <- ldply(trees, foo, metric = method)  # calculate metric for each tree
      colless_df = cbind(colless_df,jclabel[i.true.order],plabel[j])
      colless_alldf = rbind(colless_alldf,colless_df)
      
    }
    i.true.order=i.true.order+1
    
  }
  colnames(colless_alldf) = c('colvalue','psi','phi')
  
  if(method=='gamma'){
    value.min <- -20
    value.max <- 5
    gap <- 3
  }else if(method=='deltar'){
    value.min <- -1.2
    value.max <- 0.5
    gap <- 0.1
  }else if(method=='sackin'){
    value.min <- -1.2
    value.max <- 40000
    gap <- 0.1
    
  }else{
    value.min <- 0
    value.max <- 0.025
    gap <- 0.005
  }
  
  colless_alldf$psi1 <- factor(colless_alldf$psi, labels = c('psi==0','psi==0.25','psi==0.5','psi==0.75','psi==1'))
  colless_alldf$phi1 <- factor(colless_alldf$phi, labels = c('phi==1','phi==10^-2',
                                                             'phi==10^-4','phi==10^-6','phi==10^-8','phi==0'))
  
  

  plotl_est[[tens1]] <-ggplot(colless_alldf, aes(x=phi,y=colvalue)) +  # plot histogram of distribution of values
    geom_boxplot(aes(fill=psi), position=position_dodge(.9),outlier.shape = NA)+ theme_tufte() +
    scale_fill_manual(values=c("#080808", "#4A225D", "#E83015","#F05E1C","#FFC408"))+
    labs(x="",y="") +
    theme(legend.position = "none",axis.title.x=element_blank(),axis.text.x=element_blank(),
          axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line()
          ) + ylim(value.min, value.max)
  if(tens1 %in% c(1,2,3)){
    plotl_est[[tens1]] <- plotl_est[[tens1]]+
    theme(axis.title.y=element_text(color="black", size=15, face="bold"),
          axis.text.y=element_text(color="black", size=15, face="bold"),
          axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line())
  }
  if(tens1 %in% c(3,6,9)){
    plotl_est[[tens1]] <- plotl_est[[tens1]]+
      theme(axis.title.x=element_text(color="black", size=15, face="bold"),
            axis.text.x=element_text(color="black", size=15, face="bold"),
            axis.ticks.x=element_line(),axis.line.x=element_line(),
            axis.ticks.y=element_line(),axis.line.y=element_line())+
      scale_x_discrete(breaks=plabel,labels = x.ticks.labels)
  }

}

legend_plot = ggplot(colless_alldf, aes(x=phi,y=colvalue)) +  # plot histogram of distribution of values
  geom_boxplot(aes(fill=psi), position=position_dodge(.9))+ theme_gdocs() +
  scale_fill_manual(values=c("#080808", "#4A225D", "#E83015","#F05E1C","#FFC408"))+
  labs(x=method,y="Frequency") +guides(fill=guide_legend(title=expression(psi)))

mylegend = g_legend(legend_plot)

sigdisp1 <- textGrob(substitute(paste(sigma[disp],'=0.1')), gp=gpar(fontsize=16),rot = 0)
sigdisp2 <- textGrob(substitute(paste(sigma[disp],'=1')), gp=gpar(fontsize=16),rot = 0)
sigdisp3 <- textGrob(substitute(paste(sigma[disp],'=10')), gp=gpar(fontsize=16),rot = 0)

sigjc1 <- textGrob(substitute(paste(sigma[JC],'=0.1')), gp=gpar(fontsize=16),rot = 90)
sigjc2 <- textGrob(substitute(paste(sigma[JC],'=1')), gp=gpar(fontsize=16),rot = 90)
sigjc3 <- textGrob(substitute(paste(sigma[JC],'=10')), gp=gpar(fontsize=16),rot = 90)

m = matrix(1:9,3,3)

grob.sigdisp <- arrangeGrob(sigdisp1, sigdisp2,sigdisp3,ncol = 3)
grob.sigjc <- arrangeGrob(sigjc1,sigjc2,sigjc3,ncol = 1)
grob2 = arrangeGrob(grobs = plotl_est, layout_matrix = m)
grob3 = textGrob("")


wholeplot=grid.arrange(grob3,grob.sigdisp,grob3,grob.sigjc,grob2,mylegend,ncol = 3,
             widths = c(1,16,2.5),heights = c(1,16))


dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
savefilename <- paste0(dir_save,method,'_dis.pdf')
ggsave(savefilename,wholeplot,width = 15,height = 10)


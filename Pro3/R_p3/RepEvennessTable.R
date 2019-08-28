library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(vegan)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/1e+07/'
plot_dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
scenario = c('LR','MR','HR')
sce.short = c('L','M','H')
x.breaks = seq(0,18,2)
max.logabund = 12
tas1=list()
base_size <- 6
q.scale=2
i.true.order = 1
for(i_n in c(1:3)){
  evenness.df = NULL
  sce = scenario[i_n]
  f.name = sce.short[i_n]
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      evenness.vec=NULL
      if(i==1){
        comb=paste0(i,i)
      }else{
        comb=paste0(i,j)
      }
      
      for(rep in c(1:100)){
        rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
        Rs = read.csv(rname,header = FALSE)
        evenness = tsallis(Rs,scales=q.scale,norm=TRUE)
        evenness.vec = rbind(evenness.vec, evenness)
      }
      mean.sim = mean(evenness.vec)
      sd.sim = sqrt(var(evenness.vec))
      
      evenness.df = rbind(evenness.df,cbind(mean.sim,sd.sim,i.true.order,j))
    }
    i.true.order=i.true.order+1
  }
  
  colnames(evenness.df) <- c('mean','sd','JCvalue','phyvalue')
  evenness.df <- as.data.frame(evenness.df)
  
  evenness.df$JCvalue1 <- factor(evenness.df$JCvalue, labels = c('psi==0','psi==0.25','psi==0.5',
                                                                 'psi==0.75','psi==1'))
  evenness.df$phyvalue1 <- factor(evenness.df$phyvalue, labels = c('sigma[phi]==0','sigma[phi]==10^2',
                                                                   'sigma[phi]==10^4','sigma[phi]==10^6',
                                                                   'sigma[phi]==10^8','sigma[phi]==Inf'))
  
  xlabels <- c(~list(psi==0),~list(psi==0.25),~list(psi==0.5),
               ~list(psi==0.75),~list(psi==1))
  ylabels <- c(~list(phi==Inf),~list(phi==10^{-2}),~list(phi==10^{-4}),
               ~list(phi==10^{-6}),~list(phi==10^{-8}),~list(phi==0))
  
  
  if(i_n == 1){
    tas1[[i_n]] <- ggplot(evenness.df, aes(JCvalue1, phyvalue1)) + geom_tile(aes(fill = mean),
                                                                             colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")+
      theme_grey(base_size = base_size) + labs(x = "", y = "") + 
      geom_text(aes(label = round(mean, 3)))+
      scale_x_discrete(expand = c(0, 0),labels=xlabels) +scale_y_discrete(expand = c(0, 0),
                                                                          labels=ylabels) + 
      theme(legend.position = "none",axis.ticks = element_blank(), 
            axis.text.x = element_text(size = base_size *2, angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = base_size *2, angle = 0, hjust = 0, colour = "black"))+
      ggtitle("Low distance")+theme(plot.title = element_text(size = 10, face = "bold"))
    
  }
  if(i_n == 2){
    tas1[[i_n]] <- ggplot(evenness.df, aes(JCvalue1, phyvalue1)) + geom_tile(aes(fill = mean),
                                    colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")+
   theme_grey(base_size = base_size) + labs(x = "", y = "") + 
      geom_text(aes(label = round(mean, 3)))+
      scale_x_discrete(expand = c(0, 0),labels=xlabels) +scale_y_discrete(expand = c(0, 0),
                                                                          labels=ylabels) + 
      theme(legend.position = "none",axis.ticks = element_blank(), 
            axis.text.x = element_text(size = base_size *2, angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_blank())+
      ggtitle("Intermediate distance")+theme(plot.title = element_text(size = 10, face = "bold"))
    
  }
  if(i_n == 3){
    tas1[[i_n]] <-   ggplot(evenness.df, aes(JCvalue1, phyvalue1)) + geom_tile(aes(fill = mean),
                                                                               colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")+
      theme_grey(base_size = base_size) + labs(x = "", y = "") + 
      geom_text(aes(label = round(mean, 3)))+
      scale_x_discrete(expand = c(0, 0),labels=xlabels) +scale_y_discrete(expand = c(0, 0),
                                                                          labels=ylabels) + 
      theme(legend.position = "none",axis.ticks = element_blank(), 
            axis.text.x = element_text(size = base_size *2, angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_blank())+
      ggtitle("High distance")+theme(plot.title = element_text(size = 10, face = "bold"))
    
  }
  
}

legend_plot = ggplot(evenness.df, aes(JCvalue1, phyvalue1)) + geom_tile(aes(fill = mean),
                                                                        colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")+
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0),labels=xlabels) +scale_y_discrete(expand = c(0, 0),breaks=seq(0,0.5,1),
                                                                      limits=c(0,1),labels=ylabels)
mylegend = g_legend(legend_plot)



m = matrix(1:3,1,3)

grob1 = textGrob("The strength of phylogenetic effect", gp=gpar(fontsize=16),rot = 90)

grob2 = arrangeGrob(grobs = tas1, layout_matrix = m,respect=TRUE)
grob4 = textGrob("The strength of J-C effect", gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")



grid.arrange(grob1,grob2,mylegend,grob3,grob4,grob3,ncol = 3, widths = c(1,24,1),heights = c(20,1))

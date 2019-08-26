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
evenness.df = NULL
max.logabund = 12
i_n=3
sce = scenario[i_n]
f.name = sce.short[i_n]
i.true.order = 1
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
      evenness = tsallis(Rs,scales=1,norm=TRUE)
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



p <- ggplot(evenness.df, aes(JCvalue1, phyvalue1)) + geom_tile(aes(fill = mean),
       colour = "white") + scale_fill_gradient(low = "white", high = "steelblue",
                                               breaks=c(0,0.5,1),labels=c("Minimum",0.5,"Maximum"),
                                               limits=c(0,1))
base_size <- 9
p+ theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",axis.ticks = element_blank(), 
       axis.text.x = element_text(size = base_size *2, angle = 0, hjust = 0, colour = "grey50"))



filenames <-  paste0(plot_dir,'abundance_dis_',f.name,'.png')
ggsave(filenames,sce.plot,width = 15, height = 10,dpi=300)


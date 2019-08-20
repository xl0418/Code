library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/1e+07/'
plot_dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
scenario = c('LR','MR','HR')
sce.short = c('L','M','H')

abund.df = NULL
max.logabund = 12
i_n=3
sce = scenario[i_n]
f.name = sce.short[i_n]
i.true.order = 1
for(i in c(1,4,2,5,3)){
  for(j in c(1:6)){
  abund=NULL
  if(i==1){
    comb=paste0(i,i)
  }else{
    comb=paste0(i,j)
  }
  
  for(rep in c(1:100)){
    rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
    Rs = read.csv(rname,header = FALSE)
    log.Rs = log10(Rs)
    freq = hist(as.numeric(log.Rs),plot=FALSE,breaks = seq(0,5.5,0.5))
    counts = freq$counts
    abund = rbind(abund, counts)
  }
  mean.sim = apply(abund,MARGIN=2,FUN=mean)
  sd.sim = sqrt(apply(abund,MARGIN=2,FUN=var))
  col.quan = length(mean.sim)
  if(col.quan<12){
    mean.sim <- c(mean.sim,matrix(0,1,12-col.quan))
    sd.sim <- c(sd.sim,matrix(0,1,12-col.quan))
    
  }
  abund.df = rbind(abund.df,cbind(mean.sim,sd.sim,i.true.order,j,c(1:12)))
  }
  i.true.order=i.true.order+1
}

colnames(abund.df) <- c('mean','sd','JCvalue','phyvalue','species')
abund.df <- as.data.frame(abund.df)

abund.df$JCvalue1 <- factor(abund.df$JCvalue, labels = c('psi==0','0.25','0.5','0.75','1'))
abund.df$phyvalue1 <- factor(abund.df$phyvalue, labels = c('sigma[phi]==0','10^2',
                                                           '10^4','10^6','10^8','Inf'))


sce.plot <- ggplot(abund.df) +
  geom_bar( aes(x=species, y=mean), stat="identity", fill="#6A4028", alpha=0.7) +
  geom_errorbar( aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="#F596AA", alpha=0.7, size=1.3)+
  facet_grid(JCvalue1~phyvalue1, 
             labeller = label_parsed)+
  theme_gdocs()+ scale_color_calc()+
  scale_x_continuous(name="Abundance (log)", breaks=seq(1,12,1),labels = c(freq$mids,5.75)) +
  scale_y_continuous(name="Number of species",breaks=seq(0,60,20))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5))

filenames <-  paste0(plot_dir,'abundance_dis_',f.name,'.png')
ggsave(filenames,sce.plot,width = 15, height = 10,dpi=300)


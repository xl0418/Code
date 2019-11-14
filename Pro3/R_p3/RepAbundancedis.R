library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
plot_dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = c('H','M','L')
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

x.breaks = seq(0,18,1)
max.logabund = 12
  
for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  i.true.order = 1
  abund.df = NULL
  
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
    abund=NULL
    comb = paste0(i,j)
    print(paste(letter.comb,i,j))
    for(rep in c(1:100)){
      if(i_n == 8 & i == 2 & j==4 & rep %in% c(55,58)){
        next
      }
      rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'R',comb,'rep',rep,'.csv')
      Rs = read.csv(rname,header = FALSE)
      log.Rs = log2(Rs)
      freq = hist(as.numeric(log.Rs),plot=FALSE,breaks = x.breaks)
      counts = freq$counts
      abund = rbind(abund, counts)
    }
    mean.sim = apply(abund,MARGIN=2,FUN=mean)
    sd.sim = sqrt(apply(abund,MARGIN=2,FUN=var))
    col.quan = length(mean.sim)
    if(col.quan<length(x.breaks)){
      mean.sim <- c(mean.sim,matrix(0,1,length(x.breaks)-col.quan))
      sd.sim <- c(sd.sim,matrix(0,1,length(x.breaks)-col.quan))
      
    }
    abund.df = rbind(abund.df,cbind(mean.sim,sd.sim,i.true.order,j,c(1:length(x.breaks))))
    }
    i.true.order=i.true.order+1
  }
  
  colnames(abund.df) <- c('mean','sd','JCvalue','phyvalue','species')
  abund.df <- as.data.frame(abund.df)
  
  abund.df$JCvalue1 <- factor(abund.df$JCvalue, labels = c('psi==0','psi==0.25','psi==0.5',
                                                           'psi==0.75','psi==1'))
  abund.df$phyvalue1 <- factor(abund.df$phyvalue, labels = c('sigma[phi]==0','sigma[phi]==10^2',
                                                             'sigma[phi]==10^4','sigma[phi]==10^6',
                                                             'sigma[phi]==10^8','sigma[phi]==Inf'))
  
  interleave <- function(x,y){
    lx <- length(x)
    ly <- length(y)
    n <- max(lx,ly)
    as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
  }
  
  d <- data.frame(x=1:10, y=rnorm(10))
  
  my_labs <- interleave(seq(1,length(x.breaks),2), "")
  my_labs = my_labs[1:19]
  
  
  sce.plot <- ggplot(abund.df) +
    geom_bar( aes(x=species, y=mean), stat="identity", fill="red", alpha=0.7) +
    geom_errorbar( aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="blue", alpha=0.7, size=1.3)+
    geom_line(aes(species, mean),size=0.8,color="blue")+
    facet_grid(JCvalue1~phyvalue1, 
               labeller = label_parsed)+
    #theme_gdocs()+ #scale_color_calc()+
    scale_x_continuous(name="Abundance (log2)", breaks=seq(1,length(x.breaks),1),labels = my_labs) +
    scale_y_continuous(name="Number of species",breaks=seq(0,60,20))+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_blank(),strip.text.x = element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"))
  
  filenames <-  paste0(plot_dir,'abundance_dis_',scefolder,'.png')
  ggsave(filenames,sce.plot,width = 15, height = 10,dpi=300)
}

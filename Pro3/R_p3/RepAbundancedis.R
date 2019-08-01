library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/barplot3d.R', echo=TRUE)
moviedir = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/'
dir = 'C:/Liang/Googlebox/Research/Project3/batchsim_results/1e+07/'
scenario = c('LR','MR','HR')
sce.short = c('L','M','H')
jclabel = c(0,0.5,1)
plabel = c(0,1e2,1e4,1e6,1e8,-1)
diversity.upperlimit = 250
diversity.lowerlimit = 80
lowcol = "#F2F0F7"
highcol = '#54278F'
tas1 = list()
abund.df = NULL
max.logabund = 12
probs = c(0,.25,.5,.75,1)
i_n=1
sce = scenario[i_n]
f.name = sce.short[i_n]
i=3
for(j in c(1:6)){
  abund=NULL
  for(rep in c(61:100)){
    comb=paste0(i,j)
    rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
    Rs = read.csv(rname,header = FALSE)
    log.Rs = log10(Rs)
    freq = hist(as.numeric(log.Rs),plot=FALSE)
    counts = freq$counts
    abund = rbind(abund, counts)
  }
  # quantile1 = apply(abund, MARGIN = 2, FUN = quantile, probs = probs)
  mean.sim = apply(abund,MARGIN=2,FUN=mean)
  sd.sim = sqrt(apply(abund,MARGIN=2,FUN=var))
  col.quan = length(mean.sim)
  if(col.quan<12){
    mean.sim <- c(mean.sim,matrix(0,1,12-col.quan))
    sd.sim <- c(sd.sim,matrix(0,1,12-col.quan))
    
  }
  abund.df = rbind(abund.df,cbind(mean.sim,sd.sim,j,c(1:12)))
}

colnames(abund.df) <- c('mean','sd','phyvalue','species')
abund.df <- as.data.frame(abund.df)


ggplot(abund.df) +
  geom_bar( aes(x=species, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  facet_grid(vars(phyvalue))+ theme_calc()+ scale_color_calc()+
  scale_x_continuous(name="Species", breaks=seq(1,12,1)) +
  scale_y_continuous(name="Frequency")



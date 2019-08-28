library(ape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
source('C:/Liang/Code/Pro3/R_p3/noplotltt95.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/Shadedplotltt95.R', echo=TRUE)
dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')
jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('0','1e2','1e4','1e6','1e8','Inf')
i_n=1
sce = scenario[i_n]
f.name = sce.short[i_n]

i.true.order = 1
windows() ## create window to plot your file
par(mfrow = c(2, 2),omi=c(0.5,0.3,0,0))
for(i in c(4,5)){
  for(j in c(3,5)){
    if(i==1){
      comb=paste0(i,i)
    }else{
      comb=paste0(i,j)
    }    
    multitreefile <- paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/multitree',f.name,comb,'.tre')
    trees <- read.tree(multitreefile)
    mltt <- noplotltt95(trees,log=FALSE,plot=FALSE,mode="mean",alpha=0.05)
    plot.ltt95(mltt,xaxis="flipped",bg=rgb(1,0,0,0.25),shaded=TRUE)
  }
  i.true.order=i.true.order+1
  
}
## ... your plotting code here ...
dev.off() 
# 
# primates<-read.nexus(file=
#                        "http://www.phytools.org/TS2018/data/10kTrees_Primates.nex")
# ltt95(primates,log=TRUE)
# primates<-lapply(primates,force.ultrametric)
# class(primates)<-"multiPhylo"
# object<-ltt95(primates,log=TRUE,plot=FALSE,mode="mean",alpha=0.05)
# par(bty="l")
# plot(object,shaded=TRUE,xaxis="flipped",bg='red')

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

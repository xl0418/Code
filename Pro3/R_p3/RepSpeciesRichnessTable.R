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
scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')
jclabel = c(0,0.5,1)
plabel = c(0,1e2,1e4,1e6,1e8,-1)
diversity.upperlimit = 250
diversity.lowerlimit = 80
lowcol = "#F2F0F7"
highcol = '#54278F'
tas1 = list()
r_df = NULL

for(i_n in c(1:3)){
  sce = scenario[i_n]
  f.name = sce.short[i_n]
  for(i in c(1:3)){
    for(j in c(1:6)){
      ns.vec=NULL
      for(rep in c(61:100)){
        comb=paste0(i,j)
        rname = paste0(dir,'spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
        events = read.csv(rname,header = FALSE)
        ns = tail(events[,2],n=1)
        ns.vec = c(ns.vec, ns)
      }
      quantile1 = quantile(ns.vec)
      r_df = rbind(r_df,c(quantile1,i,j,i_n))
    }
  }
}
z = r_df[,c(1:5)]
z[,1] = min(z)

# Figure size
r3dDefaults$windowRect <- c(0, 100, 800, 800) 

# Animate the 3d plot and save as a gif
open3d()
barplot3d(z,alpha=0.6,mode='m5')
if (!rgl.useNULL())
  # play3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20)
  movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20,dir=moviedir,fps=10)


## Save your images to files if you wish
barplot3d(z,alpha=0.6,mode='m5')
rgl.snapshot(filename=paste0(moviedir,"example.png"))


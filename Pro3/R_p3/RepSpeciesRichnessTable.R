library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
# source('C:/Liang/Code/Pro3/R_p3/barplot3d.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/multi3dbar.R', echo=TRUE)
moviedir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
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

r_df = NULL

for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb = paste0(i,j)
      ns.vec=NULL
      for(rep in c(1:100)){

        rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'psi',i,'s_phi',j,'rep',rep,'Ltable.csv')
        L.table = read.csv(rname,header = FALSE)
        ns = nrow(L.table)
        ns.vec = c(ns.vec, ns)
      }
      quantile1 = quantile(ns.vec)
      r_df = rbind(r_df,c(quantile1,i,j,i_n))
    }
  }
}
z = r_df[,c(1:5)]
z[,1] = min(z[,1]) # min(z)

# Figure size
r3dDefaults$windowRect <- c(0, 100, 800, 800) 

# Animate the 3d plot and save as a gif
open3d()
barplot3d(z,alpha=0.6)
legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
         col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))

if (!rgl.useNULL())
  # play3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20)
  movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20,dir=moviedir,fps=10)


## Save your images to files if you wish
barplot3d(z,alpha=0.6)
legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
         col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))



testz = z
# Test 15*18 9 groups
color.plate =brewer.pal(n = 9, name = 'YlOrRd')
  # c('#DB4D6D','#F596AA','#FEDFE1',
  #               '#CB4042','#EB7A77','#F19483',
  #               '#D0104C','#B5495B','#E87A90'
  #               )
# Figure size
r3dDefaults$windowRect <- c(0, 200, 1600, 1600) 

# Animate the 3d plot and save as a gif
open3d()
barplot3d(testz,group.dim=c(15,18),alpha=0.5,barcolors = color.plate,
          y.intercept=c(35,150,250,400),gap = 0.2,scalexy = 50)
# legend3d("topright", legend = c('Low spatiality','Intermediate spatiality','High spatiality'), pch = 12,
#          col = c("#8CD790","#EFDC05","#30A9DE"), cex=1.5, inset=c(0.02))
rgl.viewpoint( theta = 28, phi = 39, fov = 60, zoom = 1, 
               scale = par3d("scale"), interactive = TRUE, 
               type = c("userviewpoint", "modelviewpoint") )
rgl.snapshot(filename=paste0(moviedir,"Species_richness.png"))

if (!rgl.useNULL())
  # play3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20)
  movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20,dir=moviedir,fps=10)

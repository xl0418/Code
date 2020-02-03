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
dir.result = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'

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





for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      species.area = NULL
      
      comb = paste0(i,j)
      for(local.scale in seq(10,333,10)){
        ns.vec=NULL
        print(paste0('i_n = ',i_n,'...','i = ', i, '; j = ',j,'; area = ',local.scale))
        for(rep in c(1:100)){
          rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'M',i,j,'rep',rep,'.csv')
          L.table = read.csv(rname,header = FALSE)
          global.matrix = as.matrix(L.table)
  
            submatrix.vec = c(0:(333 %/% local.scale - 1))*local.scale+1
            for(row.num in submatrix.vec){
              local.grid = global.matrix[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
              local.richness = length(unique(as.vector(as.matrix(local.grid))))
              ns.vec = c(ns.vec, local.richness)
            }
            
          }
        quantile1 = quantile(ns.vec)
        species.area = rbind(species.area,c(quantile1,i,j,i_n,local.scale))
        
      }
      colnames(species.area) = c('0','25','50','75','100','i','j','i_n','area')
      data.save.name = paste0(dir.result,'/speciesarea/',letter.comb,i,j,'.Rda')
      save(species.area,file=data.save.name)
      # load(data.save.name)
    }
  }
}


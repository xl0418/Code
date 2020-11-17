library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
# source('E:/Code/Pro3/R_p3/barplot3d.R', echo=TRUE)
source('E:/Code/Pro3/R_p3/multi3dbar.R', echo=TRUE)
moviedir = 'E:/Googlebox/Research/Project3/replicate_sim_9sces_results/'
dir = 'E:/Googlebox/Research/Project3/replicate_sim_9sces/'
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
plot.letter = c('(A)', '(B)', '(C)', '(D)')
backward.time.vec = c(0, 1e4, 1e5, 5e5, 1e6)
# backward.time = 1e3 # 1e4, 1e5, 1e6
color.sce = c('red','purple','green','yellow')
scale.vec = rev(c(333,111,55))

# for( sce.rich in c(1:length(scale.vec))){
sce.rich = 1
r_df1 = NULL
r_df2 = NULL
r_df3 = NULL
r_df4 = NULL
r_df5 = NULL

local.scale = scale.vec[sce.rich]


for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb = paste0(i,j)
      ns.vec1=NULL
      ns.vec2=NULL
      ns.vec3=NULL
      ns.vec4=NULL
      ns.vec5=NULL
      
      for(rep in c(1:100)){
        rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'M',i,j,'rep',rep,'.csv')
        lname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'psi',i, 's_phi',j,'rep',rep,'Ltable.csv')
        spatial.table = read.csv(rname,header = FALSE)
        l.table = read.csv(lname, header = FALSE)
        speciation.time = l.table[,1]
        
        species.id.protracted1 = which(speciation.time < backward.time.vec[1])
        global.matrix1 = as.matrix(spatial.table)
        global.matrix1[global.matrix1 %in% species.id.protracted1] <- 0
        
        species.id.protracted2 = which(speciation.time < backward.time.vec[2])
        global.matrix2 = as.matrix(spatial.table)
        global.matrix2[global.matrix2 %in% species.id.protracted2] <- 0
        
        species.id.protracted3 = which(speciation.time < backward.time.vec[3])
        global.matrix3 = as.matrix(spatial.table)
        global.matrix3[global.matrix3 %in% species.id.protracted3] <- 0
        
        species.id.protracted4 = which(speciation.time < backward.time.vec[4])
        global.matrix4 = as.matrix(spatial.table)
        global.matrix4[global.matrix4 %in% species.id.protracted4] <- 0
        
        species.id.protracted5 = which(speciation.time < backward.time.vec[5])
        global.matrix5 = as.matrix(spatial.table)
        global.matrix5[global.matrix5 %in% species.id.protracted5] <- 0
        
        submatrix.vec = c(0:(333 %/% local.scale - 1))*local.scale+1
        for(row.num in submatrix.vec){
          local.grid1 = global.matrix1[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
          local.richness1 = length(unique(as.vector(as.matrix(local.grid1))))
          ns.vec1 = c(ns.vec1, local.richness1)
          
          local.grid2 = global.matrix2[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
          local.richness2 = length(unique(as.vector(as.matrix(local.grid2))))
          ns.vec2 = c(ns.vec2, local.richness2)
          
          local.grid3 = global.matrix3[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
          local.richness3 = length(unique(as.vector(as.matrix(local.grid3))))
          ns.vec3 = c(ns.vec3, local.richness3)
          
          local.grid4 = global.matrix4[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
          local.richness4 = length(unique(as.vector(as.matrix(local.grid4))))
          ns.vec4 = c(ns.vec4, local.richness4)
          
          local.grid5 = global.matrix5[row.num:(row.num+local.scale-1),row.num:(row.num+local.scale-1)]
          local.richness5 = length(unique(as.vector(as.matrix(local.grid5))))
          ns.vec5 = c(ns.vec5, local.richness5)
        }
        
      }
      quantile11 = quantile(ns.vec1)
      quantile12 = quantile(ns.vec2)
      quantile13 = quantile(ns.vec3)
      quantile14 = quantile(ns.vec4)
      quantile15 = quantile(ns.vec5)
      
      r_df1 = rbind(r_df1,c(quantile11,i,j,i_n))
      
      quantile12 = quantile(ns.vec2)
      r_df2 = rbind(r_df2,c(quantile12,i,j,i_n))
      
      quantile13 = quantile(ns.vec3)
      r_df3 = rbind(r_df3,c(quantile13,i,j,i_n))
      
      quantile14 = quantile(ns.vec4)
      r_df4 = rbind(r_df4,c(quantile14,i,j,i_n))
      
      quantile15 = quantile(ns.vec5)
      r_df5 = rbind(r_df5,c(quantile15,i,j,i_n))
    }
  }
}

result.matrix.list = list(r_df1, r_df2, r_df3, r_df4, r_df5)
for ( matrix.index in seq(1,5)){
  r_df = result.matrix.list[[matrix.index]]
  z = r_df[,c(1:5)]
  z[,1] = min(z[,1]) # min(z)
  
  testz = z
  # Test 15*18 9 groups
  # color.plate =brewer.pal(n = 9, name = 'YlOrRd')
  color.gradation <- colorRampPalette(c(color.sce[sce.rich], "white"))
  color.plate <- color.gradation(9)
  # c('#DB4D6D','#F596AA','#FEDFE1',
  #               '#CB4042','#EB7A77','#F19483',
  #               '#D0104C','#B5495B','#E87A90'
  #               )
  # Figure size
  r3dDefaults$windowRect <- c(0, 200, 1600, 1600) 
  
  # Animate the 3d plot and save as a gif
  open3d()
  barplot3d(testz,group.dim=c(15,18),alpha=0.7,barcolors = color.plate,
            y.intercept=c(0,20,40,60,80),gap = 0.2,scalexy = 10)
  legend3d("topleft", legend = paste0(plot.letter[sce.rich],' Scale ',local.scale), 
           cex=1.5, inset=c(0.02))
  rgl.viewpoint( theta = 45, phi = 45, fov = 60, zoom = 1, 
                 scale = par3d("scale"), interactive = TRUE, 
                 type = c("userviewpoint", "modelviewpoint") )
  rgl.snapshot(filename=paste0(moviedir,"Species_richness_scale_protracted_", backward.time.vec[matrix.index], "_local_",local.scale,".png"))
  # if (!rgl.useNULL())
  #   # play3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20)
  #   movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), duration = 20,dir=moviedir,fps=10)
}
  
  


  

# 

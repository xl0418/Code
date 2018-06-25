library(DDD)
library(MASS)
library(rgl)
library(stringr)
library(phytools)
library(plotrix)
source(paste0('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/L2phylo_loc.R'))
source(paste0('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/splitAt.R'))

#cycle = T: open tree label that is tree tips will be labeled by cycles; cycle = F: close
#ls = 1: legend, crown age and titles will be drawn

plottree <- function(age = 15, dropextinct = T,file = filename, ls = 1,cycle = T){
      maps_loc_new = list()
      maps_loc_color = list()
      
     load(file = file)
      L = result$L
      track_loc = result$track_loc
      phy = L2phylo_loc(L,loc = 1,dropextinct = dropextinct)
      lin = as.numeric(phy$tip.label)
      phy$tip.label = L[lin,5]
     
      
      if(dropextinct == F) {break_point = round(unname(as.numeric(result$brts[,1])),digits = 6)}
      #branching time as the split positions
      
      
      if(dropextinct == T){brts=unname(sort(branching.times(phy),decreasing = T))
      break_point = round(age-brts,digits = 6)}
      
      lin_sur_position=which(L[,4]==-1)
      lin_sur_index = which(lin %in% lin_sur_position)
      j=0
      #add the present time to the survived lineages
      for(k in lin[lin_sur_index]){
        track_loc[[k]] = c(track_loc[[k]],age)
      }
      
      #split the list consisting of vectors that contain the states of species and duration time. 
      
      for(i in lin){
        j=j+1
        names_new=names(track_loc[[i]] )[-length(names(track_loc[[i]] ))]
        maps_loc_new[[j]] = diff( track_loc[[i]] )
        names(maps_loc_new[[j]]) = names_new
        names(maps_loc_new[[j]])[length(names(maps_loc_new[[j]]))] = c(paste0(names(maps_loc_new[[j]])[length(names(maps_loc_new[[j]]))],".",i))
        
        pos = which(round( track_loc[[i]],digits = 6) %in% break_point)[-1]
        maps_loc_color = c(maps_loc_color,splitAt(maps_loc_new[[j]],pos))
      }
      map_loc_real_color = unique(maps_loc_color)
     
      for(c in 1:length(map_loc_real_color)){
        names(map_loc_real_color[[c]])[length(map_loc_real_color[[c]])] = c(paste(trunc( as.numeric(names(map_loc_real_color[[c]])[length(map_loc_real_color[[c]])]))))
      }
      cols<-setNames(c("blue","red","green"),c("1","2","3"))
      phy$maps = map_loc_real_color
      
      #add legend for the colors
      if(ls == 1)  {
        plotSimmap(phy,cols,offset = 10,ylim=c(-4,Ntip(phy)),fsize = 1,mar = c(4,1,1,1),split.vertical=TRUE)
        cols_legend = cols 
      names(cols_legend) = c("   Endemic","   Endemic","   Widespread")
      add.simmap.legend(colors=cols_legend,prompt=FALSE,x=0,y=-2,vertical=FALSE)
      axis(1)
      #add the title
      title(xlab="Time from the root")
      }
      #add the x-axis
      else plotSimmap(phy,cols,offset = 10,fsize = 1,split.vertical=TRUE)
      
      #handy function: Click where you want to draw the legend
      # add.simmap.legend(colors=cols,prompt=TRUE, vertical=FALSE)
      
      #loc info on the tips
      if(dropextinct == T & cycle == T){
        phy_loc = result$loctable[lin,]
        
        obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
        x.tip<-obj$xx[1:obj$Ntip]
        y.tip<-obj$yy[1:obj$Ntip]
        
        x_cor_vector1 = x.tip+1*strwidth("W")
        for(i in 1:length(x.tip)){
          x_cor = x_cor_vector1[i]
          y_cor = y.tip[i]
          if(phy_loc[,1][i] == 1) draw.circle(x_cor,y_cor,nv=1000,radius=0.08,col="purple")
          else       draw.circle(x_cor,y_cor,nv=1000,radius=0.08)
          
        }
        x_cor_vector2 = x.tip+2*strwidth("W")
        for(i in 1:length(x.tip)){
          x_cor = x_cor_vector2[i]
          y_cor = y.tip[i]
          if(phy_loc[,2][i] == 1) draw.circle(x_cor,y_cor,nv=1000,radius=0.08,col="purple")
          else       draw.circle(x_cor,y_cor,nv=1000,radius=0.08)    }
        
      }
}

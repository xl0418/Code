library("ggplot2")
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))

tas1 = list()
tas2 = list()
tas3 = list()

Ptable = list()
kpri = c(20,30,40)
for(i_n in c(1:3)){
  pp_group_all = NULL
  
  tens_vec = c(1:4)
if(kpri[i_n] == 20) powergroup = c(1:4)
if(kpri[i_n] == 40) powergroup = c(5:8)
if(kpri[i_n] == 30) powergroup = c(9:12)


for(i_m in tens_vec){
  tens = powergroup[i_m]
 
  j = 2
  p = NA
  power = NA
  
  a = 100+10*(tens-1)+1
  b = 100+10*(tens-1)+9
  #for mu = 0
  for(num in c(a:b)){
    if(kpri[i_n] == 20){
      #K'=20
      if(num %in% c(101:109)) num1 = 101
      else if (num %in% c(111:119)) num1 = 111
      else if (num %in% c(121:129)) num1 = 121
      else if (num %in% c(131:139)) num1 = 131
    }
    
    if(kpri[i_n] == 40){
      # K'=40
      if(num %in% c(141:149)) num1 = 141
      else if (num %in% c(151:159)) num1 = 151
      else if (num %in% c(161:169)) num1 = 161
      else if (num %in% c(171:179)) num1 = 171
    }
    if(kpri[i_n] == 30){
      # K'=40
      if(num %in% c(181:189)) num1 = 181
      else if (num %in% c(191:199)) num1 = 191
      else if (num %in% c(201:209)) num1 = 201
      else if (num %in% c(211:219)) num1 = 211
    }
    
    
    if(kpri[i_n] == 20) file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    if(kpri[i_n] == 40)file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    if(kpri[i_n] == 30)file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    
    print(num)
    
    for(i in 1:100){
      # print(i)
      
      file = paste(file1,"/",num,"analysis",i,".Rdata",sep = "")
      
      load(file = file)
      i <- 100*(num-num1)+i
      
      
      p[i] <- pvalue
      power[i] <- poweroftest
      
      
      
    }
    
    
  }
  
  
  group = rep(c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000),each = 100)
  
  p2 = cbind(group, p , power)
  p2 = as.data.frame(p2)
  
  mu_vec = c(0,0.1,0.2,0.4)
  
  pp_group = aggregate(p2[, 2:3], list(p2$group), mean)
  pp_group = cbind(pp_group,mu_vec[i_m])
  colnames(pp_group)=c("group","p","power","mu")
  pp_group_all = rbind(pp_group_all,pp_group)
  
}
Ptable[[i_n]] = pp_group_all
if(i_n == 1){
tas1[[i_n]] <-    ggplot( pp_group_all , aes(factor(mu), factor(group)  ) ) + 
  geom_tile(aes(fill = power), color = "lightblue" ) +  scale_fill_gradient(low = "white",     high = "steelblue",limits=c(0,1))+
  theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
        axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
 ylab("")+xlab("") +ggtitle("K' = 20")
}
if(i_n == 2){
  tas1[[i_n]] <-    ggplot( pp_group_all , aes(factor(mu), factor(group)  ) ) + 
    geom_tile(aes(fill = power), color = "lightblue" ) +  scale_fill_gradient(low = "white",     high = "steelblue",limits=c(0,1))+
    theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
          axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
    xlab("")+ylab("")+ggtitle("K' = 20 vs. 40")
}
if(i_n == 3){
  tas1[[i_n]] <-    ggplot( pp_group_all , aes(factor(mu), factor(group)  ) ) + 
    geom_tile(aes(fill = power), color = "lightblue" ) +  scale_fill_gradient(low = "white",     high = "steelblue",limits=c(0,1))+
    theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
          axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
    xlab("")+ylab("")+ggtitle("K' = 40")
}

}

legend_plot = ggplot( pp_group_all , aes(factor(mu), factor(group)  ) ) + 
  geom_tile(aes(fill = power), color = "lightblue" ) +  scale_fill_gradient(low = "white",     high = "steelblue",limits=c(0,1))+
  theme(axis.ticks = element_blank(),axis.text.x=element_text(size = 12,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y=element_text(size = 12,angle=90,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
        axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  xlab("Extinction")+ylab("dispersal")+ggtitle("K'=20")
mylegend = g_legend(legend_plot)



m = matrix(1:3,1,3)

grob1 = textGrob("Dispersal rate", gp=gpar(fontsize=16),rot = 90)

grob2 = arrangeGrob(grobs = tas1, layout_matrix = m,respect=TRUE)
grob4 = textGrob("Extinction rate", gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")



grid.arrange(grob1,grob2,mylegend,grob3,grob4,grob3,ncol = 3, widths = c(1,24,1),heights = c(20,1))

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
source('C:/Liang/Code/Pro3/R_p3/multipleplot.R', echo=TRUE)
dir = 'C:/Liang/Googlebox/Research/Project3/simdata_1e+07newpara/1e+07/'
scenario = c('LR','MR','HR')
jclabel = c(0,.2,.4,.6,.8,1)
plabel = c(0,.1,.3,1,3,10)
diversity.upperlimit = 1500
diversity.lowerlimit = 80
lowcol = "#F2F0F7"
highcol = '#54278F'
tas1 = list()

for(i_n in c(1:3)){
  sce = scenario[i_n]
  
  r_df = NULL
  for(i in c(1:6)){
    for(j in c(1:6)){
      comb=paste0(i,j)
      rname = paste0(dir,sce,comb,'.Rdata')
      source(rname)
      r_df_g = cbind(ncol(R),plabel[i],jclabel[j])
      r_df = rbind(r_df,r_df_g)
      
    }
  }
  df = as.data.frame(r_df)
  colnames(df) = c('Richness','phi','psi')
  df$psi = as.character(df$psi)
  df$phi = as.character(df$phi)
  df$phi = factor(df$phi,levels = plabel )
  
  if(i_n == 1){
    tas1[[i_n]] <-    ggplot( df , aes(factor(psi), factor(phi)  ) ) + 
      geom_tile(aes(fill = Richness), color = "lightblue" ) +  scale_fill_gradient(low = lowcol,     high = highcol,limits=c(diversity.lowerlimit,diversity.upperlimit))+
      theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(angle = 90,size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
            axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
      ylab("")+xlab("") +ggtitle("Low spatial strength")
  }
  if(i_n == 2){
    tas1[[i_n]] <-    ggplot( df , aes(factor(psi), factor(phi)  ) ) + 
      geom_tile(aes(fill = Richness), color = "lightblue" ) +  scale_fill_gradient(low = lowcol,     high = highcol,limits=c(diversity.lowerlimit,diversity.upperlimit))+
      theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(angle = 90,size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
            axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
      xlab("")+ylab("")+ggtitle("Mediate spatial strength")
  }
  if(i_n == 3){
    tas1[[i_n]] <-    ggplot( df , aes(factor(psi), factor(phi)  ) ) + 
      geom_tile(aes(fill = Richness), color = "lightblue" ) +  scale_fill_gradient(low = lowcol,     high = highcol,limits=c(diversity.lowerlimit,diversity.upperlimit))+
      theme(legend.position="none",axis.ticks = element_blank(),axis.text.x=element_text(angle = 90,size = 12,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
            axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
      xlab("")+ylab("")+ggtitle("High spatial strength")
  }
  
}

legend_plot = ggplot( df , aes(factor(psi), factor(phi)  ) ) + 
  geom_tile(aes(fill = Richness), color = "lightblue" ) +  scale_fill_gradient(low = lowcol,     high = highcol,limits=c(diversity.lowerlimit,diversity.upperlimit))+
  theme(axis.ticks = element_blank(),axis.text.x=element_text(size = 12,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y=element_text(size = 12,angle=90,hjust=0.5,vjust=0.5),panel.background = element_blank(),plot.margin = unit(c(1,0.5,0.5,1), "cm"),
        axis.title.x=element_text(size = 14,hjust=0.5,vjust=1.5),axis.title.y=element_text(size = 14,angle=90,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  xlab("Extinction")+ylab("dispersal")+ggtitle("High spatial strength")
mylegend = g_legend(legend_plot)



m = matrix(1:3,1,3)

# grob1 <- arrangeGrob(K_grob, lambda_grob,mu_grob,ncol = 1)
grob1 = textGrob("Phylogenetic strength", gp=gpar(fontsize=16),rot = 90)

grob2 = arrangeGrob(grobs = tas1, layout_matrix = m,respect=TRUE)
grob4 = textGrob("Jazen-Connell strength", gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")



# grid.arrange(g2,g4,heights = c(1,32),left = textGrob("Number of lineages", rot = 90, vjust = 1.5))
# grid.arrange(grob1,grob2,grob3,grob4,ncol = 2,widths = c(1,16),heights = c(16,1))

grid.arrange(grob1,grob2,mylegend,grob3,grob4,grob3,ncol = 3, widths = c(1,24,2),heights = c(10,1))

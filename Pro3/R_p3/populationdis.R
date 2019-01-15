library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
source(paste0(getwd(),'/g_legend.R'))
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)
source('C:/Liang/Code/Pro3/R_p3/multipleplot.R', echo=TRUE)
dir = 'C:/Liang/Code/Pro3/data1_20181220/'

r_df = NULL
for(i in c(1:6)){
  for(j in c(1:6)){
    comb=paste0(i,j)
    rname = paste0(dir,'HR',comb,'.Rdata')
    source(rname)
    r_df_g = cbind(t(log(R)),i,j)
    r_df = rbind(r_df,r_df_g)
    
  }
}
psi_vec = c(0,.001,.01,.1,.5,1)
phi_vec = psi_vec
popdf = as.data.frame(r_df)
colnames(popdf) = c('Abundance','phi','psi')
popdf$psi = as.character(popdf$psi)
popdf$phi = as.character(popdf$phi)


# Option 1: histogram plot
dislist = list()
for(count in c(1:6)){
  dfpart = popdf[popdf['psi'] == count,]
  
  dislist[[count]] = ggplot(dfpart, aes(x=Abundance)) + 
    geom_histogram(aes(x=Abundance,fill = phi),binwidth=0.7, position="dodge")+
    #geom_density(alpha=0.55)+
    theme(legend.position="none",axis.text.x=element_text(size = 12,hjust=1,vjust=0.5),
          axis.text.y=element_text(size = 12,hjust=0.5,vjust=0.5),panel.background = element_blank(),
          axis.title.y=element_text(size = 14,margin = margin(t = 0, r = 10, b = 0, l = 0)))+
          scale_fill_manual(values = rev(brewer.pal(6,"Purples")))+
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))+
          ylab( substitute(psi == A, list(A = psi_vec[count])))+xlab("")
}

legend_plot = ggplot(dfpart, aes(x=Abundance)) + 
  geom_histogram(aes(x=Abundance,fill = phi),color = 'black',binwidth=0.5, position="dodge")+
  scale_fill_manual(values = rev(brewer.pal(6,"Purples")),labels = phi_vec)+
  guides(fill = guide_legend(title=expression(phi)))
  
mylegend = g_legend(legend_plot)



m = matrix(1:6,6,1)

grob1 = textGrob("Number of species", gp=gpar(fontsize=16),rot = 90)

grob2 = arrangeGrob(grobs = dislist, layout_matrix = m)
grob4 = textGrob("(log) Abundance", gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")
grob5 = textGrob("High spatial strength", gp=gpar(fontsize=16))

grid.arrange(grob3,grob5,grob3,grob1,grob2,mylegend,grob3,grob4,grob3,ncol = 3, widths = c(1,24,2),heights = c(1,20,1))

# layout_matrix = gridmatrix



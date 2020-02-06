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

x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16



for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  p = list()
  count1 = 1
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb = paste0(i,j)
      data.load.name = paste0(dir.result,'/speciesareamean/',letter.comb,i,j,'mean.Rda')
      load(data.load.name)
      species.area.df <- as.data.frame(species.area)

  
      area = species.area.df$area
      df_lineage = rbind(c(rep(1,5),i,j,i_n,1),species.area.df)
      
      df_lineage$`0` = log(df_lineage$`0`)
      df_lineage$`25` = log(df_lineage$`25`)
      df_lineage$`50` = log(df_lineage$`50`)
      df_lineage$`75` = log(df_lineage$`75`)
      df_lineage$`100` = log(df_lineage$`100`)
      df_lineage$area = log(df_lineage$area)
      df_min_max = data.frame(id = "min_max", value = 1, x = c(df_lineage$area,rev(df_lineage$area)), y = c(df_lineage$'0',rev(df_lineage$'100')))
      df_0025 = data.frame(id = "0025", value = 2, x = c(df_lineage$area,rev(df_lineage$area)), y = c(df_lineage$'25',rev(df_lineage$'75')))
      df_mean = data.frame(id = 'mean',y = df_lineage$`50`,x = df_lineage$area)
      df_lineage_all = rbind(df_min_max,df_0025)
      
      
      p[[count1]] <- ggplot(df_mean, aes(x = x, y = y)) +
        theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",
              panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.text.x=element_blank(),axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
              ,plot.margin=unit(c(0,0,0,0),"cm"))+
        geom_polygon(data = df_min_max, aes(  group = id),fill = "gray70", alpha = 0.8)+
        geom_polygon(data = df_0025, aes( group = id),fill = "gray27", alpha = 0.8)+
        geom_line()+
        ylab("")+xlab("")+
        coord_cartesian(xlim=c(log(1),log(330)),ylim=c(0,log(1000))) + 
        scale_y_continuous(breaks = log(c(1,10,100)),labels = c('1','10','100'))+
        scale_x_continuous(breaks = log(c(1,10,100)),labels = c('1','10','100'))
      if(count1 %in% c(1:6)){
        p[[count1]] <- p[[count1]]+theme(axis.text.y=element_text(angle=90,size = y_label_fontsize))
      }
      if(count1 %in% seq(6,30,6)){
        p[[count1]] <- p[[count1]]+theme(axis.text.x=element_text(size = x_label_fontsize))
      }
      count1 = count1+1
    }
  }
  m = matrix(1:30,ncol = 5)
  
  
  psi_labels = c('0','0.25','0.5','0.75','1')
  
  
  phi1 <- textGrob(expression(phi == 1), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  phi2 <- textGrob(expression(phi == 10^-2), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  phi3 <- textGrob(expression(phi == 10^-4), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  phi4 <- textGrob(expression(phi == 10^-6), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  phi5 <- textGrob(expression(phi == 10^-8), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  phi6 <- textGrob(expression(phi == 0), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  
  psi1 <- textGrob(substitute(paste(psi," = ",nn),list(nn = psi_labels[1])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  psi2 <- textGrob(substitute(paste(psi," = ",nn),list(nn = psi_labels[2])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  psi3 <- textGrob(substitute(paste(psi," = ",nn),list(nn = psi_labels[3])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  psi4 <- textGrob(substitute(paste(psi," = ",nn),list(nn = psi_labels[4])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  psi5 <- textGrob(substitute(paste(psi," = ",nn),list(nn = psi_labels[5])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  
  
  row_titles <- arrangeGrob(phi1,phi2,phi3,phi4,phi5,phi6,ncol = 1)
  column_titles <- arrangeGrob(psi1,psi2,psi3,psi4,psi5,ncol = 5)
  
  
  label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)
  
  g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
  g_ltt5 = textGrob("area", gp=gpar(fontsize=x_title_fontsize),rot = 0)
  g_ltt1 = textGrob("")
  g_a = textGrob("(a)")
  g_b = textGrob("(b)")
  
  ltt.sce <- grid.arrange(g_ltt1,column_titles,g_ltt1,label,g_ltt4,row_titles,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,16,3),heights = c(1,36,1))
  
  dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
  savefilename <- paste0(dir_save,scefolder,'_species_area_log_new.pdf')
  ggsave(savefilename,ltt.sce,width = 15,height = 10)
}


  
      
      
      
library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
source(paste0(getwd(),'/g_legend.R'))
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = c('H','M','L')
scenario = NULL
sce.short.comb.vec = NULL
method = 'AWMIPD'
for(i.letter in sce.short){
  for(j.letter in sce.short){
    sce.folder = paste0('sce',i.letter,j.letter)
    scenario = c(scenario,sce.folder)
    sce.short.comb = paste0(i.letter,j.letter)
    sce.short.comb.vec = c(sce.short.comb.vec,sce.short.comb)
  }
}

scale.vec = c(333,222,111,55)
# for( sce.rich in c(1:length(scale.vec))){
local.scale = 333
max.D.interval = c()


x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16
# 
# for(i_n in c(1:9)){
#   print(paste0("Processing Scenario ",i_n,"..."))
#   plots_df = list()
#   plot.num = 1
#   scefolder = scenario[i_n]
#   letter.comb = sce.short.comb.vec[i_n]
#   sce = scenario[i_n]
#   
#   for(i in c(1,4,2,5,3)){
#     for(j in c(1:6)){
#       comb = paste0(i,j)
#       rep = 1
#       Dname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'D',i,j,'rep',rep,'.csv')
#       D.table = read.csv(Dname,header = FALSE)
#       D.matrix = as.matrix(D.table)
#       total.dvalues = rowSums(D.matrix)
#       D.interval = (max(total.dvalues)-min(total.dvalues))
#       max.D.interval = c(max.D.interval,D.interval)
#     }
#   }
# }
# max.D = max(max.D.interval)
# 


for(i_n in c(1:9)){
  print(paste0("Processing Scenario ",i_n,"..."))
  plots_df = list()
  plot.num = 1
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb = paste0(i,j)
      rep = 1
      rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'M',i,j,'rep',rep,'.csv')
      L.table = read.csv(rname,header = FALSE)
      Dname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'D',i,j,'rep',rep,'.csv')
      D.table = read.csv(Dname,header = FALSE)
      Rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'R',i,j,'rep',rep,'.csv')
      R.table = read.csv(Rname,header = FALSE)
      global.matrix = as.matrix(L.table)
      D.matrix = as.matrix(D.table)
      if(method == 'abs'){
        total.dvalues = rowSums(D.matrix)
        D.normalized = (total.dvalues-min(total.dvalues))/(max(total.dvalues)-min(total.dvalues))
      }else if (method == 'MIPD'){
        IPD.matrix = 1/D.matrix
        diag(IPD.matrix) = 0
        total.dvalues = rowSums(IPD.matrix)/(nrow(D.matrix)-1)
        D.normalized = (total.dvalues-min(total.dvalues))/(max(total.dvalues)-min(total.dvalues))
      }else if (method == 'AWMIPD'){
        AD.matrix = sweep(D.matrix, MARGIN=2, as.matrix(R.table/local.scale^2), `*`)
        IPD.matrix = 1/AD.matrix
        diag(IPD.matrix) = 0
        total.dvalues = rowSums(IPD.matrix) * as.matrix(R.table/local.scale^2)
        
        D.normalized = (total.dvalues-min(total.dvalues))/(max(total.dvalues)-min(total.dvalues))
      }else{
        print('Provide method...')
        break
      }
      local.x = c(1:local.scale)
      local.y = c(1:local.scale)
      distribution.data = expand.grid(X=local.x,Y=local.y)
      distribution.data$Z = global.matrix[cbind(distribution.data$X,distribution.data$Y)]
      distribution.data$D = D.normalized[distribution.data$Z+1]
      
      plots_df[[plot.num]] <- ggplot(distribution.data, aes(X, Y, fill= D)) + geom_tile()+
        theme(legend.position = '',axis.text = element_blank(),axis.ticks = element_blank(),
              panel.background = element_blank())+
        xlab("")+ylab("") + scale_fill_gradient2(low="#005CAF",mid = 'white',
                                                 high="#D0104C",midpoint=0.5)
      plot.num = plot.num +1
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
  
  g_ltt4 = arrangeGrob(grobs = plots_df, layout_matrix = m)
  g_ltt5 = textGrob("Species distribution", gp=gpar(fontsize=x_title_fontsize),rot = 0)
  g_ltt1 = textGrob("")

  
  ltt.sce <- grid.arrange(g_ltt1,column_titles,g_ltt1,g_ltt1,g_ltt4,
                          row_titles,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,16,3),heights = c(1,24,1))
  
  dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
  savefilename <- paste0(dir_save,scefolder,'_species_dis_AWMIPD.eps')
  ggsave(savefilename,ltt.sce,width = 15,height = 10)
  
}


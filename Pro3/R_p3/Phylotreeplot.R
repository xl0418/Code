library(ggplot2)
library(gridExtra)
library(grid)
library(DDD)
library(ggtree)
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

age = 10000000
x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

for(i_n in c(1:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  
  p = list()
  count1 = 1
  
  j = 2
  combinations = NULL
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb.temp <- paste0(i,j)
      combinations <- c(combinations,comb.temp)
    }
  }
  
  for(num in combinations){
    L3 = NULL
    fre = NULL
    L = list()
    data = NA
    print(num)
    multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,num,'/','multitree',letter.comb,num,'.tre')
    
    trees <- read.tree(multitreefile)
    single.tree <- trees[[15]]
    
    p[[count1]] <- ggtree(single.tree,layout = "circular")+xlim(0,age)
    count1 = count1+1
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
  
  
  label = textGrob("",gp=gpar(fontsize=y_title_fontsize), rot = 90)
  
  g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
  g_ltt5 = textGrob("", gp=gpar(fontsize=x_title_fontsize),rot = 0)
  g_ltt1 = textGrob("")
  g_a = textGrob("(a)")
  g_b = textGrob("(b)")
  
  ltt.sce <- grid.arrange(g_ltt1,column_titles,g_ltt1,label,g_ltt4,row_titles,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,16,3),heights = c(1,36,1))
  
  dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
  savefilename <- paste0(dir_save,scefolder,'_trees_cir.pdf')
  ggsave(savefilename,ltt.sce,width = 15,height = 10)
}
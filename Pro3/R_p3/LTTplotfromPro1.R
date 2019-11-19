library(ggplot2)
library(gridExtra)
library(grid)
library(DDD)
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

for(i_n in c(1:7)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]

  p = list()
  count1 = 1

  j = 2
  combinations = NULL
  for(i in c(1:5)){
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
    for (i in 1:length(trees)) {
      tes = trees[[i]]
      brts= -unname(sort(branching.times(tes),decreasing = T))
      data0 = cbind(i,brts,c(2:(length(brts)+1)))
      if(max(brts)<0) data0 = rbind(data0,c(i,0,data0[nrow(data0),3]))
      data = rbind(data, data0)
      
    }
    data <- data[-1,]
    time = data[order(data[,2]),2]
    timeu = unique(time)
    data_lineage = timeu
    for(i in 1:length(trees)){
      tes = trees[[i]]
      
      brts= -unname(sort(branching.times(tes),decreasing = T))
      M1 = match(brts,timeu)
      M1[1] = 1
      M11 = diff(M1)
      M13 = length(timeu)-max(M1)+1
      M12 = c(M11,M13)
      N1 = rep(2:(length(brts)+1),M12)
      data_lineage = cbind(data_lineage,N1)
    }
    x = data_lineage[,1]
    z = data_lineage[,2:length(trees)+1]
    
    data_average_z <- apply(z, 1, median)
    data_q0.025_z <- apply(z, 1 , quantile, 0.025)
    data_q0.25_z <- apply(z, 1, quantile, 0.25)
    data_q0.75_z <- apply(z, 1, quantile, 0.75)
    data_q0.975_z <- apply(z, 1, quantile, 0.975)
    data_lower_z <- apply(z,1,min)
    data_upper_z <- apply(z,1,max)
    
    lineage_stat = cbind(x,log(data_average_z),log(data_q0.025_z),log(data_q0.25_z),log(data_q0.75_z),log(data_q0.975_z),log(data_lower_z),log(data_upper_z))
    colnames(lineage_stat) = c("time", "median","0.025","0.25","0.75","0.975","min","max")
    
    time = min(lineage_stat[,1])
    df_lineage = data.frame(lineage_stat)
    df_min_max = data.frame(id = "min_max", value = 1, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$min,rev(df_lineage$max)))
    df_0025 = data.frame(id = "0025", value = 2, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$X0.025,rev(df_lineage$X0.975)))
    df_025 = data.frame(id = "025", value = 3, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$X0.25,rev(df_lineage$X0.75)))
    df_lineage_all = rbind(df_min_max,df_025,df_0025)

      
    p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
      theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",
            panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
            axis.text.x=element_blank(),axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
            ,plot.margin=unit(c(0,0,0,0),"cm"))+
      geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
      geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
      geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab("")+
      coord_cartesian(xlim=c(-age,0),ylim=c(log(2),log(600))) + scale_y_continuous(breaks = c(log(2),log(10),log(50),log(400)),labels = c(2,10,50,400))+
       scale_x_continuous(breaks = -rev(seq(0,1e7,1e7/5)),labels = c('10','8','6','4','2','0'))
    if(count1 %in% c(1:6)){
      p[[count1]] <- p[[count1]]+theme(axis.text.y=element_text(angle=90,size = y_label_fontsize))
    }
    if(count1 %in% seq(6,30,6)){
      p[[count1]] <- p[[count1]]+theme(axis.text.x=element_text(size = x_label_fontsize))
    }
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
  
  
  label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)
  
  g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
  g_ltt5 = textGrob("Time", gp=gpar(fontsize=x_title_fontsize),rot = 0)
  g_ltt1 = textGrob("")
  g_a = textGrob("(a)")
  g_b = textGrob("(b)")
  
  ltt.sce <- grid.arrange(g_ltt1,column_titles,g_ltt1,label,g_ltt4,row_titles,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,16,3),heights = c(1,36,1))
  
  dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'
  savefilename <- paste0(dir_save,scefolder,'_ltt.pdf')
  ggsave(savefilename,ltt.sce,width = 15,height = 10)
}
library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
library(ggtree)
library(ape)
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

interleave <- function(x,y){
  lx <- length(x)
  ly <- length(y)
  n <- max(lx,ly)
  as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
}

x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

dispersal_title <- rep(c('high dispersal', 'intermediate dispersal', 'low dispersal'), each = 3)
interaction_title <- rep(c('high interaction distance', 'intermediate interaction distance', 'low interaction distance'), 3)
psi_title <- c('0', '0.5', '1', '0.25', '0.75')
phi_title <- c('0', '-2', '-4', '-6', '-8', 0)


count1 = 1
p = list()
# AWMIPD mode
pdmode = 'exp'
a = 1
rep.sample = 18

plot.combination = rbind(c(1, 3, 3),
                         c(1, 4, 4),
                         c(8, 3, 2))

col_labels = NULL

for(plot.comb in c(1:nrow(plot.combination))){
  i_n = plot.combination[plot.comb,1]
  i = plot.combination[plot.comb,2]
  j = plot.combination[plot.comb,3]
  comb = paste0(i,j)
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  print(paste0('i_n = ',i_n,'; comb = ', comb))
  
  multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')
  
  trees <- read.tree(multitreefile)
  single.tree <- trees[[rep.sample]]
  
  p[[count1]] <- ggtree(single.tree,layout = "circular") #+xlim(0,age)
  
  count1 = count1+1
  
  # plot SAR
  data.load.name = paste0(dir.result,'/speciesareamean/',letter.comb,comb,'mean.Rda')
  load(data.load.name)
  species.area.df <- as.data.frame(species.area)
  area = species.area.df$area
  df_lineage_sar = rbind(c(rep(1,5),i,j,i_n,1),species.area.df)
  
  df_lineage_sar$`0` = log(df_lineage_sar$`0`)
  df_lineage_sar$`25` = log(df_lineage_sar$`25`)
  df_lineage_sar$`50` = log(df_lineage_sar$`50`)
  df_lineage_sar$`75` = log(df_lineage_sar$`75`)
  df_lineage_sar$`100` = log(df_lineage_sar$`100`)
  df_lineage_sar$area = log(df_lineage_sar$area)
  df_min_max_sar = data.frame(id = "min_max", value = 1, x = c(df_lineage_sar$area,rev(df_lineage_sar$area)), y = c(df_lineage_sar$'0',rev(df_lineage_sar$'100')))
  df_0025_sar = data.frame(id = "0025", value = 2, x = c(df_lineage_sar$area,rev(df_lineage_sar$area)), y = c(df_lineage_sar$'25',rev(df_lineage_sar$'75')))
  df_mean_sar = data.frame(id = 'mean',y = df_lineage_sar$`50`,x = df_lineage_sar$area)
  df_lineage_all_sar = rbind(df_min_max_sar,df_0025_sar)
  
  
  p[[count1]] <- ggplot(df_mean_sar, aes(x = x, y = y)) +
    theme(legend.position="none",
          panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
          ,plot.margin=unit(c(0,0,0,0),"cm"))+
    geom_polygon(data = df_min_max_sar, aes(  group = id),fill = "gray70", alpha = 0.8)+
    geom_polygon(data = df_0025_sar, aes( group = id),fill = "gray27", alpha = 0.8)+
    geom_line()+
    coord_cartesian(xlim=c(log(1),log(330)),ylim=c(0,log(1000))) + 
    scale_y_continuous(name = "No. of species",breaks = log(c(1,10,100)),labels = c('1','10','100'))+
    scale_x_continuous(name = "Area",breaks = log(c(1,10,100)),labels = c('1','10','100'))+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
          axis.text.x=element_text(size = x_label_fontsize))
  
  count1 = count1+1
  
  # SAD
  x.breaks = seq(0,17,1)
  
  abund=NULL
  for(rep.no in c(1:100)){
    rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'R',comb,'rep',rep.no,'.csv')
    Rs = read.csv(rname,header = FALSE)
    log.Rs = log2(Rs)
    freq = hist(as.numeric(log.Rs),plot=FALSE,breaks = x.breaks)
    counts = freq$counts
    abund = rbind(abund, counts)
  }
  mean.sim = apply(abund,MARGIN=2,FUN=mean)
  sd.sim = sqrt(apply(abund,MARGIN=2,FUN=var))
  col.quan = length(mean.sim)
  if(col.quan<length(x.breaks)){
    mean.sim <- c(mean.sim,matrix(0,1,length(x.breaks)-col.quan))
    sd.sim <- c(sd.sim,matrix(0,1,length(x.breaks)-col.quan))
    
  }
  abund.df = cbind(mean.sim,sd.sim,c(1:length(x.breaks)))
  colnames(abund.df) <- c('mean','sd','species')
  abund.df <- as.data.frame(abund.df)
  
  my_labs <- interleave(seq(1,length(x.breaks),2), "")
  my_labs = my_labs[1:18]
  
  
  p[[count1]]  <- ggplot(abund.df) +
    geom_bar( aes(x=species, y=mean),width = 0.6, stat="identity", fill="red", alpha=0.7) +
    geom_errorbar( aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="blue", alpha=0.7, size=1.3)+
    geom_line(aes(species, mean),size=0.8,color="blue")+
    #theme_gdocs()+ #scale_color_calc()+
    scale_x_continuous(name="Abundance (log2)", breaks=seq(1,length(x.breaks),1),labels = my_labs) +
    scale_y_continuous(name="Frequency",breaks=seq(0,60,20))+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_blank(),strip.text.x = element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"))
  
  count1 = count1+1
  
  # AWMPID
  rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'M',comb,'rep',rep.sample,'.csv')
  L.table = read.csv(rname,header = FALSE)
  Dname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'D',comb,'rep',rep.sample,'.csv')
  D.table = read.csv(Dname,header = FALSE)
  Rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'R',comb,'rep',rep.sample,'.csv')
  R.table = read.csv(Rname,header = FALSE)
  global.matrix = as.matrix(L.table)
  D.matrix = as.matrix(D.table) * 2 / 10^7
  
  # phylogenetic distance
  if (pdmode == 'inv') {
    ID = 1 / D.matrix
    diag(ID) = 1
    AD.matrix = sweep(ID, MARGIN = 2, as.matrix(as.numeric(R.table)), `*`)
  } else if (pdmode == 'exp') {
    ID = exp(- a * D.matrix)
    AD.matrix = sweep(ID, MARGIN = 2, as.matrix(as.numeric(R.table)), `*`)
  }
  
  total.dvalues = rowSums(AD.matrix) * as.matrix( as.numeric(R.table))
  
  D.normalized =   total.dvalues/sum(total.dvalues)
  D.normalized = (D.normalized-min(D.normalized))/(max(D.normalized)-min(D.normalized))
  local.x = c(1:333)
  local.y = c(1:333)
  distribution.data = expand.grid(X=local.x,Y=local.y)
  distribution.data$Z = global.matrix[cbind(distribution.data$X,distribution.data$Y)]
  distribution.data$D = D.normalized[distribution.data$Z+1]
  
  p[[count1]] <- ggplot(distribution.data, aes(X, Y, fill= D)) + geom_tile()+
    theme(legend.position = '',axis.text = element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank())+
    xlab("")+ylab("") + scale_fill_gradient2(low="#005CAF",mid = 'green',
                                             high="#D0104C",midpoint=0.5)
  
  count1 = count1+1
  
  # LTT plot
  age = 10000000
  data = NULL
  for (i in 1:length(trees)) {
    tes = trees[[i]]
    brts= -unname(sort(branching.times(tes),decreasing = T))
    data0 = cbind(i,brts,c(2:(length(brts)+1)))
    if(max(brts)<0) data0 = rbind(data0,c(i,0,data0[nrow(data0),3]))
    data = rbind(data, data0)
    
  }
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
    theme(legend.position="none",
          panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")
          ,plot.margin=unit(c(0,0,0,0),"cm"))+
    geom_polygon(data = df_min_max, aes(  group = id),fill = "gray80", alpha = 0.8)+
    geom_polygon(data = df_0025, aes( group = id),fill = "gray50", alpha = 0.8)+
    geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+
    coord_cartesian(xlim=c(-age,0),ylim=c(log(2),log(600))) + scale_y_continuous(name="No. of species",breaks = c(log(2),log(10),log(50),log(400)),labels = c(2,10,50,400))+
    scale_x_continuous(name="Time",breaks = -rev(seq(0,1e7,1e7/5)),labels = c('10','8','6','4','2','0'))+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),axis.text.x=element_text(size = x_label_fontsize))
  
  
  
  count1 = count1+1
  col_labels <- c(col_labels,
                  paste0('SPJC ', dispersal_title[i_n], ' &\n', interaction_title[i_n]))
}

m = matrix(1:15,ncol = 3)



row_labels = c('Tree', 'SAR','SAD','SPD','LTT')

phi1 <- textGrob(row_labels[1], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi2 <- textGrob(row_labels[2], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi3 <- textGrob(row_labels[3], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi4 <- textGrob(row_labels[4], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi5 <- textGrob(row_labels[5], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

psi11 <- textGrob(col_labels[1],
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi12 <- textGrob(col_labels[2], 
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi13 <- textGrob(col_labels[3],
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))


psi21 <- textGrob(bquote(psi == .(psi_title[plot.combination[1, 2]]) ~ phi == 10^.(phi_title[plot.combination[1, 3]])),
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi22 <- textGrob(bquote(psi == .(psi_title[plot.combination[2, 2]]) ~ phi == 10^.(phi_title[plot.combination[2, 3]])), 
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi23 <- textGrob(bquote(psi == .(psi_title[plot.combination[3, 2]]) ~ phi == 10^.(phi_title[plot.combination[3, 3]])),
                 gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

row_titles <- arrangeGrob(phi1,phi2,phi3,phi4,phi5,ncol = 1)
column_titles1 <- arrangeGrob(psi11,psi12,psi13,ncol = 3)
column_titles2 <- arrangeGrob(psi21,psi22,psi23,ncol = 3)


# label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)

g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
g_ltt5 = textGrob("area", gp=gpar(fontsize=x_title_fontsize),rot = 0)
g_ltt1 = textGrob("")
g_a = textGrob("(a)")
g_b = textGrob("(b)")

ltt.sce <- grid.arrange(column_titles1,g_ltt1, column_titles2,g_ltt1, g_ltt4,row_titles,ncol = 2,widths = c(20,1),heights = c(2,1,25))

dir_save <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/newmix/'
savefilename <- paste0(dir_save,'fullmixing_example.pdf')
ggsave(savefilename,ltt.sce,width = 15,height = 10)








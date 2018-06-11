library(ggplot2)
library(gridExtra)
# library(nLTT)
library(gridExtra)
# library(DDD)
library(grid)
# source('~/Googlebox/R/MMgit/multiplot.R')

p = list()

Scenario = 2      
locationtype = 3  #1: location 1; 2: location 2; 3 Both location

x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

age = 15
count1 = 1
count2 = 0
lambda = 0.8
k = 20
lambda_allo = 0.2




phi_group = c(1:5)
psi_group = c(1:5)
i= 1
j = 1
mu_vec = rep(psi,each = 5)
for(phiindex in phi_group){
  for(psiindex in psi_group){
  Ntable = NULL
  time = NULL
  L = list()
  Ntable1 = NULL
  Ntable2 = NULL
  x_no2 = NULL
  time_list = NULL
  
  for (seed in 1:100) {
    
    file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",phiindex,"phi",psiindex,"psi","sim",seed,".Rdata",sep = "")
    load(file = file)
    if(max(result$t)<age){
      time = c(result$t,age)
    }
    else time = c(result$t[-length(result$t)],age)
    time_list = c(time_list,time)
    
    
  }
  time = time_list[order(time_list)]
  timeu = unique(time)
  number_loc1 = timeu
  number_loc2 = timeu
  for(i in 1:100){
    
    
    if(Scenario == 2) file = paste("/Volumes/Liang/Research/data/cluster_Kpri_40/",j,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    if(Scenario == 1) file = paste("/Volumes/Liang/Research/data/cluster_Kpri_20/",j,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    if(Scenario == 3) file = paste("/Volumes/Liang/Research/data/cluster_Kpri_20V40/",j,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    
    load(file = file)
    Ntable = rbind(result$Ntable,result$Ntable[nrow(result$Ntable),])
   
    if(locationtype == 1)     Ntable_loc1 = Ntable[,1]+Ntable[,3]
    if(locationtype == 2)     Ntable_loc1 = Ntable[,2]+Ntable[,3]
    if(locationtype == 3)        Ntable_loc1 = Ntable[,1]+Ntable[,3]+Ntable[,2]
    
    if(max(result$t)<age){
      time_each = c(result$t,age)
    }
    else time_each = c(result$t[-length(result$t)],age)
    M1 = match(time_each,timeu)
    M1[1] = 1
    M11 = diff(M1)
    M13 = length(timeu)-max(M1)+1
    M12 = c(M11,M13)
    N1 = rep(Ntable_loc1,M12)
    
    number_loc1 = cbind(number_loc1,N1)
    
    
    x_no2 = c(x_no2, result$Ntable[nrow(result$Ntable),2])
  }
  
  x = median(x_no2)
  
  x_t1 = number_loc1[,1]-age
  z_l1 = number_loc1[,2:101]
  
  data_average_z_l1 <- apply(z_l1, 1, median)
  data_q0.025_z_l1 <- apply(z_l1, 1 , quantile, 0.025)
  data_q0.25_z_l1 <- apply(z_l1, 1, quantile, 0.25)
  data_q0.75_z_l1 <- apply(z_l1, 1, quantile, 0.75)
  data_q0.975_z_l1 <- apply(z_l1, 1, quantile, 0.975)
  data_lower_z_l1 <- apply(z_l1,1,min)
  data_upper_z_l1 <- apply(z_l1,1,max)
 
  # log number of lineages
  species_stat1 = cbind(x_t1,log(data_average_z_l1),log(data_q0.025_z_l1),log(data_q0.25_z_l1),log(data_q0.75_z_l1),log(data_q0.975_z_l1),log(data_lower_z_l1),log(data_upper_z_l1))
  
  colnames(species_stat1) = c("time", "median","0.025","0.25","0.75","0.975","min","max")

  time = min(species_stat1[,1])
  df_species1= data.frame(species_stat1)
  df_min_max_1 = data.frame(id = "min_max", value = 1, x = c(df_species1$time,rev(df_species1$time)), y = c(df_species1$min,rev(df_species1$max)))
  df_0025_1 = data.frame(id = "0025", value = 2, x = c(df_species1$time,rev(df_species1$time)), y = c(df_species1$X0.025,rev(df_species1$X0.975)))
  df_025_1 = data.frame(id = "025", value = 3, x = c(df_species1$time,rev(df_species1$time)), y = c(df_species1$X0.25,rev(df_species1$X0.75)))
  df_lineage_all_1 = rbind(df_min_max_1,df_025_1,df_0025_1)
  mu_vec = rep(c(0,0.1,0.2,0.4),each = 8)
  m0_vec = rep(c(0,0.05,0.1,0.15,0.3,0.5,1,5),4)
  
  if(Scenario == 1 && locationtype == 1)  line_value = 20  #
  if(Scenario == 1 && locationtype == 2)  line_value = 20  #
  if(Scenario == 1 && locationtype == 3)  line_value = 40  #
  if(Scenario == 2 && locationtype == 1)  line_value = 40  #
  if(Scenario == 2 && locationtype == 2)  line_value = 40  #
  if(Scenario == 2 && locationtype == 3)  line_value = 80  #
  if(Scenario == 3 && locationtype == 1)  line_value = 20  #
  if(Scenario == 3 && locationtype == 2)  line_value = 40  #
  if(Scenario == 3 && locationtype == 3)  line_value = 60  #
  
  # allo_term = 
  k1= (lambda-mu_vec[count1])*k/lambda
  K_real =( lambda-mu_vec[count1]-m0_vec[count1]*x/k1+sqrt((m0_vec[count1]*x/k1 -lambda+mu_vec[count1])^2 + 4*m0_vec[count1]*lambda*x/k1))/(2*lambda/k1)
  # print(K_real)
  
  if(count1 %in% c(18,27,36)){
    
    p[[count1]] <- ggplot(df_min_max_1, aes(x = x, y = y)) +
      # scale_y_continuous(breaks = 20,labels = 20)+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y = element_blank(),legend.position="none",
            panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
            axis.text.x=element_text(size = x_label_fontsize), axis.line.x = element_line(color="black"))+
      geom_polygon(data = df_min_max_1, aes(  group = id),fill = "light gray", alpha = 0.8)+
      geom_polygon(data = df_0025_1, aes( group = id),fill = "dark gray", alpha = 0.8)+
      geom_polygon(data = df_025_1, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab(substitute(paste(mu,"=",nn),list(nn = mu_vec[count1])))+
      geom_line(data = df_species1, aes(time, median), col = "black")+#+ggtitle(substitute(paste(mu,"=",nn),list(nn = mu_vec[count1])))+
      geom_hline(yintercept = log(line_value), linetype = "dashed")+
      # geom_hline(yintercept = log(K_real), linetype = "dotted")+
      coord_cartesian(ylim = c(0, 5))+
      theme(plot.margin=unit(c(0,0,0,0),"cm"))#+scale_y_continuous(breaks = sort(c(seq(0, 5, length.out=5), 0.8)))
  }
  
  else if(count1 %in% c(1:8)){
    
    p[[count1]] <- ggplot(df_min_max_1, aes(x = x, y = y)) +
      # scale_y_continuous(breaks = 20,labels = 20)+
      theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background = element_blank(),
            axis.text.y=element_text(size = y_label_fontsize), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line.y = element_line(color="black"))+
      geom_polygon(data = df_min_max_1, aes(  group = id),fill = "light gray", alpha = 0.8)+
      geom_polygon(data = df_0025_1, aes( group = id),fill = "dark gray", alpha = 0.8)+
      geom_polygon(data = df_025_1, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab(substitute(paste(mu,"=",nn),list(nn = mu_vec[count1])))+
      geom_line(data = df_species1, aes(time, median), col = "black")+ 
      # geom_hline(yintercept = log(K_real), linetype = "dotted")+#+ggtitle(substitute(paste(mu,"=",nn),list(nn = mu_vec[count1])))
      geom_hline(yintercept = log(line_value), linetype = "dashed")+
      coord_cartesian(ylim = c(0, 5),xlim=c(-age,0))+scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))+
      theme(plot.margin=unit(c(0,0,0,0),"cm"))
  }
  else if(count1 == 9)
  {
    p[[count1]] <- ggplot(df_min_max_1, aes(x = x, y = y)) +
      # scale_y_continuous(breaks = 20,labels = 20)+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.text.y=element_text(size = y_label_fontsize),axis.text.x=element_text(size = x_label_fontsize), axis.line.y = element_line(color="black"))+
      geom_polygon(data = df_min_max_1, aes(  group = id),fill = "light gray", alpha = 0.8)+
      geom_polygon(data = df_0025_1, aes( group = id),fill = "dark gray", alpha = 0.8)+
      geom_polygon(data = df_025_1, aes( group = id), fill = "gray27", alpha = 0.8)+
      geom_line(data = df_species1, aes(time, median), col = "black")+ 
      # geom_hline(yintercept = log(K_real), linetype = "dotted")+
      geom_hline(yintercept = log(line_value), linetype = "dashed")+coord_cartesian(ylim = c(0, 5))+scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))+
      theme(plot.margin=unit(c(0,0,0,0),"cm"))
  }
  else
  {
    p[[count1]] <- ggplot(df_min_max_1, aes(x = x, y = y)) +
      # scale_y_continuous(breaks = 20,labels = 20)+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position="none",
            panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
      geom_polygon(data = df_min_max_1, aes(  group = id),fill = "light gray", alpha = 0.8)+
      geom_polygon(data = df_0025_1, aes( group = id),fill = "dark gray", alpha = 0.8)+
      geom_polygon(data = df_025_1, aes( group = id), fill = "gray27", alpha = 0.8)+
      geom_line(data = df_species1, aes(time, median), col = "black")+
      geom_hline(yintercept = log(line_value), linetype = "dashed")+ 
      # geom_hline(yintercept = log(K_real), linetype = "dotted")+
      coord_cartesian(ylim = c(0, 5))+
      theme(plot.margin=unit(c(0,0,0,0),"cm"))
  }
  count1 = count1+1
  }
}
# 
# par(mar=c(0.5, 0.5, 0.2, 0.2), mfrow=c(8,2),oma = c(4, 4, 0.2, 0.2))
# for(i in 1:16){
#   print(p[[i]])
# }
m = matrix(1:36,ncol = 4)
mu_vec = c(0,0.1,0.2,0.4)
Mig = c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000)



sg1 <- textGrob(substitute(paste(mu," = ",nn),list(nn = mu_vec[1])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
sg2 <- textGrob(substitute(paste(mu," = ",nn),list(nn = mu_vec[2])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

sg3 <- textGrob(substitute(paste(mu," = ",nn),list(nn = mu_vec[3])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

sg4 <- textGrob(substitute(paste(mu," = ",nn),list(nn = mu_vec[4])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
# 
# col_name = list()
# for(i in 1:8){
#   col_name[[i]]=textGrob(substitute(paste(M[0],"=",nn),list(nn = Mig[i])), gp=gpar(fontsize=15, fontface=3L))
# }

mg1 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[1])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg2 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[2])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg3 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[3])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg4 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[4])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg5 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[5])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg6 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[6])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg7 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[7])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg8 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[8])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
mg9 <- textGrob(substitute(paste("M"[0]," = ",nn),list(nn = Mig[9])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)
# g3 = arrangeGrob(grobs = col_name,ncol = 1)
g_stt4 = arrangeGrob(grobs = p, layout_matrix = m)
g_stt3 <- arrangeGrob(mg1,mg2,mg3,mg4,mg5,mg6,mg7,mg8,mg9,ncol = 1)
g_stt5 = textGrob("Time", gp=gpar(fontsize=x_title_fontsize),rot = 0)
g_stt1 = textGrob("")
g_stt2 <- arrangeGrob(sg1,sg2,sg3,sg4,ncol = 4)


file_add = paste0("~/Googlebox/Research/Project1/Revision/results/STT_S",Scenario,"L",locationtype,".pdf")
grid.arrange(g_stt1,g_stt2,g_stt1,label,g_stt4,g_stt3,g_stt1,g_stt5,g_stt1,ncol = 3,widths = c(1,16,3),heights = c(1,36,1))

# g = grid.arrange(g_stt1,g_stt2,g_stt1,label,g_stt4,g_stt3,g_stt1,g_stt5,g_stt1,ncol = 3,widths = c(1,16,3),heights = c(1,36,1))
# ggsave(file_add, plot = g)

# try(dev.off())



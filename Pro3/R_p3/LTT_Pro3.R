library(ggplot2)
library(gridExtra)
library(nLTT)
library(grid)
library(DDD)
source('C:/Liang/Googlebox/Research/Project3/Project3_R/event2L.R', echo=TRUE)

phi = c(0,0.001,0.01,0.1,1)
psi = phi

p = list()

Scenario = 2  

count1 = 1

age = 15
x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

 
# num_group = c(c(101:109),c(111:119),c(121:129),c(131:139))+(Scenario-1)*40
phi_group = c(1:5)
psi_group = c(1:5)
  i= 1
  j = 1
  mu_vec = rep(psi,each = 5)
  for(i in phi_group){
    for(j in psi_group){
    L3 = NULL
    fre = NULL
    L = list()
    data = NULL
    print(count1)
    for (seed in 1:100) {
      file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",i,"phi",j,"psi","sim",seed,".Rdata",sep = "")
      load(file = file)
      age = result$log$sT[length(result$log$sT)]
      
      if(result$events$ns[length(result$events$ns)] == 1){
        brts = c(-age,0)
        num_spe = c(1,1)
        data0 = cbind(i,brts,num_spe)
        data = rbind(data, data0)
      }else{
      L = event2L(result = result)
      tes = DDD::L2phylo(L,dropextinct = T)
      brts= -unname(sort(branching.times(tes),decreasing = T))
      if(brts[1] == -age){
        data0 = cbind(i,brts,c(2:(length(brts)+1)))
      }else{
        
        data0 = cbind(i,brts,c(2:(length(brts)+1)))
        data0 = rbind(c(i,-age,1),data0)
      }
     
      if(max(brts)<0) data0 = rbind(data0,c(i,0,data0[nrow(data0),3]))
      data = rbind(data, data0)
      }
    }
    
    # data <- data[-1,]
    time = data[order(data[,2]),2]
    timeu = unique(time)
    data_lineage = timeu
    for (seed in 1:100) {
      file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",i,"phi",j,"psi","sim",seed,".Rdata",sep = "")
      load(file = file)
      if(result$events$ns[length(result$events$ns)] == 1){
        data_lineage = cbind(data_lineage,1)
      }else{
      L = event2L(result = result)
      
      tes = DDD::L2phylo(L,dropextinct = T)
      brts= -unname(sort(branching.times(tes),decreasing = T))
      if(brts[1] == -age){
      M1 = match(brts,timeu)
      M1[1] = 1
      M11 = diff(M1)
      M13 = length(timeu)-max(M1)+1
      M12 = c(M11,M13)
      N1 = rep(2:(length(brts)+1),M12)
      data_lineage = cbind(data_lineage,N1)
      }else{
        brts = c(-age,brts)
        M1 = match(brts,timeu)
        M1[1] = 1
        M11 = diff(M1)
        M13 = length(timeu)-max(M1)+1
        M12 = c(M11,M13)
        N1 = rep(1:(length(brts)),M12)
        data_lineage = cbind(data_lineage,N1)
      }
      }
    }
    x = data_lineage[,1]
    z = data_lineage[,2:101]
    
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
   
    if(count1 %in% c(10,15,20,25)){
      
      p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
        theme(axis.title.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",
              panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.text.x=element_text(size = x_label_fontsize),axis.line.x = element_line(color="black"),plot.margin=unit(c(0,0,0,0),"cm"))+
        geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
        geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
        geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab(substitute(paste(psi,"=",nn),list(nn = mu_vec[count1])))+
        geom_line(data = df_lineage, aes(time, median), col = "black")+coord_cartesian(xlim=c(-age,0),ylim=c(0,log(80)))   #+ggtitle(substitute(paste(mu,"=",nn),list(nn = mu_vec[count1])))
     p[[count1]]=p[[count1]]+ scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))

    }
    
    else if(count1 %in% c(1:4)){
      
      p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
        theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
              panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.text.y=element_text(size = y_label_fontsize),axis.line.y = element_line(color="black"),plot.margin=unit(c(0,0,0,0),"cm"))+
        geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
        geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
        geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab(substitute(paste(psi,"=",nn),list(nn = mu_vec[count1])))+
        geom_line(data = df_lineage, aes(time, median), col = "black")+coord_cartesian(xlim=c(-age,0),ylim=c(0,log(80)))
        p[[count1]]=p[[count1]]+ scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))
    }
    else if(count1 == 5)
    {
      p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",
              panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.text.y=element_text(size = y_label_fontsize),axis.text.x=element_text(size = x_label_fontsize),axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"),plot.margin=unit(c(0,0,0,0),"cm"))+
        geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
        geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
        geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+
        geom_line(data = df_lineage, aes(time, median), col = "black")+coord_cartesian(xlim=c(-age,0),ylim=c(0,log(80)))
      p[[count1]]=p[[count1]]+ scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))
    }
    else
    {
      p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.text.x = element_blank(),legend.position="none",
              panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
        geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
        geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
        geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+
        geom_line(data = df_lineage, aes(time, median), col = "black")+coord_cartesian(xlim=c(-age,0),ylim=c(0,log(80)))
      p[[count1]]=p[[count1]]+ scale_y_continuous(breaks = c(log(2),log(5),log(10),log(20),log(40),log(80)),labels = c(2,5,10,20,40,80))
    }
    count1 = count1+1
    }
    }
  
  m = matrix(1:25,ncol = 5)

  
  
  mu_vec = psi
  # K_vec = c(20,40)
  Mig = phi
  
  
  
  sg1 <- textGrob(substitute(paste(phi," = ",nn),list(nn = mu_vec[1])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  sg2 <- textGrob(substitute(paste(phi," = ",nn),list(nn = mu_vec[2])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  
  sg3 <- textGrob(substitute(paste(phi," = ",nn),list(nn = mu_vec[3])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  
  sg4 <- textGrob(substitute(paste(phi," = ",nn),list(nn = mu_vec[4])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  sg5 <- textGrob(substitute(paste(phi," = ",nn),list(nn = mu_vec[5])), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
  
  # K_col1 <- textGrob(substitute(paste("K = ",nn),list(nn = K_vec[1])), gp=gpar(fontsize=15, fontface=3L))
  # K_col2 <- textGrob(substitute(paste("K = ",nn),list(nn = K_vec[2])), gp=gpar(fontsize=15, fontface=3L))
  
  mg1 <- textGrob(substitute(paste(psi," = ",nn),list(nn = Mig[1])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
  mg2 <- textGrob(substitute(paste(psi," = ",nn),list(nn = Mig[2])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
  mg3 <- textGrob(substitute(paste(psi," = ",nn),list(nn = Mig[3])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
  mg4 <- textGrob(substitute(paste(psi," = ",nn),list(nn = Mig[4])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
  mg5 <- textGrob(substitute(paste(psi," = ",nn),list(nn = Mig[5])), gp=gpar(fontsize=mig_title_fontsize, fontface=3L),hjust = 0,x = unit(0.1, "npc"))
  
  g_ltt2 <- arrangeGrob(sg1,sg2,sg3,sg4,sg5,ncol = 5)
  
  
  label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)
  g_ltt4 = arrangeGrob(grobs = p, layout_matrix = m)
  g_ltt3 <- arrangeGrob(mg1,mg2,mg3,mg4,mg5,ncol = 1)
  g_ltt5 = textGrob("Time", gp=gpar(fontsize=x_title_fontsize),rot = 0)
  g_ltt1 = textGrob("")
  g_a = textGrob("(a)")
  g_b = textGrob("(b)")
  
  
  file_add = paste0("C:/Liang/Googlebox/Research/Project3/Results/LTT_pro3.pdf")
  grid.arrange(g_ltt1,g_ltt2,g_ltt1,label,g_ltt4,g_ltt3,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,20,3),heights = c(1,20,1))
  # 
  # g = grid.arrange(g_ltt1,g_ltt2,g_ltt1,label,g_ltt4,g_ltt3,g_ltt1,g_ltt5,g_ltt1,ncol = 3,widths = c(1,16,1),heights = c(1,36,1))
  # ggsave(file_add, plot = g)
  # C:\Liang\Googlebox\Research\Project3\Results
  
  
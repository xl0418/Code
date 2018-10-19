library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
source('~/Googlebox/Research/Project1/R_pro1/g_legend.R', echo=TRUE)
estplot = function(sce){
plottype = 4
plotl_est = list()


x_label_fontsize = 10
y_label_fontsize = 10
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

tens1 = 1
tens2 = 0


if(sce == 1)  tens1_vec = c(1:4) # Scenario 1
if(sce == 2)  tens1_vec = c(5:8)  # Scenario 2
if(sce == 3)  tens1_vec = c(9:12)  # Scenario 3

for(tens in tens1_vec){
  tens2 = tens2+1
  
  K = NULL
  lambda = NULL
  mu = NULL
  j = 2
  N_tips = NULL
  endmc = 1000
  # par(mfrow = c(2,3),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
  lambda_vec = 0.8
  
  mu_vec = rep(c(0,0.1,0.2,0.4),each = 3)
  
  K_vec = rep(c(20,40),each = 4)
  M0_vec = rep(c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000),each = 4)
  
  
  a = 100+10*(tens-1)+1
  b = 100+10*(tens-1)+9
  g20g = c(221:224)
  g40g = c(231:234)
  g20 = g20g[tens2]
  g40 = g40g[tens2]
  
  for(num in c(c(a:b),g20,g40) ){
    file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    
    
    if(num %in% g20g){
      file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    }
    if(num %in% g40g){
      file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
      
    }
    
    # file1 = paste(getwd(),"/Desktop/data/Modeltestgroup",j,sep = "")
    print(num)
    
    
    
    for(i in 1:100){
      # print(i)
      file = paste(file1,"/",num,"analysis",i,".Rdata",sep = "")
      load(file = file)
      # print(paste("pvalue",i," = ",pvalue,"; poweroftest",i," = ",poweroftest,sep = "")) 
      
      lambda1 = out$lambda_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$lambda_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$lambda_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      mu1 = out$mu_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$mu_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$mu_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      K1 = 1E+120 * (opt == 1) + pmin(1E+120,out$K_DD1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K_DD2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
      K_est = lambda1*K1/(lambda1-mu1)
      K = c(K,K_est)
      lambda = c(lambda,lambda1)
      mu = c(mu, mu1)
      
    }
    
  }
  
  
  
  M0 = rep(c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000,1020,1040),each = 100000)
  
  p_K1 = cbind(M0,K)
  
  p_K3 = matrix(p_K1, ncol = 2)
  p_K4 = as.data.frame(p_K3)
  colnames(p_K4) = c("M0", "K")
  
  p_lambda1 = cbind(M0, lambda)
  p_lambda3 = matrix(p_lambda1, ncol = 2)
  
  p_lambda4 = as.data.frame(p_lambda3)
  colnames(p_lambda4) = c("M0", "lambda")
  
  p_mu1 = cbind(M0, mu)
  p_mu3 = matrix(p_mu1, ncol = 2)
  p_mu4 = as.data.frame(p_mu3)
  colnames(p_mu4) = c("M0", "mu")
 
  p_K4_1 = round(p_K4,2)
  p_K4_1 = p_K4_1[p_K4_1$K<100,]
  p_lambda4_1 = round(p_lambda4,2)
  p_lambda4_1 = p_lambda4_1[p_lambda4_1$lambda < 20,]
  
  
  p_mu4_1 = round(p_mu4,2)
  p_mu4_1 = p_mu4_1[p_mu4_1$mu < 20,]
  batch = c(rep(c(1),9),2,3)
  
  #Real quantiles
  K_quantile = aggregate(K ~ M0, p_K4_1,quantile, probs = c(0.05,0.25,0.5,0.75,0.95))
  lambda_quantile = aggregate(lambda ~ M0, p_lambda4_1,quantile, probs = c(0.05,0.25,0.5,0.75,0.95))
  mu_quantile = aggregate(mu ~ M0, p_mu4_1,quantile, probs = c(0.05,0.25,0.5,0.75,0.95))
  
  K_quantile = cbind(K_quantile, batch)
  lambda_quantile = cbind(lambda_quantile, batch)
  mu_quantile = cbind(mu_quantile, batch)
  
  
  #for K est.
  if(tens1 == 1){
    plotl_est[[tens1]] <-ggplot(K_quantile,aes(x=factor(M0), ymin=K[,1], lower=K[,2],middle=K[,3], upper=K[,4],ymax=K[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
      xlab("")+ylab("")+ 
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5),plot.title = element_text(size = 16),axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = sort(c(seq(0, 80, length.out=5), c(20,40))))+
      scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
    if(plottype == 4) { plotl_est[[tens1]]=  plotl_est[[tens1]]+ ggtitle(substitute(paste(mu," = ",nn),list(nn = mu_vec[tens1])))#+ geom_hline(aes(yintercept=20),linetype = "dashed",color = "black") #+geom_hline(aes(yintercept=40),linetype = "dashed",color = "black")
    } else if(plottype == 2) {plotl_est[[tens1]]=  plotl_est[[tens1]]+ ggtitle(substitute(paste("K = ",nn),list(nn = K_vec[tens1])))+ geom_hline(aes(yintercept=40),linetype = "dashed",color = "black")+geom_hline(aes(yintercept=80),linetype = "dashed",color = "black")}
    
    
  }
  else {
    plotl_est[[tens1]] <- ggplot(K_quantile,aes(x=factor(M0), ymin=K[,1], lower=K[,2],middle=K[,3], upper=K[,4],ymax=K[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
      xlab("")+ylab("")+ 
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5),plot.title = element_text(size = 16),axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = sort(c(seq(0, 80, length.out=5), c(20,40))))+
    scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
    if(plottype == 4) { plotl_est[[tens1]]=  plotl_est[[tens1]]+ ggtitle(substitute(paste(mu," = ",nn),list(nn = mu_vec[tens1])))#+ geom_hline(aes(yintercept=20),linetype = "dashed",color = "black") #+geom_hline(aes(yintercept=40),linetype = "dashed",color = "grey")
    } else if(plottype == 2) {plotl_est[[tens1]]=  plotl_est[[tens1]]+ ggtitle(substitute(paste("K = ",nn),list(nn = K_vec[tens1])))+ geom_hline(aes(yintercept=40),linetype = "dashed",color = "black")+geom_hline(aes(yintercept=80),linetype = "dashed",color = "grey")}
    
    
  }
  
  
  # for lambda est.
  tens1 = tens1+1
  if(tens1 == 2){
    plotl_est[[tens1]] <- ggplot(lambda_quantile,aes(x=factor(M0), ymin=lambda[,1], lower=lambda[,2],middle=lambda[,3], upper=lambda[,4],ymax=lambda[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(), axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5),axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      xlab("")+ylab("")+geom_hline(yintercept = 0.8,linetype = "dashed",color = "black")+coord_cartesian(ylim = c(0, 3))+scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3),labels = c(0,"",1,"",2,"",3) )+
    scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
  }
  else{
    plotl_est[[tens1]] <- ggplot(lambda_quantile,aes(x=factor(M0), ymin=lambda[,1], lower=lambda[,2],middle=lambda[,3], upper=lambda[,4],ymax=lambda[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5),axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      xlab("")+ylab("")+geom_hline(yintercept = 0.8,linetype = "dashed",color = "black")+coord_cartesian(ylim = c(0, 3))+scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3),labels = c(0,"",1,"",2,"",3) )+
    scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
  }
  
  
  #for mu est.
  tens1 = tens1+1
  if(tens1 == 3){
    plotl_est[[tens1]] <-ggplot(mu_quantile,aes(x=factor(M0), ymin=mu[,1], lower=mu[,2],middle=mu[,3], upper=mu[,4],ymax=mu[,5],fill = factor(batch),color = factor(batch))) + 
      geom_boxplot(stat="identity")+xlab("")+
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5), axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      ylab("")+geom_hline(yintercept = mu_vec[tens1],linetype = "dashed",color = "black")+coord_cartesian(ylim = c(0, 1))+scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,"",0.4,"",0.8,"") )+
    scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
  }
  else{
    plotl_est[[tens1]] <- ggplot(mu_quantile,aes(x=factor(M0), ymin=mu[,1], lower=mu[,2],middle=mu[,3], upper=mu[,4],ymax=mu[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
      xlab("")+
      theme(axis.text.x=element_text(size = x_label_fontsize,angle=90,hjust=1,vjust=0.5),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize,angle=90,hjust=0.5,vjust=0.5), axis.line.y = element_line(color="black"),panel.background = element_blank(),legend.position="none")+
      ylab("")+
      coord_cartesian(ylim = c(0, 1))+scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,"",0.4,"",0.8,"") )+
    scale_x_discrete(labels=c("1020" = "0","1040" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
    if(plottype == 4) { plotl_est[[tens1]]=  plotl_est[[tens1]]+ geom_hline(yintercept = mu_vec[tens1],linetype = "dashed",color = "black")
    } else if(plottype == 2) {plotl_est[[tens1]]=  plotl_est[[tens1]]+  geom_hline(yintercept = 0,linetype = "dashed",color = "black")}
    
  }
  
  tens1 = tens1+1
  
  
}

legend_plot = ggplot(mu_quantile,aes(x=factor(M0), ymin=mu[,1], lower=mu[,2],middle=mu[,3], upper=mu[,4],ymax=mu[,5],fill = factor(batch),color = factor(batch))) + geom_boxplot(stat="identity")+
  theme(legend.position="right")+
  scale_fill_manual(guide=FALSE,values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+
  scale_color_manual(name = "Scenarios",labels = c("S1", "S4", "S5"),values=c("#2171B5", "#EF6548","#DE77AE"))+
  guides(color=guide_legend(override.aes=list(fill=c("#9ECAE1", "#FC8D59","#FDE0EF"))))
mylegend = g_legend(legend_plot)

lambda_grob <- textGrob(substitute(paste(lambda[0])), gp=gpar(fontsize=16),rot = 90)
mu_grob <- textGrob(substitute(paste(mu)), gp=gpar(fontsize=16),rot = 90)
K_grob <- textGrob(substitute(paste("K'")), gp=gpar(fontsize=16),rot = 90)

m = matrix(1:12,3,4)

grob1 <- arrangeGrob(K_grob, lambda_grob,mu_grob,ncol = 1)
grob2 = arrangeGrob(grobs = plotl_est, layout_matrix = m)
grob4 = textGrob("Dispersal rate", gp=gpar(fontsize=16),rot = 0)
grob3 = textGrob("")


grid.arrange(grob1,grob2,mylegend,grob3,grob4,grob3,ncol = 3,widths = c(1,16,2.5),heights = c(16,1))

}


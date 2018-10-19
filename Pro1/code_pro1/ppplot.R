library(ggplot2)
library(grid)
library("RColorBrewer")
source(paste0(getwd(),'/g_legend.R'))
ppplot = function(sce){
plotl = list()


x_label_fontsize = 10
y_label_fontsize = 10
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16

if(sce == 1)  tens1_vec = c(1:4) # Scenario 1
if(sce == 2)  tens1_vec = c(5:8)  # Scenario 2
if(sce == 3)  tens1_vec = c(9:12)  # Scenario 3

tens1 = 1
tens2 = 0
yb = 0.05
p = NULL
power = NULL
g40g = c(231:234)
g20g = c(221:224)
for(tens in tens1_vec){
  tens2 = tens2+1
  j = 2
  p = NA
  power = NA
  
  g40 = g40g[tens2]
  g20 = g20g[tens2]
  a = 100+10*(tens-1)+1
  b = 100+10*(tens-1)+9
  
  for(num in c(c(a:b),g20,g40) ){
  
    file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
    
    if(num %in% g20g){
      file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
      
    }
    if(num %in% g40g){
      file1 = paste(getwd(),j,"Modeltestgroup15age",num,sep = "")
      
    }
    
    print(num)
    
    for(i in 1:100){
      
      file = paste(file1,"/",num,"analysis",i,".Rdata",sep = "")
      
      load(file = file)
    
      p <- c(p,pvalue)
      power <- c(power,poweroftest)
      
      
      
    }
  
  }
  
  
  group = rep(c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000,10000,10020),each = 100)
  batch = c(rep(c(1),900),rep(c(2),100),rep(c(3),100))
  
  p2 = cbind(group, p , power,batch)
  p2 = as.data.frame(p2)
  
  mu_vec = rep(c(0,0.1,0.2,0.4),each = 2)
  
  
  if(tens1 == 1){
    plotl[[tens1]] <- ggplot(p2, aes(x=factor(group), y=p,fill = factor(batch),color = factor(batch) )) +geom_boxplot(outlier.shape=NA)+xlab("")+ylab("p-value")+#geom_hline(yintercept = 0.05)
      geom_hline(aes(yintercept=yb),linetype = "dashed") +coord_cartesian(ylim = c(0, 1)) + 
      theme(axis.text.x=element_text(size = x_label_fontsize, angle=90,hjust=1,vjust=0.5),legend.position="none",panel.background = element_blank(),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),plot.title = element_text(size = 16))+
      scale_y_continuous(breaks = sort(c(seq(0,1, length.out=5), yb)))+scale_x_discrete(labels=c("10000" = "0","10020" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))+ 
      ggtitle(substitute(paste(mu," = ",nn),list(nn = mu_vec[tens1])))
  }
  else {
    plotl[[tens1]] <- ggplot(p2, aes(x=factor(group), y=p,fill = factor(batch),color = factor(batch) )) +geom_boxplot(outlier.shape=NA)+xlab("")+ylab("")+#geom_hline(yintercept = 0.05)
      geom_hline(aes(yintercept=yb),linetype = "dashed") +coord_cartesian(ylim = c(0, 1)) +  
      theme(axis.text.x=element_text(size = x_label_fontsize, angle=90,hjust=1,vjust=0.5),legend.position="none",panel.background = element_blank(),
            axis.title=element_blank(), axis.text.y=element_text(size = y_label_fontsize), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),plot.title = element_text(size = 16))+
      scale_y_continuous(breaks = sort(c(seq(0,1, length.out=5), yb)))+scale_x_discrete(labels=c("10000" = "0","10020" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))+
      ggtitle(substitute(paste(mu," = ",nn),list(nn = mu_vec[tens1])))
    
    
  } 
  
  
  tens1 = tens1 + 1
  if(tens1 == 2){
    plotl[[tens1]] <- ggplot(p2, aes(x=factor(group), y=power,fill = factor(batch),color = factor(batch) )) +geom_boxplot(outlier.shape=NA)+xlab("")+ylab("power")+coord_cartesian(ylim = c(0, 1))+
      theme(axis.text.x=element_text(size = x_label_fontsize, angle=90,hjust=1,vjust=0.5),legend.position="none",panel.background = element_blank(),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),plot.title = element_text(size = 16))+scale_x_discrete(labels=c("10000" = "0","10020" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
  }
  
  else {
    plotl[[tens1]] <- ggplot(p2, aes(x=factor(group), y=power,fill = factor(batch),color = factor(batch) )) +geom_boxplot(outlier.shape=NA)+xlab("")+ylab("")+coord_cartesian(ylim = c(0, 1))+  
      theme(axis.text.x=element_text(size = x_label_fontsize, angle=90,hjust=1,vjust=0.5),legend.position="none",panel.background = element_blank(),
            axis.title=element_blank(),axis.text.y=element_text(size = y_label_fontsize),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"))+scale_x_discrete(labels=c("10000" = "0","10020" = "0"))+scale_fill_manual(values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+scale_color_manual(values=c("#2171B5", "#EF6548","#DE77AE"))
  }
  
  tens1 = tens1 + 1
  
  
}

legend_plot = ggplot(p2, aes(x=factor(group), y=power,fill = factor(batch),color = factor(batch) )) +geom_boxplot(outlier.shape=NA)+xlab("")+ylab("")+coord_cartesian(ylim = c(0, 1))+  
  theme(legend.position="right")+
  scale_fill_manual(guide=FALSE,values=c("#9ECAE1", "#FC8D59","#FDE0EF"))+
  scale_color_manual(name = "Scenarios",labels = c("S3", "S4", "S5"),values=c("#2171B5", "#EF6548","#DE77AE"))+
  guides(color=guide_legend(override.aes=list(fill=c("#9ECAE1", "#FC8D59","#FDE0EF"))))
mylegend = g_legend(legend_plot)


pvalue_grob <- textGrob("p-value", gp=gpar(fontsize=16),rot = 90)
power_grob <- textGrob("Power", gp=gpar(fontsize=16),rot = 90)

m_pp = matrix(1:8,2,4)
grob_pp1 <- arrangeGrob(pvalue_grob,power_grob,ncol = 1)
grob_pp2 = arrangeGrob(grobs = plotl, layout_matrix = m_pp)
grob_pp4 = textGrob("Dispersal rate", gp=gpar(fontsize=16),rot = 0)
grob_pp3 = textGrob("")


grid.arrange(grob_pp1,grob_pp2,mylegend,grob_pp3,grob_pp4,grob_pp3,ncol = 3,widths = c(1,16,2),heights = c(8,1))


}


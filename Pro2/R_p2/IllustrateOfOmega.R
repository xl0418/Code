library(ggthemes)
library(DDD)
library(ggtree)


omega <- function(deltaz,alpha){
  contri <- alpha*deltaz*exp(-alpha*deltaz^2)
  return(contri)
}


omega.contribution = c()
alpha.range = c(0.01,0.1,0.5,1)  #seq(0,1,0.01)
deltaz.range = seq(-10,10,0.1)
for(alpha in alpha.range){
  omega.temp <-  omega(deltaz = deltaz.range, alpha = alpha)
  omega.data <-  cbind(deltaz.range,omega.temp,alpha)
  omega.contribution <- rbind(omega.contribution,omega.data)
}

omega.con.df = as.data.frame(omega.contribution)
names(omega.con.df) <- c('deltaz','values','alpha')
linetype = rev(c('solid','dotdash','dashed','dotted'))

omega.plot <- ggplot(data=omega.con.df, aes(x=deltaz, y=values,group=alpha))+
  geom_line(aes(linetype=as.factor(alpha),color = as.factor(alpha)),size = 1.2)+
  scale_color_manual(values=rev(c("#FEBC44", "#EF6174", "#373150", "#629182")),name = expression(alpha))+
  scale_linetype_manual(values=linetype,name = expression(alpha))+
  theme_bw()+ xlab(expression(mu[i]-mu[j]))+ylab(expression(paste('Competitive repulsion ',alpha,'(',mu[i]-mu[j],')',e^{-alpha(mu[i]-mu[j])^2})))+
  theme(text = element_text(size=20),legend.position = "none")+
  scale_x_continuous( breaks = seq(-10,10,5))+
  annotate("text", x = 2, y = 0.35, label = expression(alpha==1),size = 6)+
  annotate("text", x = 3, y = 0.22, label = expression(alpha==0.5),size = 6)+
  annotate("text", x = 4, y = 0.15, label = expression(alpha==0.1),size = 6)+
  annotate("text", x = 8, y = 0.08, label = expression(alpha==0.01),size = 6)
  
phi.plot.name = 'C:\\Liang\\Googlebox\\Research\\Project2\\BaleenWhales\\result_cluster\\results_ms_posterior\\OMEGA.pdf'
ggsave(omega.plot,width = 8,height = 6)

L.example = matrix(0, 4,4)
L.example[,1] = c(20,20,15,5)
L.example[,2] = c(0,1,2,3)
L.example[,3] = c(1:4)
L.example[,4] = rep(-1,4)
L.example[,1] = L.example[,1] * 100
dd <- data.frame(taxa  = c("Species 1 (400)",'Species 2 (100)','Species 3 (10)','Species 4 (800)'),
                 country = c(1:4))
row.names(dd) <- NULL

phy.example = DDD::L2phylo(L.example)
phy.example$tip.label = c("Species 1 (400)",'Species 2 (100)','Species 3 (10)','Species 4 (800)')
tree.plot = ggtree(phy.example)+theme_tree2()+geom_tiplab(size = 5,hjust = -0.2)+
  ggplot2::xlim(0, 2500) + ggplot2::xlab('Time') + theme(text = element_text(size=20))
p <- tree.plot %<+% dd + geom_tippoint(aes(color=factor(country)),size = 10,alpha=1)+
  scale_colour_manual(values = c("#D0104C", "#FFC408", "#0089A7", "#227D51"))

tree.name = 'C:\\Liang\\Googlebox\\Research\\Project3\\replicate_sim_9sces_results\\treeplot.pdf'
ggsave(tree.name,width = 8,height = 6)



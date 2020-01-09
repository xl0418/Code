library(ggthemes)
library(DDD)
library(ggtree)
abundance.vec = c(800,10,100,400)
phylogenetic.matrix = matrix(c(0,5,15,20,
                               5,0,15,20,
                               15,15,0,20,
                               20,20,20,0),nrow = 4)
phi = 1


phylo.spatial.effect <- function(abundance,phymatrix,phi){
  contri <- exp(-phi*phymatrix^2)%*%abundance / sum(abundance)
}

phi.contribution = c()
phi.range = seq(0,0.1,0.001)

for(phi in phi.range){
  phi.temp <-  phylo.spatial.effect(abundance = abundance.vec,phymatrix=phylogenetic.matrix,
                                    phi = phi)
  phi.contribution <- rbind(phi.contribution,t(phi.temp))
}

phi.con.df = data.frame()
species.index = rev(c(1:4))
for(i in c(1:ncol(phi.contribution))){
  phi.con.df = rbind(phi.con.df,cbind(phi.range,phi.contribution[,i],species.index[i]))
}
names(phi.con.df) <- c('phi','values','species')
linetype = c('solid', 'dashed','twodash','dotted')
phi.plot <- ggplot(data=phi.con.df, aes(x=phi, y=values,group=species))+
            geom_line(aes(linetype=as.factor(species),color = as.factor(species)),size = 1.2)+
  scale_color_manual(values=c("#D0104C", "#FFC408", "#0089A7", "#227D51"),name = "Species")+
  scale_linetype_manual(values=linetype,name = "Species")+
  theme_bw()+ xlab(expression(phi))+ylab('CNDD contribution')+theme(text = element_text(size=20))

phi.plot.name = 'C:\\Liang\\Googlebox\\Research\\Project3\\replicate_sim_9sces_results\\CNDD.pdf'
ggsave(phi.plot.name,width = 8,height = 6)

L.example = matrix(0, 4,4)
L.example[,1] = c(20,20,15,5)
L.example[,2] = c(0,1,2,3)
L.example[,3] = c(1:4)
L.example[,4] = rep(-1,4)
dd <- data.frame(taxa  = c("Species 1 (400)",'Species 2 (100)','Species 3 (10)','Species 4 (800)'),
                 country = c(1:4))
row.names(dd) <- NULL

phy.example = DDD::L2phylo(L.example)
phy.example$tip.label = c("Species 1 (400)",'Species 2 (100)','Species 3 (10)','Species 4 (800)')
tree.plot = ggtree(phy.example)+theme_tree2()+geom_tiplab(size = 5,hjust = -0.2)+
  ggplot2::xlim(0, 25) + ggplot2::xlab('Time') + theme(text = element_text(size=20))
p <- tree.plot %<+% dd + geom_tippoint(aes(color=factor(country)),size = 10,alpha=1)+
  scale_colour_manual(values = c("#D0104C", "#FFC408", "#0089A7", "#227D51"))

tree.name = 'C:\\Liang\\Googlebox\\Research\\Project3\\replicate_sim_9sces_results\\treeplot.pdf'
ggsave(tree.name,width = 8,height = 6)



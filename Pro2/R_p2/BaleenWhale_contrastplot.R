library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)

source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/Est/'

emdata = read.nexus(emdatadir)

dropextinct = T
baleenwhale = phylo2L(emdata,error = 1e-5)
L_ext = baleenwhale$L
extantspecieslabel = baleenwhale$ESL

extantspecieslabel <- c("B.mysticetus", "E.australis"       
                        , "E.glacialis", "E.japonica"        
                        , "B.acutorostrata" ,"B.bonaerensis"  
                        , "B.borealis" ,"B.brydei"       
                        , "B.edeni" , "B.omurai"       
                        , "B.musculus" , "B.physalus"     
                        , "M.novaeangliae" , "E.robustus"     
                        , "C.marginata" )
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
phylo_test$tip.label <- extantspecieslabel
phylo_test$node.label <- c(1:14)
# plot(phylo_test,show.tip.label = TRUE,show.node.label = TRUE)

# Extract the order names of simulated species
prunedL = phylo2L(phylo_test)
pl = prunedL$L
phylo_pl = DDD::L2phylo(pl)
# plot(phylo_pl,show.tip.label = TRUE,show.node.label = TRUE)

rank.sim.species <- as.numeric(gsub("t","", phylo_pl$tip.label))
sim.species.col.names = c()
for(i in c(1:15)){
  sim.species.col.names[rank.sim.species[i]] <- extantspecieslabel[i]
}




# empirical data
fileemp_name = paste0(dir,'emp.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]


obsZ_pic = 10^obsZ_emp[,2]
names(obsZ_pic) = extantspecieslabel

# empirical contrast
empirical_pic<-pic(obsZ_pic,phylo_test)


# calculate sim contrast
mode.label = c('TVP','TV','TVM')


simfile_TVP = paste0(dir,'predictsimTVP.csv')
simfile_TV = paste0(dir,'predictsimTV.csv')
simfile_TVM = paste0(dir,'predictsimTVM.csv')


# predict simulations
predictZ_TVP = read.csv(simfile_TVP)
predictZ_TV = read.csv(simfile_TV)
predictZ_TVM = read.csv(simfile_TVM)




sort = 0

if(sort == 0){
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)

  dimnames(predictZ_matrix_TVP)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TV)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TVM)[[2]] = sim.species.col.names
  
  
}else{
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  TVP_Z_order = t(apply(predictZ_matrix_TVP,1,order))
  TVM_Z_order = t(apply(predictZ_matrix_TVM,1,order))
  
  predictZ_matrix_TVP=t(apply(predictZ_matrix_TVP,1,sort))
  predictZ_matrix_TV=t(apply(predictZ_matrix_TV,1,sort))
  predictZ_matrix_TVM=t(apply(predictZ_matrix_TVM,1,sort))

  dimnames(predictZ_matrix_TVP)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TV)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TVM)[[2]] = sorted.species.labels
  
}


samplesize = nrow(predictZ_matrix_TVP)


TVP.pic <- c()
TV.pic <- c()
TVM.pic <- c()

scaled = TRUE

for(size in c(1:samplesize)){
  TVP.pic <- rbind(TVP.pic,pic(predictZ_matrix_TVP[size,],phylo_test, scaled = scaled))
  TV.pic <- rbind(TV.pic,pic(predictZ_matrix_TV[size,],phylo_test, scaled = scaled))
  TVM.pic <- rbind(TVM.pic,pic(predictZ_matrix_TVM[size,],phylo_test, scaled = scaled))
  
}
restructured_col <- c(2,3,4,1,5,7,6,13,14,8,9,10,11,12)
TVP.pic <- TVP.pic[,restructured_col]
TV.pic <- TV.pic[,restructured_col]
TVM.pic <- TVM.pic[,restructured_col]



d_all_TVP = as.data.frame(as.table(TVP.pic))[,2:3]
colnames(d_all_TVP) = c('betspecies','contrast')
d_all_TV = as.data.frame(as.table(TV.pic))[,2:3]
colnames(d_all_TV) = c('betspecies','contrast')
d_all_TVM = as.data.frame(as.table(TVM.pic))[,2:3]
colnames(d_all_TVM) = c('betspecies','contrast')



d_meanemp = data.frame(betspecies=c(1:14), contrast=empirical_pic)


plot_tree <- ggtree(phylo_test,size=.5) +geom_tiplab(size=3.5,fontface="bold") #+xlim(0,80)


plot_sepboxplt_TVP <- facet_plot(plot_tree, panel="TVP", data=d_all_TVP, geom_boxploth, 
                                 mapping = aes(x=contrast, group=label),color='#3ac569',fill= '#cff0da',outlier.colour = NULL)  + theme_tree2()

p_finalTVP <- facet_plot(plot_sepboxplt_TVP+xlim_tree(40), panel="TVP", data=d_meanemp, geom_point, 
                         mapping = aes(x=contrast, group=label ),color = 'black')


plot_sepboxplt_TV <- facet_plot(p_finalTVP, panel="TV", data=d_all_TV, geom_boxploth, 
                                mapping = aes(x=contrast, group=label ),color='#f9320c',fill='#f1bbba',outlier.colour = NULL )  + theme_tree2()

p_finalTV <- facet_plot(plot_sepboxplt_TV+xlim_tree(40), panel="TV", data=d_meanemp, geom_point, 
                        mapping = aes(x=contrast, group=label ),color = 'black')


plot_sepboxplt_TVM <- facet_plot(p_finalTV, panel="TVM", data=d_all_TVM, geom_boxploth, 
                                 mapping = aes(x=contrast, group=label),color='#6a60a9',fill='#dedcee' ,outlier.colour = NULL )  + theme_tree2()

p_finalTVM <- facet_plot(plot_sepboxplt_TVM+xlim_tree(40), panel="TVM", data=d_meanemp, geom_point, 
                         mapping = aes(x=contrast, group=label ),color = 'black')+
  theme(strip.background = element_rect(colour="white", fill="white"), 
        strip.text.x = element_text(size=12, face="bold"),
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))


p_finalTVM
lbs <- c(Tree = "Phylogenetic tree\nof baleen whales", TVP = "Trait evolution \n+ population dynamics",
         TV = "Trait evolution", TVM ="Trait evolution \n+ metabolism dynamics")
facet_labeller(p_finalTVM, lbs)


savefile = paste0(dir,'predictcontrastimage_TVP_TV_TVM_sorted',count,'.png')
ggsave(savefile,p_finalTVM)




# phylomorphospace(phylo_test,cbind(obsZ_pic,predictZ_matrix_TVP[1,]),label="off",node.size=c(0,1))





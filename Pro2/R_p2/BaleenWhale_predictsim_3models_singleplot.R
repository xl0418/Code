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
plot(phylo_test,show.tip.label = TRUE)


# empirical data
fileemp_name = paste0(dir,'emp.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]

mode.label = c('TVP','TV','TVM')


simfile_TVP = paste0(dir,'predictsimTVP.csv')
simfile_TV = paste0(dir,'predictsimTV.csv')
simfile_TVM = paste0(dir,'predictsimTVM.csv')

# predict simulations
predictZ_TVP = read.csv(simfile_TVP)
predictZ_TV = read.csv(simfile_TV)
predictZ_TVM = read.csv(simfile_TVM)


sort = 1

if(sort == 0){
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  
  
}else{
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  predictZ_matrix_TVP=t(apply(predictZ_matrix_TVP,1,sort))
  predictZ_matrix_TV=t(apply(predictZ_matrix_TV,1,sort))
  predictZ_matrix_TVM=t(apply(predictZ_matrix_TVM,1,sort))
}

samplesize = nrow(predictZ_matrix_TVP)
dimnames(predictZ_matrix_TVP)[[2]] = sorted.species.labels
dimnames(predictZ_matrix_TV)[[2]] = sorted.species.labels
dimnames(predictZ_matrix_TVM)[[2]] = sorted.species.labels


d_all_TVP = as.data.frame(as.table(predictZ_matrix_TVP))[,2:3]
colnames(d_all_TVP) = c('species','traitall')
d_all_TV = as.data.frame(as.table(predictZ_matrix_TV))[,2:3]
colnames(d_all_TV) = c('species','traitall')
d_all_TVM = as.data.frame(as.table(predictZ_matrix_TVM))[,2:3]
colnames(d_all_TVM) = c('species','traitall')


obsZ_mean = 10^(obsZ_emp[,2])



d_meanemp = data.frame(species=species_label, trait=obsZ_mean)


plot_tree <- ggtree(phylo_test)+geom_tiplab(size=3.5) #+xlim(0,80)


plot_sepboxplt_TVP <- facet_plot(plot_tree, panel="TL of TVP", data=d_all_TVP, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#99CCFF"))

p_finalTVP <- facet_plot(plot_sepboxplt_TVP+xlim_tree(40), panel="TL of TVP", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'red')


plot_sepboxplt_TV <- facet_plot(p_finalTVP, panel="TL of TV", data=d_all_TV, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#99CCFF"))

p_finalTV <- facet_plot(plot_sepboxplt_TV+xlim_tree(40), panel="TL of TV", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'red')


plot_sepboxplt_TVM <- facet_plot(p_finalTV, panel="TL of TVM", data=d_all_TVM, geom_boxploth, 
                                mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#99CCFF"))

p_finalTVM <- facet_plot(plot_sepboxplt_TVM+xlim_tree(40), panel="TL of TVM", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'red')


# p_finalTP <- p_finalTP+ggtitle(count)+
#   theme(plot.title = element_text(hjust = 0.5))


p_finalTVM

savefile = paste0(dir,'predictimage_TVP_TV_TVM',count,'.png')
ggsave(savefile,p_finalTVM)










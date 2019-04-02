library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)
os = Sys.info()['sysname']
if(os == 'Darwin'){
  source('~/Documents/GitHub/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('~/Documents/GitHub/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
  dir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/'
  
}else{
  source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
  dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'
}
emdata = read.nexus(emdatadir)

dropextinct = T
baleenwhale = phylo2L(emdata,error = 1e-5)
L_ext = baleenwhale$L
extantspecieslabel = baleenwhale$ESL
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
phylo_test$tip.label <- extantspecieslabel
plot(phylo_test,show.tip.label = TRUE)

# empirical data
fileemp_name = paste0(dir,'trait_labels.csv')
obsZ_emp = read.csv(fileemp_name)
species_label = obsZ_emp[,1]
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]


filepredict_name =  paste0(dir,'predictsim.csv')
predictZ = read.csv(filepredict_name)



predictZ_matrix = as.matrix(predictZ)
predictZ_mean = colMeans(predictZ_matrix)

samplesize = nrow(predictZ_matrix)
dimnames(predictZ_matrix)[[2]] = sorted.species.labels

d_all = as.data.frame(as.table(predictZ_matrix))[,2:3]
colnames(d_all) = c('species','traitall')

obsZ_mean = obsZ_emp[,2]

d_meanemp = data.frame(species=species_label, trait=obsZ_mean)
d_meansim = data.frame(species=species_label, traitdot=predictZ_mean)


plot_tree <- ggtree(phylo_test)+geom_tiplab(size=2.5)

# plot_dottips = plot_tree %<+% d_meanemp+ geom_tippoint(aes(size=trait))




plot_sepdots = facet_plot(plot_tree, panel="dot", data=d_meanemp, geom=geom_point, aes(x=trait), color='firebrick')+ theme_tree2()

plot_sepboxplt <- facet_plot(plot_tree, panel="Trait", data=d_all, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#D1B6E1"))

p_final <- facet_plot(plot_sepboxplt+xlim_tree(35), panel="Trait", data=d_meanemp, geom_point, 
           mapping = aes(x=trait, group=label ),color = 'red')

p_final










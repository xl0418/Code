library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
os = Sys.info()['sysname']
if(os == 'Darwin'){
  source('~/Documents/GitHub/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('~/Documents/GitHub/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
  dir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/'
  
}else{
source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
dir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/'
}
emdata = read.tree(emdatadir)
# prune tree by phytools
phy_prune = fancyTree(emdata, type="droptip",tip = getExtinct(emdata),cex = 0.7)

# trait data 
setwd('C:/Liang/Code/Pro2/data')
fileunsort_name = 'unsort.csv'
obsZ_unsort = read.csv(fileunsort_name)
obsZ_matrix = as.matrix(obsZ_unsort)
obsZ_mean = colMeans(obsZ_matrix)
species_label = phy_prune$tip.label
samplesize = nrow(obsZ_matrix)
dimnames(obsZ_matrix)[[2]] = species_label


d_all = as.data.frame(as.table(obsZ_matrix))[,2:3]
colnames(d_all) = c('species','traitall')

nor_trait = (obsZ_mean- range(obsZ_mean)[1])/(range(obsZ_mean)[2]-range(obsZ_mean)[1])

d_mean1 = data.frame(species=phy_prune$tip.label, trait=obsZ_mean)
d_mean2 = data.frame(species=phy_prune$tip.label, traitdot=obsZ_mean)

plot_tree <- ggtree(phy_prune)

plot_dottips = plot_tree %<+% d_mean1+ geom_tippoint(aes(size=trait))

plot_sepdots = facet_plot(plot_dottips, panel="dot", data=d_mean2, geom=geom_point, aes(x=traitdot), color='firebrick')+ theme_tree2()

plot_sepboxplt <- facet_plot(plot_sepdots, panel="2", data=d_all, geom_boxploth, 
                 mapping = aes(x=traitall, group=label,color = traitall))  + theme_tree2()

plot_sepboxplt






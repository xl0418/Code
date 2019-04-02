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
fileemp_name = paste0(dir,'trait_labels.csv')
obsZ_emp = read.csv(fileemp_name)
species_label = obsZ_emp[,1]

# sort=1
# if(sort==0){
#   fileunsort_name = 'unsort.csv'
#   obsZ_read = read.csv(fileunsort_name)
#   species_label = phy_prune$tip.label
# }else{
#   fileunsort_name = 'sort.csv'
#   obsZ_read = read.csv(fileunsort_name)
#   fileemp_name = 'emp.csv'
#   obsZ_emp = read.csv(fileemp_name)
#   species_label = obsZ_emp[,1]
# }

obsZ_mean = obsZ_emp[,2]

d_mean_emp = data.frame(species=species_label, traitdot=obsZ_mean)


plot_tree <- ggtree(phylo_test)+geom_tiplab(size=2)

plot_dottips = plot_tree %<+% d_mean_emp+ geom_tippoint(aes(size=traitdot))




plot_sepdots = facet_plot(plot_dottips, panel="dot", data=d_mean2, geom=geom_point, aes(x=traitdot), color='firebrick')+ theme_tree2()

plot_sepboxplt <- facet_plot(plot_dottips, panel="Trait", data=d_all, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#D1B6E1"))
# geom_phylopic(image="6fe75b0b-1488-4193-8523-f240c2d59575", color="#cff0da", alpha = .1, size=Inf)


plot_sepboxplt






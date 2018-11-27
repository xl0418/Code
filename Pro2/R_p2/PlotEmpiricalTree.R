library(phytools)
library(ggtree)
library(treeio)
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

p = ggtree(phy_prune)


trait_value = sample(1:100, 32, replace = FALSE, prob = NULL)
d = data.frame(label=phy_prune$tip.label, trait=trait_value)
p = p %<+% d
x <- as.treedata(p)

ggtree(x) + geom_tiplab(align=T, offset=.005) + geom_tippoint(aes(size=trait)) + xlim(0, 100) 
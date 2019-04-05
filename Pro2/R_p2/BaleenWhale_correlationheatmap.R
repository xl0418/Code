library(ape)
library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)
os = Sys.info()['sysname']

source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'

emdata = read.nexus(emdatadir)

baleenwhale = phylo2L(emdata,error = 1e-5)
L_ext = baleenwhale$L
extantspecieslabel = baleenwhale$ESL
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
phylo_test$tip.label <- extantspecieslabel
plot(phylo_test,show.tip.label = TRUE)


hc <- as.hclust(phylo_test) #Compulsory step as as.dendrogram doesn't have a method for phylo objects.
dend <- as.dendrogram(hc)
plot(dend, horiz=TRUE)
mat <- matrix(rnorm(15*15),nrow=15, dimnames=list(extantspecieslabel, extantspecieslabel)) #Some random data to plot
ord.mat <- mat[extantspecieslabel,extantspecieslabel]
heatmap(ord.mat, Rowv=dend, Colv=dend)

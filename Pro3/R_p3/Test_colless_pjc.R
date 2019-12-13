function (tree, norm = NULL) 
{
  if (identical(tree, NULL)) {
    stop("invalid tree", "\n")
  }
  norm.yule <- function(ICN, tree) {
    leaf_nb <- nrow(tree$merge) + 1
    EICN <- leaf_nb * log(leaf_nb) + (0.57721566 - 1 - log(2)) * 
      leaf_nb
    IC <- (ICN - EICN)/(leaf_nb)
    IC
  }
  norm.pda <- function(ICN, tree) {
    leaf_nb <- nrow(tree$merge) + 1
    IC <- ICN/((leaf_nb)^(3/2))
    IC
  }
  tmp = smaller.clade.spectrum(tree)
  res = sum(abs(tmp[, 1] - 2 * tmp[, 2]))
  if (identical(norm, NULL) == TRUE) {
    return(res)
  }
  if (norm == "pda") {
    return(norm.pda(res, tree))
  }
  if (norm == "yule") {
    return(norm.yule(res, tree))
  }
  stop("Incorrect argument for 'norm'")
}


library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
library(phyloTop)
source('C:/Liang/Code/Pro3/R_p3/deltar.R', echo=TRUE)
method = 'colless'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = rev(c('H','M','L'))
scenario = NULL
sce.short.comb.vec = NULL
for(i.letter in sce.short){
  for(j.letter in sce.short){
    sce.folder = paste0('sce',i.letter,j.letter)
    scenario = c(scenario,sce.folder)
    sce.short.comb = paste0(i.letter,j.letter)
    sce.short.comb.vec = c(sce.short.comb.vec,sce.short.comb)
  }
}

jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('1','1e-2','1e-4','1e-6','1e-8','0')
tmplist <- list()
num = 1
for(tens1 in c(1)){
  scefolder = scenario[tens1]
  letter.comb = sce.short.comb.vec[tens1]
  i.true.order = 1
  print(scefolder)
  
  colless_alldf = data.frame()
  for(i in c(1,4,2,5,3)){
    for(j in c(1)){
      comb=paste0(i,j)
      
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')
      
      trees <- read.tree(multitreefile)
      plot(trees[[1]])
      xx <- as.treeshape(trees[[1]])  # convert to apTreeshape format
      tmp = smaller.clade.spectrum(xx)
      tmplist[[num]] = tmp
      col.tree = apTreeshape::colless(xx)  # calculate colless' metric
      print(col.tree)
      num = num +1
    }
  }
}
library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
library(phyloTop)
library(quantreg)

source('C:/Liang/Code/Pro3/R_p3/deltar.R', echo=TRUE)
method = 'beta'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
dir_save = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces_results/'

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

# Compute the colless values and gamma values for a given phylo tree.
foo <- function(x, metric = "colless") {
  if (metric == "colless") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    num.tips <- x$Nnode+1
    apTreeshape::colless(xx, "yule")  #*2/(num.tips-1)/(num.tips-2)  # calculate colless' metric
    #apTreeshape::colless(xx)*2/(num.tips-1)/(num.tips-2)  # calculate colless' metric
  } else if (metric == "gamma") {
    gammaStat(x)
  } else if (metric == "deltar") {
    deltar(x)
  } else if (metric == "sackin") {
    sackin.phylo(x)
  } else if (metric == "beta") {
    maxlik.betasplit(x)$max_lik
  } else stop("metric should be one of colless or gamma")
}

jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('1','1e-2','1e-4','1e-6','1e-8','0')
x.ticks.labels = c(
  ~ "1",
  ~ list(10^{-2}),
  ~ list(10^{-4}),
  ~ list(10^{-6}),
  ~ list(10^{-8}),
  ~ "0"
)
plotl_est <- list()

for(tens1 in c(1:9)){
  scefolder = scenario[tens1]
  letter.comb = sce.short.comb.vec[tens1]
  i.true.order = 1
  print(scefolder)
  
  colless_alldf = data.frame()
  for(i in c(1,4,2,5,3)){
    for(j in c(1:6)){
      comb=paste0(i,j)
      
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')
      
      trees <- read.tree(multitreefile)
      colless_df <- ldply(trees, foo, metric = method)  # calculate metric for each tree
      colless_df = cbind(colless_df,jclabel[i.true.order],plabel[j])
      colless_alldf = rbind(colless_alldf,colless_df)
      
    }
    i.true.order=i.true.order+1
    
  }
  colnames(colless_alldf) = c('colvalue','psi','phi')
  

  
  colless_alldf$psi1 <- factor(colless_alldf$psi, labels = c('psi==0','psi==0.25','psi==0.5','psi==0.75','psi==1'))
  colless_alldf$phi1 <- factor(colless_alldf$phi, labels = c('phi==1','phi==10^-2',
                                                             'phi==10^-4','phi==10^-6','phi==10^-8','phi==0'))
  
  indexfile <- paste0(dir_save,'index/Index',method,scefolder,'.Rda')
  save(colless_alldf,file = indexfile)
  
 
  
}

library(DDD)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(viridis)
library("RColorBrewer")
library(grid)
library(gridExtra)
library(ggtree)
library(ape)
library(phytools)
library(plyr)
library(apTreeshape)
library(phyloTop)
library(quantreg)
library(treebase)
library(ggsci)
source('E:/Code/Pro3/R_p3/g_legend.R', echo=TRUE)
source('E:/Code/Pro3/R_p3/deltar.R', echo=TRUE)

dir = 'E:/Googlebox/Research/Project3/replicate_sim_9sces/'
dir.result = 'E:/Googlebox/Research/Project3/replicate_sim_9sces_results/'


# Load phylogenetic tree; Input 1
tree <- read.tree('E:/Googlebox/Research/Project3/BCI_data/control.tre')
# Read the BCI data
mydatapath = 'E:/Googlebox/Research/Project3/BCI_data/'
if (!exists('bci.full7')) {
  attach(paste(mydatapath, 'bci.full7.rdata', sep = ''))
}

full_data <- bci.full7
# Extract the species mnemonic and the geographic info; Input 2
bci_alive_data <-
  subset(full_data,
         status = 'A',
         select = c('sp', 'gx', 'gy'))
sp_code <- unique(bci_alive_data[, 1])

# Get full species names
load('E:/Googlebox/Research/Project3/BCI_data/bci.spptable.Rdata')

sp_name <- NULL
genus_name <- NULL
for (i in c(1:length(sp_code))) {
  index <- which(bci.spptable$sp == sp_code[i])
  if (length(index) == 0) {
    sp_name <-
      c(sp_name, 0)
    genus_name <-
      c(genus_name, 0)
  } else {
    sp_name <-
      c(sp_name, bci.spptable$Species[index])
    genus_name <-
      c(genus_name, bci.spptable$Genus[index])
  }
}

species_name <- paste0(genus_name, '_', sp_name)

# Subsetting the species; Input 4
species_in <- subset(species_name, (species_name %in% tree$tip.label))
di_tree <- multi2di(tree)
bci_tree <- keep.tip(di_tree, species_in)
# bci_tree = multi2di(bci_tree)

# Compute the colless values and gamma values for a given phylo tree.
foo <- function(x, metric = "colless") {
  
  if (metric == "colless") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    num.tips <- x$Nnode+1
    #apTreeshape::colless(xx, "yule")  #*2/(num.tips-1)/(num.tips-2)  # calculate colless' metric
    apTreeshape::colless(xx)*2/(num.tips-1)/(num.tips-2)  # calculate colless' metric
  } else if (metric == "colless_yule") {
    xx <- as.treeshape(x)  # convert to apTreeshape format
    num.tips <- x$Nnode+1
    apTreeshape::colless(xx, "yule")
  } else if (metric == "gamma") {
    gammaStat(x)
  } else if (metric == "deltar") {
    deltar(x)
  } else if (metric == "sackin") {
    sackin.phylo(x)
  } else if (metric == "beta") {
    maxlik.betasplit(x,size.bootstrap = 1000)$max_lik
  } else stop("metric should be one of colless or gamma")
  
  
}

deltar_BCI = foo(bci_tree, metric = "deltar")
gamma_BCI = foo(bci_tree, metric = "gamma")
colless_BCI = foo(bci_tree, metric = "colless")
beta_BCI = foo(bci_tree, metric = "beta")
# ltt(bci_tree,show.tree=TRUE,log.lineages=F,log="y")

# tree<-pbtree(n=266)
# ebTree<-phytools:::ebTree
# gamma<-function(r,tree,g) 
#   (g-ltt(ebTree(tree,r),plot=FALSE)$gamma)^2
# fit<-optimize(gamma,c(-10,10),tree=tree,g=5,tol=1e-12)
# ltt(ebTree(tree,fit$minimum),show.tree=TRUE,log.lineages=F,log="y",
#     lwd=2)


# Empirical index
emp.data.file <- 'E:/Googlebox/Research/Project3/treebase/treebase.rda'
load(emp.data.file)

beta.index = c()
colless.index = c()
gamma.index = c()
deltar.index = c()

reconstruct.tree.fun <- function(tree){
  tree = multi2di(tree)
  L = DDD::phylo2L(tree)
  if(length(which(L[,4]==-1))<3){
    return(NA)
  }else{
    reconstruct.tree = DDD::L2phylo(L)
    return(reconstruct.tree)
  }
}

data.size = length(treebase)
for(treenum in c(1:data.size)){
  print(treenum)
  tree = treebase[[treenum]]
  
  reconstruct.tree = tryCatch(reconstruct.tree.fun(tree), error=function(err) NA)
  
  beta.index = c(beta.index,tryCatch(foo(reconstruct.tree,metric = 'beta'), error=function(err) NA))
  colless.index = c(colless.index,tryCatch(foo(reconstruct.tree,metric = 'colless'), error=function(err) NA))
  gamma.index = c(gamma.index,tryCatch(foo(reconstruct.tree,metric = 'gamma'), error=function(err) NA))
  deltar.index = c(deltar.index,tryCatch(foo(reconstruct.tree,metric = 'deltar'), error=function(err) NA))
  
  
}

# color 
mypal = pal_npg("nrc", alpha = 0.7)(4)


na.beta.pos = which(is.na(beta.index))
na.colless.pos = which(is.na(colless.index))
na.gamma.pos = which(is.na(gamma.index))
na.deltar.pos = which(is.na(deltar.index))

invalid.index = unique(c(na.beta.pos,na.colless.pos,na.gamma.pos,na.deltar.pos))
whole.index = c(1:data.size)

valid.index = whole.index[is.na(pmatch(whole.index,invalid.index))]

quantiles.values = c(0.05,0.5,0.95)

quantiles.treebase.beta = quantile(beta.index[valid.index],probs = quantiles.values)
quantiles.treebase.colless = quantile(colless.index[valid.index],probs = quantiles.values)
quantiles.treebase.gamma = quantile(gamma.index[valid.index],probs = quantiles.values)
quantiles.treebase.deltar = quantile(deltar.index[valid.index],probs = quantiles.values)




jclabel = c('0','0.25','0.5','0.75','1')
plabel = c('1','1e-2','1e-4','1e-6','1e-8','0')

beta.jclabel = c('0','0.5','1','0.25','0.75')

x.ticks.labels = c(
  ~ "1",
  ~ list(10^{-2}),
  ~ list(10^{-4}),
  ~ list(10^{-6}),
  ~ list(10^{-8}),
  ~ "0"
)

sce.short = c('H','M','L')
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

x_label_fontsize = 12
y_label_fontsize = 12
mu_title_fontsize = 16
mig_title_fontsize = 16
x_title_fontsize = 16
y_title_fontsize = 16
line.size = 0.7

count1 = 1
treeshape_plot = list()
group.vec = c('g1','g2','g3')
plot.combination = list()
plot.combination[[1]] = rbind(c(1,1,1),c(5,1,1),c(9,1,1))
plot.combination[[2]] = rbind(c(2,4,3),c(2,5,3),c(2,3,5))
plot.combination[[3]] = rbind(c(6,4,3),c(6,5,3),c(6,3,5))
plot.combination[[4]] = rbind(c(9,4,3),c(9,5,3),c(9,3,5))

method1 = 'deltar'
method2 = 'gamma'
method3 = 'colless'
method4 = 'beta'


for(plot.comb in plot.combination){
  deltar_df_col = data.frame()
  gamma_df_col = data.frame()
  colless_df_col = data.frame()
  beta_df_col = data.frame()
  trees.to.plot = list()
  for(row.plot.comb in c(1:nrow(plot.comb))){
    i_n = plot.comb[row.plot.comb,1]
    i = plot.comb[row.plot.comb,2]
    j = plot.comb[row.plot.comb,3]
    comb = paste0(i,j)
    scefolder = scenario[i_n]
    letter.comb = sce.short.comb.vec[i_n]
    sce = scenario[i_n]
    print(paste0('i_n = ',i_n,'; comb = ', comb))
    
    multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')
    
    i.true.order = 1
    trees <- read.tree(multitreefile)
    trees.to.plot[[row.plot.comb]] = trees[[18]]
    deltar_df <- ldply(trees, foo, metric = method1)  # calculate metric for each tree
    deltar_df = cbind(deltar_df,group.vec[row.plot.comb],method1)
    colnames(deltar_df) = c('value','group','method')
    deltar_df_col = rbind(deltar_df_col,deltar_df)
    
    gamma_df <- ldply(trees, foo, metric = method2)  # calculate metric for each tree
    gamma_df = cbind(gamma_df,group.vec[row.plot.comb],method2)
    colnames(gamma_df) = c('value','group','method')
    gamma_df_col = rbind(gamma_df_col,gamma_df)
    
    colless_df <- ldply(trees, foo, metric = method3)  # calculate metric for each tree
    colless_df = cbind(colless_df,group.vec[row.plot.comb],method3)
    colnames(colless_df) = c('value','group','method')
    colless_df_col = rbind(colless_df_col,colless_df)
    
    indexfile <- paste0(dir.result,'index/Index',method4,scefolder,'.Rda')
    load(file = indexfile)
    beta_df = colless_alldf[colless_alldf$psi==beta.jclabel[i] & colless_alldf$phi==plabel[j] ,]$colvalue
    beta_df = cbind(as.data.frame(beta_df),group.vec[row.plot.comb],method4)
    colnames(beta_df) = c('value','group','method')
    beta_df_col = rbind(beta_df_col,beta_df)    
    
    i.true.order=i.true.order+1
      
    
  }

  # Plot trees
  treeshape_plot[[count1]] <- grid.arrange(ggtree(trees.to.plot[[1]],layout = "circular"),
                          ggtree(trees.to.plot[[2]],layout = "circular"),
                          ggtree(trees.to.plot[[3]],layout = "circular"),
                          ncol = 3,widths = c(1,1,1))
            
  count1 = count1+1
  
  # Plot deltar
  value.min <- -1.2
  value.max <- 1
  treeshape_plot[[count1]] <- ggplot(deltar_df_col, aes(x=group,y=value)) +  # plot histogram of distribution of values
    geom_boxplot(color="black",aes(fill=group), position=position_dodge(.9))+ theme_bw() +
    scale_fill_manual(values=rep("grey", 3))+
    labs(x="",y="") +
    theme(legend.position = "none",axis.title.x=element_blank(),
          axis.title.y=element_blank(),axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line()
    ) +
    geom_hline(yintercept=quantiles.treebase.deltar[1],
              linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.deltar[2],
               linetype="solid", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.deltar[3],
               linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=deltar_BCI,
               linetype="twodash", color = "red",size = line.size+0.5)+
    scale_x_discrete(labels = c())+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
          axis.text.x=element_text(size = x_label_fontsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(value.min,value.max),breaks = seq(-1,1,1),
                       labels = c(-1,0,1))
  
  count1 = count1+1
  
  # Plot gamma
  value.min <- -15
  value.max <- 6
  treeshape_plot[[count1]] <-ggplot(gamma_df_col, aes(x=group,y=value)) +  # plot histogram of distribution of values
    geom_boxplot(color = 'black',aes(fill=group), position=position_dodge(.9))+ theme_bw() +
    scale_fill_manual(values=rep("grey", 3))+
    labs(x="",y="") +
    theme(legend.position = "none",axis.title.x=element_blank(),axis.text.x=element_blank(),
          axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line()
    ) +geom_hline(yintercept=quantiles.treebase.gamma[1],
                                              linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.gamma[2],
               linetype="solid", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.gamma[3],
               linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=gamma_BCI,
               linetype="twodash", color = "red",size = line.size+0.5)+
    scale_x_discrete(labels = c())+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
          axis.text.x=element_text(size = x_label_fontsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(value.min,value.max),breaks = seq(-15,5,5))
  
  
  
  count1 = count1+1
  
  # Plot colless
  value.min <- 0
  value.max <- 1
  treeshape_plot[[count1]] <-ggplot(colless_df_col, aes(x=group,y=value)) +  # plot histogram of distribution of values
    geom_boxplot(color='black',aes(fill=group), position=position_dodge(.9))+ theme_bw() +
    scale_fill_manual(values=rep("grey", 3))+
    labs(x="",y="") +
    theme(legend.position = "none",axis.title.x=element_blank(),axis.text.x=element_blank(),
          axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line()
    ) +geom_hline(yintercept=quantiles.treebase.colless[1],
                                              linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.colless[2],
               linetype="solid", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.colless[3],
               linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=colless_BCI,
               linetype="twodash", color = "red",size = line.size+0.5)+
    scale_x_discrete(labels = c())+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
          axis.text.x=element_text(size = x_label_fontsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(value.min,value.max),breaks = seq(0,1,0.5),labels = c('0','0.5','1'))
  
  
  count1 = count1+1
  
  # Plot beta
  value.min <- -2
  value.max <- 0
  treeshape_plot[[count1]] <-ggplot(beta_df_col, aes(x=group,y=value)) +  # plot histogram of distribution of values
    geom_boxplot(color = 'black',aes(fill=group), position=position_dodge(.9))+ theme_bw() +
    scale_fill_manual(values=rep("grey", 3))+
    labs(x="",y="") +
    theme(legend.position = "none",axis.text.x=element_blank(),
          axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_line(),axis.line.x=element_line(),
          axis.ticks.y=element_line(),axis.line.y=element_line()
    ) +geom_hline(yintercept=quantiles.treebase.beta[1],
                                              linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.beta[2],
               linetype="solid", color = "black",size = line.size)+
    geom_hline(yintercept=quantiles.treebase.beta[3],
               linetype="dashed", color = "black",size = line.size)+
    geom_hline(yintercept=beta_BCI,
               linetype="twodash", color = "red",size = line.size+0.5)+
    theme(axis.text.y=element_text(angle=90,size = y_label_fontsize),
          axis.text.x=element_text(size = x_label_fontsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(limits = c(value.min,value.max),breaks = seq(-2,0,1),labels = c('-2','-1','0'))
  
  
  
  if(count1 == 5){
    treeshape_plot[[count1]] = treeshape_plot[[count1]]  +
      scale_x_discrete(name = 'Dispersal',labels = c('High','Intermediate','Low'))
      
  }else if(count1 %in% c(10,15,20))
    treeshape_plot[[count1]] = treeshape_plot[[count1]]  +
    scale_x_discrete(name = expression(psi),labels = c('0.25','0.75','1'))

  count1 = count1+1
}


m = matrix(1:20,ncol = 4)


col_labels = c('Neutral','SPJC high dispersal \n& intermediate interaction distance',
               'SPJC intermediate dispersal \n& low interaction distance','SPJC low dispersal \n& interaction distance')
row_labels = c('Tree', 'SAR','SAD','Colless','LTT')

phi1 <- textGrob(row_labels[1], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi2 <- textGrob(expression(paste(Delta,r)), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi3 <- textGrob(expression(gamma), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi4 <- textGrob(expression(I[m]), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
phi5 <- textGrob(expression(beta), gp=gpar(fontsize=mu_title_fontsize, fontface=3L))

psi1 <- textGrob(col_labels[1], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi2 <- textGrob(col_labels[2], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi3 <- textGrob(col_labels[3], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))
psi4 <- textGrob(col_labels[4], gp=gpar(fontsize=mu_title_fontsize, fontface=3L))


row_titles <- arrangeGrob(phi1,phi2,phi3,phi4,phi5,ncol = 1)
column_titles <- arrangeGrob(psi1,psi2,psi3,psi4,ncol = 4)


# label = textGrob("Number of lineages",gp=gpar(fontsize=y_title_fontsize), rot = 90)

g_ltt4 = arrangeGrob(grobs = treeshape_plot, layout_matrix = m)
g_ltt5 = textGrob("area", gp=gpar(fontsize=x_title_fontsize),rot = 0)
g_ltt1 = textGrob("")
g_a = textGrob("(a)")
g_b = textGrob("(b)")

ltt.sce <- grid.arrange(column_titles,g_ltt1,g_ltt4,row_titles,ncol = 2,widths = c(20,1),heights = c(2,25))

dir_save <- 'E:/Googlebox/Research/Project3/replicate_sim_9sces_results/'
savefilename <- paste0(dir_save,'mixing_treeshape.png')
ggsave(savefilename,ltt.sce,width = 15,height = 10)





library(DDD)
dir = 'C:/Liang/Googlebox/Research/Project2/treesim_newexp/example'

for(i in c(15:22)){
  filename = paste0(dir,i,'/Treedata.Rdata')
  pngname = paste0(dir,i,'/testtree',i,'.png')
  load(filename)
  phy = DDD::L2phylo(result$L,dropextinct = T)
  png(filename=pngname)
  plot(phy)
  dev.off()
}
library(DDD)
dir = 'C:/Liang/Googlebox/Research/Project2/treesim_newexp/example'

for(i in c(1:22)){
  filename = paste0(dir,i,'/Treedata.Rdata')
  pngname = paste0(dir,i,'/tree.png')
  load(filename)
  phy = DDD::L2phylo(result$L)
  png(filename=pngname)
  plot(phy)
  dev.off()
}
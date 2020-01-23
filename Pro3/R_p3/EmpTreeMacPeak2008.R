library("readxl")
library(ape)
data.wd = 'c:/Liang/Googlebox/Research/Project3/treedata_McPeak2008/'
data.dir = 'c:/Liang/Googlebox/Research/Project3/treedata_McPeak2008/treedata.xls'
my_data <- read_excel(data.dir)
data.length = c(1:246)
computed.gamma = as.matrix(my_data[data.length,17])

multi.phylo = my_data[data.length,23]

phylomatrix = as.matrix(multi.phylo)
for(count in data.length){
  cat(phylomatrix[count], file = paste0(data.wd,"/tree",count,".tre"), sep = "\n")
}
tree.owls <- read.tree(paste0(data.wd,"/tree",count,".tre"))

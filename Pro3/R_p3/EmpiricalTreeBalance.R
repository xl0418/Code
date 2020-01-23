library(treebase)
data.file = "c:/Liang/Googlebox/Research/Project3/treebase/treebase1000_2.rda"
cache_treebase(file = data.file,
               pause1 = 3, pause2 = 3, attempts = 2, max_trees = 1000,
               only_metadata = FALSE, save = TRUE)

Delphinus = search_treebase('Delphinus', by='taxon')

# df treebase
load('c:/Liang/Googlebox/Research/Project3/treebase/treebase100.rda')

class(treebase[[1]])
treebase[[1]]


data(treebase)
treebase

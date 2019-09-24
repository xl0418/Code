library(DDD)
dir <- 'C:/Liang/jc_standalone/x64/Release/'
filename <- paste0(dir,'Dneu.csv')
dmatrix = read.csv(filename,header = FALSE)
dmatrix = as.matrix(dmatrix)
lowertri.dmatrix = dmatrix[lower.tri(dmatrix)]
hist(lowertri.dmatrix)
min(lowertri.dmatrix)



dir <- 'C:\\Liang\\Googlebox\\Research\\Project3\\replicate_neutral\\sce1\\1e+09\\spatialpara1e+09NE11\\'
eventtablefile <- paste0(dir,'NEevent11rep100.csv')
eventtable <- read.csv(eventtablefile,header = FALSE)
turnover <- 1e9
colnames(eventtable) = c('T','ns','x','y','sp','ancestor')

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)
L.neutral <- event2L(eventtable,turnover)
phy = DDD::L2phylo(L.neutral,dropextinct = T)
plot(phy,show.tip.label = FALSE)

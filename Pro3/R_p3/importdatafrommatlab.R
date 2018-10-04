library(DDD)
dir = 'C:/Liang/Code/Pro3/10_9simdata/'

source(paste0(dir,'events00.Rdata'))
source(paste0(dir,'turnover00.Rdata'))
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)

events = as.data.frame(t(events))
colnames(events) = c('T','ns','x','y','sp','ancestor')
turnover
result = list(turnover = turnover, events = events)
L = event2L(result, mode = 'Matlab')
phy = DDD::L2phylo(L,dropextinct = F)
plot(phy)


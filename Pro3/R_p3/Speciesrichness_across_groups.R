library(DDD)
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)

dir = 'C:/Liang/Code/Pro3/10_9simdata/'

richness = matrix(0,5,5)
for(i in c(0:4)){
  for(j in c(0:4)){
    comb=paste0(i,j)
    eventname = paste0(dir,'events',comb,'.Rdata')
    source(eventname)
    source(turnovername)
    
    events = as.data.frame(t(events))
    colnames(events) = c('T','ns','x','y','sp','ancestor')
    turnover
    result = list(turnover = turnover, events = events)
    L = event2L(result, mode = 'Matlab')
    no_species = length(which(L[,4]==-1))
    richness[i+1,j+1] = no_species
  }
}




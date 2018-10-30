library(DDD)
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)

dir = 'C:/Liang/Code/Pro3/10_9simdata/'
comb='00'
eventname = paste0(dir,'events',comb,'.Rdata')
turnovername = paste0(dir,'turnover',comb,'.Rdata')
Dname = paste0(dir,'D',comb,'.Rdata')
Rname =  paste0(dir,'R',comb,'.Rdata')
source(eventname)
source(turnovername)
source(Dname)
source(Rname)
events = as.data.frame(t(events))
colnames(events) = c('T','ns','x','y','sp','ancestor')
turnover
result = list(turnover = turnover, events = events)
L = event2L(result, mode = 'Matlab')
if(is.vector(L[which(L[,4]==-1),])){
  print('Only one species exist!')  
}else{
  phy = DDD::L2phylo(L,dropextinct = T)
  plot(phy)
}





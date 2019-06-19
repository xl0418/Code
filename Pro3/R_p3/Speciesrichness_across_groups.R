library(DDD)
source('C:/Liang/Code/Pro3/R_p3/event2L.R', echo=TRUE)

dir = 'C:/Liang/Googlebox/Research/Project3/1e+06/'

richness = matrix(0,6,6)
for(i in c(1:6)){
  for(j in c(1:6)){
    comb=paste0(i,j)
    eventname = paste0(dir,'Mevent',comb,'.csv')
    turnovername = paste0(dir,'Mturnover',comb,'.Rdata')
    events = read.csv(eventname,header = FALSE)
    # source(turnovername)
    colnames(events) = c('T','ns','x','y','sp','ancestor')
    # L = event2L(events, turnover)
    no_species = tail(events$ns,n=1)
    richness[i,j] = no_species
  }
}




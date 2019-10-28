library(DDD)
library(ape)
library(parallel)
n.cores <- detectCores()

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)

neu.eventtablefile <- 'C:/Liang/jc_standalone/x64/Release/NEevent11rep100.csv'

neu.eventtable <- read.csv(neu.eventtablefile,header = FALSE)
turnover <- 1e9
L <- event2L(neu.eventtable, turnover=1e9)

ev2l <- function(eventtable){
  L <- event2L(eventtable, turnover=1e9+1e7)
  phy = DDD::L2phylo(L,dropextinct = T)
  return(phy)
}
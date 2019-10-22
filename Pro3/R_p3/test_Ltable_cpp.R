library(DDD)
library(ape)
library(parallel)
n.cores <- detectCores()

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)

dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_LM/'

neu.eventtablefile <- paste0(dir,'neutraldata/neuevent.csv')
neu.eventtable <- read.csv(neu.eventtablefile,header = FALSE)
turnover <- 1e9+1e7

ev2l <- function(eventtable){
  L <- event2L(eventtable, turnover=1e9+1e7)
  phy = DDD::L2phylo(L,dropextinct = T)
  return(phy)
}
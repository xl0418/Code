
# library(R.matlab)
# eventfile = "C:/Users/xl041/Desktop/Standalone/x64/Release/events.Rdata"
source('C:/Users/xl041/Desktop/Standalone/x64/Release/events.Rdata') 

source('C:/Users/xl041/Desktop/Standalone/x64/Release/turnover.Rdata') 

events = as.data.frame(events)
colnames(events) = c('T','ns','x','y','sp','ancestor')
turnover
result = list(turnover = turnover, events = events)
L = event2L(result, mode = 'Matlab')
phy = DDD::L2phylo(L,dropextinct = F)
plot(phy)

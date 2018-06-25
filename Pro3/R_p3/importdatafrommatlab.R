
source('C:/Liang/Code/Pro3/C++/Standalone/x64/Release/events2.Rdata')
source('C:/Liang/Code/Pro3/C++/Standalone/x64/Release/turnover1.Rdata')

events = as.data.frame(t(events))
colnames(events) = c('T','ns','x','y','sp','ancestor')
turnover
result = list(turnover = turnover, events = events)
L = event2L(result, mode = 'Matlab')
phy = DDD::L2phylo(L,dropextinct = T)
plot(phy)


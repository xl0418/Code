library(DDD)
source('C:/Liang/Googlebox/Research/Project3/Project3_R/event2L.R', echo=TRUE)
replicate = 5
for (phi in 1:replicate){
 
  for (psi in 1:replicate){
    file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",phi,"phi",psi,"psi","sim100.Rdata",sep = "")
    # file = paste("/home/p274981/cluster/",n,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    load(file = file)
    L = event2L(result = result)
    phy = DDD::L2phylo(L,dropextinct = F)
    
    # plot(phy)
    #in case the graph is too big.
    pdf(paste("C:/Liang/Googlebox/Research/Project3/simdata/",phi,"phi",psi,"psi","sim100.pdf",sep = ""))
    plot(phy)
    dev.off()
    }
}




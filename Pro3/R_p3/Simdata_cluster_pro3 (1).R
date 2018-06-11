library(POJCE)

# simulate data sets for m loc
phi = c(0,0.001,0.01,0.1,1)
psi = phi
sA = 0.5
sD = 0.5
ticks = 20000
L_size = 25
v = 0.001
# source("/home/p274981/MMM_sim_cluster.R")
  # mainDir = "/home/p274981/Proj2/"
  mainDir = "C:/Liang/Googlebox/Research/Project3/"
  subDir = "simdata"
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
 
  
replicate = length(phi)
  for (i in 1:replicate){
    phisetup = phi[i]
    for (j in 1:replicate){
      psisetup = psi[j]
    for(seed in 1:100){
  
 
    result = pojce(phi=phisetup, psi=psisetup, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
    file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",i,"phi",j,"psi","sim",seed,".Rdata",sep = "")
    # file = paste("/home/p274981/cluster/",n,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    
    # setwd(file.path(mainDir, subDir))
    # file = paste("/Users/mac/Dropbox/R/cluster/Modeltestgroupmu/out",i,"sim.Rdata",sep = "")
    save(result,file = file)
    # print(seed)
    }
    }
  }



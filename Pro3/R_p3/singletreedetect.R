phi_group = c(1)
psi_group = c(1:5)
for(phiindex in phi_group){
  for(psiindex in psi_group){
    for (seed in 1:100) {
    
      file = paste("C:/Liang/Googlebox/Research/Project3/simdata/",phiindex,"phi",psiindex,"psi","sim",seed,".Rdata",sep = "")
      load(file = file)
      if(result$events$ns[length(result$events$ns)] == 1){
        print(c(phiindex,psiindex,seed))
      }
    }
  }
}
      
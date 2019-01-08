phi_vec = c(0,0.001,0.01,0.1,0.5,1)
psi_vec = c(0,0.001,0.01,0.1,0.5,1)
L = 333
v=0.0001

scenario = 2
if(scenario == 0){
  sA=0.1
  sB=0.1
  sce = "L"
  scefile = "low"
}else if(scenario == 1){
  sA=10
  sB=10
  sce = "M"
  scefile = "med"
  
}else if(scenario == 2){
  sA=1000
  sB=1000
  sce = "H"
  scefile = "high"
  
}else{
  print("Pls provide scenario in (0,1,2)!")
}

ticks=1000000000
log=1e8
count = 0
append = TRUE
for(i in c(1:6)){
  for(j in c(1:6)){
    count = count +1
    # str = paste0('L=',L,' v=',v,' phi=',phi_vec[i], ' psi=',psi_vec[j],
    #              ' sA=',sA, ' sB=',sB, ' ticks=',ticks, ' log=',log,' continue=neutral.m', ' file=Lphi',i,'psi',j
    #              ,'.m')
    str = sprintf('L=%i v=%.4f phi=%.4f psi=%.4f sA=%.1f sB=%.1f ticks=%i log=%i continue=neutral.m file=%sphi%ipsi%i.m',
                  L,v,phi_vec[i],psi_vec[j],sA,sB,ticks, log,sce,i,j)
    filename<-paste0("C:/Liang/Googlebox/Research/Project3/simdata_190103/spatialpara",scefile,".txt")
    if(count == 1){
      # write(str, file=filename,append=FALSE)
      cat(str,'\n', file=filename, append=FALSE, sep='')
    }else{
      cat(str ,'\n', file=filename, append=append, sep='')
      print(str)
    }
  }
}

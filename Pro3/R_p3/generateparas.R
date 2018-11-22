phi_vec = c(0,0.001,0.01,0.1,0.5,1)
psi_vec = c(0,0.001,0.01,0.1,0.5,1)
L = 333
v=0.0001
sA=50
sB=50
ticks=1e8
log=1000000
count = 0
append = TRUE
for(i in c(1:6)){
  for(j in c(1:6)){
    count = count +1
    str = paste0('L=',L,' v=',v,' phi=',phi_vec[i], ' psi=',psi_vec[j],
                 ' sA=',sA, ' sB=',sB, ' ticks=',ticks, ' log=',log, ' file=phi',i,'psi',j
                 ,'.m')
    filename<-"C:/Liang/Googlebox/Research/Project3/simdata_181122/spatialpara.txt"
    if(count == 1){
      # write(str, file=filename,append=FALSE)
      cat(str,'\n', file=filename, append=FALSE, sep='')
    }else{
      cat(str ,'\n', file=filename, append=append, sep='')
      print(str)
    }
  }
}

psi_vec = c(0,0.5,1,0.25,0.75)
sig_phi = c(sqrt(2)/2,sqrt(2)*1e2/2,sqrt(2)*1e4/2,sqrt(2)*1e6/2,sqrt(2)*1e8/2,-1)
L = 333
v=0.0001
disp_vec = c(0.1,1,10)
spar_vec = c(0.1,1,10)
name.spar=c('L','M','H')
batch.size=10
dir = "C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces"
dir.folder = paste0(dir,"/1e+07")
dir.create(file.path(dir.folder), showWarnings = FALSE)
ticks=10000000
log=1e8
for(dis.ind in c(1:3)){
  for(spar.ind in c(1:3)){
  # scenario = 0
    s_spar=spar_vec[spar.ind]
    s_disp=disp_vec[dis.ind]
    sce = paste0(name.spar[dis.ind],name.spar[spar.ind])
    dir.scefolder = paste0(dir,"/sce",spar.ind,dis.ind)
    dir.folder = paste0(dir.scefolder,"/1e+07")
    dir.create(file.path(dir.scefolder), showWarnings = FALSE)
    dir.create(file.path(dir.folder), showWarnings = FALSE)
    
    dir.simfolder = paste0(dir.scefolder,"/sim_scripts")
    dir.create(file.path(dir.simfolder), showWarnings = FALSE)
    append = TRUE
    for(i in c(1:5)){
      for(j in c(1:6)){
        subDir = paste0("spatialpara",formatC(ticks),sce,i,j)
        dir.create(file.path(dir.folder, subDir), showWarnings = FALSE)
        setwd(file.path(dir.folder, subDir))
        count = 0
        for(rep.sim in c(1:100)){
          
          
          str = sprintf('L=%i v=%.4f Psi=%.4f s_phi=%.4f s_spar=%.1f s_disp=%.1f ticks=%i seed=%i file=%s/%s/%spsi%is_phi%irep%i.m',
                        L,v,psi_vec[i],sig_phi[j],s_spar,s_disp,ticks,rep.sim,formatC(ticks),subDir,sce,i,j,rep.sim)
          if(count%%batch.size == 0){
            no.batch = count%/%batch.size+1
            filename<-paste0(dir.simfolder,"/spatialpara",formatC(ticks),sce,i,j,"batch",no.batch,".txt")
            
            # write(str, file=filename,append=FALSE)
            cat(str,'\n', file=filename, append=FALSE, sep='')
          }else{
            cat(str ,'\n', file=filename, append=append, sep='')
            print(str)
          }
          count = count +1
          
        }
      }
    }
  
  }
}


psi_vec = c(0,0.5,1,0.25,0.75)
sig_phi = c(0,1e2,1e4,1e6,1e8,-1)
L = 333
v=0.0001
batch.size=10
dir = "C:/Liang/Googlebox/Research/Project3/replicate_sim_smallbatch"
ticks=10000000
log=1e8
for(scenario in c(0:2)){
  # scenario = 0
  if(scenario == 0){
    s_spar=0.1
    s_disp=0.1
    sce = "L"
    scefile = "low"
  }else if(scenario == 1){
    s_spar=10
    s_disp=10
    sce = "M"
    scefile = "med"
    
  }else if(scenario == 2){
    s_spar=1000
    s_disp=1000
    sce = "H"
    scefile = "high"
    
  }else{
    print("Pls provide scenario in (0,1,2)!")
  }
  
  append = TRUE
  for(i in c(4:length(psi_vec))){
    for(j in c(1:length(sig_phi))){
      subDir = paste0("spatialpara",formatC(ticks),sce,i,j)
      dir.create(file.path(dir, subDir), showWarnings = FALSE)
      setwd(file.path(dir, subDir))
      count = 0
      for(rep.sim in c(1:100)){
        
        
        str = sprintf('L=%i v=%.4f Psi=%.4f s_phi=%.4f s_spar=%.1f s_disp=%.1f ticks=%i seed=%i continue=neutral.m file=%s/%s/%spsi%is_phi%irep%i.m',
                      L,v,psi_vec[i],sig_phi[j],s_spar,s_disp,ticks,rep.sim,formatC(ticks),subDir,sce,i,j,rep.sim)
        if(count%%batch.size == 0){
          no.batch = count%/%batch.size+1
          filename<-paste0(dir,"/spatialpara",formatC(ticks),sce,i,j,"batch",no.batch,".txt")
          
          # write(str, file=filename,append=FALSE)
          cat(str,'\n', file=filename, append=FALSE, sep='')
        }else{
          cat(str ,'\n', file=filename, append=append, sep='')
          print(str)
        }
        count = count +1
        
      }
      
      #       # generate bash files
      #       bashfilename = paste0(dir,"/spatialsim",scefile,formatC(ticks),sce,i,j,".sh")
      #       bashsetting = sprintf('#!/bin/bash
      # #SBATCH --time=0:10:00
      # #SBATCH --partition=gelifes
      # #SBATCH --ntasks=1
      # #SBATCH --nodes=1
      # #SBATCH --cpus-per-task=5
      # #SBATCH --mem=12GB 
      # #SBATCH --job-name=%s%i%isim
      # #SBATCH --output=PJC%s%i%i.log
      # #SBATCH --mail-type=FAIL,TIME_LIMIT 
      # #SBATCH --mail-user=xl0418@gmail.com 
      # ./jc batch=%s.txt
      #   ',sce,i,j,sce,i,j,subDir
      #       )
      #       cat(bashsetting,file=bashfilename)
    }
  }
  
}


psi_vec = c(0)
sig_phi = c(0)
L = 333
v=0.0001
batch.size=10
dir = "C:/Liang/Googlebox/Research/Project3/replicate_neutral"


ticks=1000000000
log=1e9
for(scenario in c(1)){
  dir.scefolder = paste0(dir,"/sce",scenario)
  dir.folder = paste0(dir.scefolder,"/1e+07")
  dir.create(file.path(dir.scefolder), showWarnings = FALSE)
  dir.create(file.path(dir.folder), showWarnings = FALSE)
  
  dir.simfolder = paste0(dir.scefolder,"/sim_scripts")
  dir.create(file.path(dir.simfolder), showWarnings = FALSE)
  # scenario = 1
  if(scenario == 1){
    s_spar=-1
    s_disp=-1
    sce = "NE"
  }else if(scenario == 2){
    s_spar=0.1
    s_disp=10
    sce = "LH"
  }else if(scenario == 3){
    s_spar=1
    s_disp=10
    sce = "MH"
  }else if(scenario == 4){
    s_spar=1
    s_disp=0.1
    sce = "ML"
  }else if(scenario == 5){
    s_spar=10
    s_disp=0.1
    sce = "HL"
  }else if(scenario == 6){
    s_spar=10
    s_disp=1
    sce = "HM"
    
  }else{
    print("Pls provide scenario in (1:6)!")
  }
  
  append = TRUE
  for(i in c(1)){
    for(j in c(1)){
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


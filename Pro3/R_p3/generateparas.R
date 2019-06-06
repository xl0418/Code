psi_vec = c(0,0.2,0.4,0.6,0.8,1)
sig_phi = c(0,1e2,1e4,1e6,1e8,-1)
L = 333
v=0.0001

ticks=100000000
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

  count = 0
  append = TRUE
  for(i in c(1:length(psi_vec))){
    for(j in c(1:length(sig_phi))){
      count = count +1
      # str = paste0('L=',L,' v=',v,' Psi=',psi_vec[i], ' s_phi=',sig_phi[j],
      #              ' s_spar=',sA, ' s_disp=',sB, ' ticks=',ticks,' continue=neutral.m', ' file=Lpsi',i,'s_phi',j
      #              ,'.m')
      str = sprintf('L=%i v=%.4f Psi=%.4f s_phi=%.4f s_spar=%.1f s_disp=%.1f ticks=%i continue=neutral.m file=%s/%spsi%is_phi%i.m',
                    L,v,psi_vec[i],sig_phi[j],s_spar,s_disp,ticks,formatC(ticks),sce,i,j)
      filename<-paste0("C:/Liang/Googlebox/Research/Project3/simdata_1e9/spatialpara",scefile,formatC(ticks),".txt")
      if(count == 1){
        # write(str, file=filename,append=FALSE)
        cat(str,'\n', file=filename, append=FALSE, sep='')
      }else{
        cat(str ,'\n', file=filename, append=append, sep='')
        print(str)
      }
    }
  }
  
  # generate bash files
  bashfilename = paste0("C:/Liang/Googlebox/Research/Project3/simdata_1e9/spatialsim",scefile,formatC(ticks),".sh")
  bashsetting = sprintf('#!/bin/bash
#SBATCH --time=9-23:59:00
#SBATCH --partition=gelifes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=12GB 
#SBATCH --job-name=%s%ssim
#SBATCH --mail-type=FAIL,TIME_LIMIT 
#SBATCH --mail-user=xl0418@gmail.com 
./jc batch=spatialparamed%s.txt
  ',formatC(ticks),scefile,formatC(ticks)
  )
  cat(bashsetting,file=bashfilename)
}


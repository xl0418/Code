phi_vec = c(0,0.1,0.3,1,3,10)
psi_vec = c(0,0.2,0.4,0.6,0.8,1)
L = 333
v=0.0001

ticks=100000000
log=1e8
scefile = 'imp'
sce = 'HR'
count = 0
append = TRUE
for(i in c(1:6)){
  for(j in c(1:6)){
    count = count +1
    # str = paste0('L=',L,' v=',v,' phi=',phi_vec[i], ' psi=',psi_vec[j],
    #              ' sA=',sA, ' sB=',sB, ' ticks=',ticks, ' log=',log,' continue=neutral.m', ' file=Lphi',i,'psi',j
    #              ,'.m')
    str = sprintf('-i L=%i v=%.4f phi=%.4f psi=%.4f ticks=%i log=%i continue=neutral.m file=%s/%sphi%ipsi%i.m',
                  L,v,phi_vec[i],psi_vec[j],ticks, log,formatC(ticks),sce,i,j)
    filename<-paste0("C:/Liang/Googlebox/Research/Project3/simdata_",formatC(ticks),"implicit/spatialpara",scefile,formatC(ticks),".txt")
    if(count == 1){
      # write(str, file=filename,append=FALSE)
      cat(str,'\n', file=filename, append=FALSE, sep='')
    }else{
      cat(str ,'\n', file=filename, append=append, sep='')
      print(str)
    }
  }
}

# generate bash files for gelifes partition
bashfilename = paste0("C:/Liang/Googlebox/Research/Project3/simdata_",formatC(ticks),"implicit/spatialsim",scefile,formatC(ticks),".sh")
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
./jc batch=spatialparaimp%s.txt
',formatC(ticks),scefile,formatC(ticks)
)
cat(bashsetting,file=bashfilename)


# generate bash files for regular partition
bashfilename = paste0("C:/Liang/Googlebox/Research/Project3/simdata_",formatC(ticks),"implicit/reguspatialsim",scefile,formatC(ticks),".sh")
bashsetting = sprintf('#!/bin/bash
#SBATCH --time=9-23:59:00
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=12GB 
#SBATCH --job-name=%s%ssim
#SBATCH --mail-type=FAIL,TIME_LIMIT 
#SBATCH --mail-user=xl0418@gmail.com 
./jc batch=spatialparaimp%s.txt
',formatC(ticks),scefile,formatC(ticks)
)
cat(bashsetting,file=bashfilename)



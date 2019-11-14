
event2L = function(events,turnover){
  colnames(events) = c('T','ns','x','y','sp','ancestor')
  ticks = turnover
  L_raw = cbind(events$T,events$ancestor+1,events$sp+1)
  
  L_raw[which(L_raw[,2]==0),2] = -1
  L_raw = rbind(c(0,0,1),c(0,1,2),L_raw)
  L_raw = cbind (L_raw,-1)
  L_raw[,1] = ticks - L_raw[,1]
  
  ext_rownum = which(L_raw[,2] ==-1)
  ext.spec.index = L_raw[ext_rownum,3]
  ext.time.vec = L_raw[ext_rownum,1]
  if(length(ext_rownum)>0){
    
    for(indicator in c(1:length(ext_rownum))){

      corrector = 2*(indicator-1)
      ext.spec = ext.spec.index[indicator]
      ext.row = ext_rownum[indicator]-corrector
      ext.spec.time = ext.time.vec[indicator]
      dead.spec.row.candidate = which(L_raw[1:(ext.row-1),3]==ext.spec)

      if(length(dead.spec.row.candidate)!=1){
        print("More species are dieing than expected (1).")
        break
      }
      possible.ext.mother = which(L_raw[1:(ext.row-1),2]==ext.spec)
      if(length(possible.ext.mother)>0){
        mother.mother = L_raw[ext.spec,2]
        L_raw[possible.ext.mother,2] = mother.mother
        
      }
      
      change.spec = which(L_raw[1:(ext.row-1),3]>ext.spec)
      L_raw[change.spec,3] = L_raw[change.spec,3]-1
      change.mot.spec = which(L_raw[1:(ext.row-1),2]>ext.spec)
      L_raw[change.mot.spec,2] = L_raw[change.mot.spec,2]-1
      if(length(which(L_raw[change.mot.spec,1]==176726339))){
        print("The target branching time is changing")
      }
      
      
      # remove the extinction information row
      L_raw = L_raw[-ext.row,]
      # remove the extinct species
      L_raw = L_raw[-dead.spec.row.candidate,]
      if(length(which(L_raw[,2]==L_raw[,3]))>=1){
        print("Born from itself!")
      }
      
      
      # 
      # progress.i = progress.i+1
      # i = progress.i/length(ext_rownum)*100
      # progress(i, progress.bar = TRUE)
      # Sys.sleep(0.01)
      # if (i > 100) cat("Done!\n")
      }
    L_ini <- L_raw[which(L_raw[,2]>=0),]
    L <- L_ini[which(L_ini[,4]==-1),]
  }

  return(L)
}

event2L = function(events,turnover){
  colnames(events) = c('T','ns','x','y','sp','ancestor')
  ticks = turnover
  L_raw = cbind(events$T,events$ancestor+1,events$sp+1)
  
  L_raw[which(L_raw[,2]==0),2] = -1
  L_raw = rbind(c(0,0,1),c(0,1,2),L_raw)
  L_raw = cbind (L_raw,-1)
  L_raw[,1] = ticks - L_raw[,1]
  
  ext_rownum = which(L_raw[,2] ==-1)
  if(length(ext_rownum)>0){
    
    for(ext.row in ext_rownum){
      L_raw[ext.row,4] = 0
      ext.spec = L_raw[ext.row,3]
      ext.spec.time = L_raw[ext.row,1]

      dead.spec.row.candidate = which(L_raw[1:(ext.row-1),3]==ext.spec)
      dead.spec.row = dead.spec.row.candidate[which(L_raw[dead.spec.row.candidate,4]==-1)]
      if(length(dead.spec.row)!=1){
        print("More species are dieing than expected (1).")
        break
      }

      L_raw[dead.spec.row,4]=L_raw[ext.row,1]
      change.spec = which(L_raw[1:(ext.row-1),3]>ext.spec)
      L_raw[change.spec,3] = L_raw[change.spec,3]-1
      change.mot.spec = which(L_raw[1:(ext.row-1),2]>ext.spec)
      L_raw[change.mot.spec,2] = L_raw[change.mot.spec,2]-1
      
      
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
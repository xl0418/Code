
event2L = function(result){
  ticks = result$log$sT[length(result$log$sT)]
  
  L_raw = cbind(result$events$T,result$events$ancestor+1,result$events$sp+1)
  L_raw[which(L_raw[,2]==0),2] = -1
  L_raw = rbind(c(0,0,1),c(0,1,2),L_raw)
  L_raw = cbind (L_raw,-1)
  L_raw[,1] = ticks - L_raw[,1]
  L_exi = L_raw[which(L_raw[,2] >=0),]
  
  
  ext_rownum = which(L_raw[,2] ==-1)
  if(length(ext_rownum)>0){
    L_ext2 = L_exi
    ext_vector = L_raw[ext_rownum,,drop = FALSE]
    ext_rownum = c(0,ext_rownum)
    L_block = NULL
    for (i in 1:nrow(ext_vector)){
      # print(i)
      if((ext_rownum[i]+1)<=(ext_rownum[i+1]-1)){
        L_block = rbind(L_block,L_raw[(ext_rownum[i]+1):(ext_rownum[i+1]-1),,drop = FALSE])
        row_index_ext = which(L_block[,3]%in% ext_vector[i,3] )
        row_ext =row_index_ext[length(row_index_ext)]
        L_block[which(L_block[,3]>ext_vector[i,3]),3] =  L_block[which(L_block[,3]>ext_vector[i,3]),3]-1
        L_exi[row_ext,4] = ext_vector[i,1]
        L_block[row_ext,3] = -1
     
      }
      else{
        ext_vector_updata = L_block[,3]
         row_index_ext = which(ext_vector_updata %in% ext_vector[i,3] )
        row_ext =row_index_ext[length(row_index_ext)]
        L_exi[row_ext,4] = ext_vector[i,1]
        L_block[row_ext,3] = -1        
        L_block[which(L_block[,3]>ext_vector[i,3]),3] =  L_block[which(L_block[,3]>ext_vector[i,3]),3]-1
        
      }
    }
    L_exi[,3] = c(1:nrow(L_exi))
    for(i in 1:nrow(L_exi)){
      match_vec = which(L_ext2[1:i,3] %in% L_exi[i,2])
      if(length(match_vec)>0 ){
        L_exi[i,2] = L_exi[match_vec[length(match_vec)],3]
      }
    }
  }
  return(L_exi)
}
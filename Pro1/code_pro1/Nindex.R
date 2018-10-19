Nindex = function(n){
Ntable_index = matrix(0,nrow= n, ncol =n)
for(i in 1:n){
  for(j in 1:(n-i)){
    Ntable_index[i,j+1] = choose(n-i,j)
  }
}
Ntable_index[,1] = 1
return(Ntable_index)
}
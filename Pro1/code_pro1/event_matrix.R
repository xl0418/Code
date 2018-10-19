event_matrix = function(n){
  num = 0
  for(i in 1:n)
  {
    num = num + choose(n,i)
  }
  E1 = matrix(0,nrow = n, ncol = num)
  col = 0
  for(i in 1:n){
  index = combn(n,i)
  for(j in 1:ncol(index)){
  E1[index[,j], j+col] = 1
  }
  col = col + ncol(index)
  }
  # E1 = Matrix(E1,sparse = TRUE)
  
  
 return(E1)
}
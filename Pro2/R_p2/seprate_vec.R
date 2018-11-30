seprate_vec = function(vec,length,leftend = 0,binwidth=0.1){
  sepvec = c()
  rightend = leftend + binwidth
  for(i in c(1:length)){
    leftend = rightend
    rightend = rightend+binwidth
    num = sum(vec > leftend& vec < rightend)
    sepvec = rbind(sepvec, c(num,(leftend+rightend)/2))
  }
  return(sepvec)
}

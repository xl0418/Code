deltar <- function(L){
  branchtime.first <- L[1,1]
  L.beginning <- length(which(L[,1]==branchtime.first))
  L.tips <- nrow(L)
  half.tree.length <- branchtime.first/2
  L.half <- L.tips-length(which(L[,1]<=half.tree.length))
  deltar.value <- (log(L.tips/L.half)-log(L.half/L.beginning))/(log(L.tips/L.half)+log(L.half/L.beginning))
  return(deltar.value)
}
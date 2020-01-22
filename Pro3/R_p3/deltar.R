deltar <- function(L){
  if(class(L)=='phylo')L <- phylo2L(L)
  branchtime.first <- L[1,1]
  L.beginning <- 2
  L.tips <- nrow(L)
  half.tree.length <- branchtime.first/2
  L.half <- L.tips-length(which(L[,1]<=half.tree.length))
  deltar.value <- (log(L.tips/L.half)-log(L.half/L.beginning))/(log(L.tips/L.half)+log(L.half/L.beginning))
  return(deltar.value)
}
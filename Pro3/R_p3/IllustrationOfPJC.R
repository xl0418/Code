abundance.vec = c(2,10,500,800)
phylogenetic.matrix = matrix(c(0,5,15,20,
                               5,0,15,20,
                               15,15,0,20,
                               20,20,20,0),nrow = 4)
phi = 1


phylo.spatial.effect <- function(abundance,phymatrix,phi){
  contri <- exp(-phi*phymatrix)%*%abundance / sum(abundance)
}

phi.contribution = c()
for(phi in seq(0,1,0.01)){
  phi.temp <-  phylo.spatial.effect(abundance = abundance.vec,phymatrix=phylogenetic.matrix,
                                    phi = phi)
  phi.contribution <- rbind(phi.contribution,t(phi.temp))
}

